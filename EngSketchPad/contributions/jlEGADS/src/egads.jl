#=
*
*              jlEGADS --- Julia version of EGADS API
*
*
*      Copyright 2011-2021, Massachusetts Institute of Technology
*      Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
*
=#
push!(LOAD_PATH,".")

__precompile__()

module egads

#import ExportAll

# Dependencies
include(ENV["jlEGADS"]*"/deps/sharedFile.jl")
include(EGADS_WRAP_PATH*"/"*FILE_EGADS_CTT)

struct ErrorMessage <: Exception
    msg::String
end

struct egadsERR_NAMES
    var  ::Int32
    name ::String
end

egadsERROR = []


# MY DATA TYPES
const str = Union{String,Char}
const int = Union{Int32,Int64}
const dbl = Union{Float32,Float64}

const vInt = Union{Vector{Int32}  ,Vector{Int64}}
const vDbl = Union{Vector{Float32},Vector{Float64}}
const vNum = Union{vInt,vDbl}
const vStr = Union{Vector{Char},Vector{String}}

const mInt = Union{Matrix{Int64},Matrix{Int32}}
const mDbl = Union{Matrix{Float64},Matrix{Float32}}
const mNum = Union{mInt,mDbl}
const mStr = Union{Matrix{String}, Matrix{Char}}

# Check if a variable || vector is empty
emptyVar(a) = all(x -> x ∈ [C_NULL,nothing,Nothing], vcat([a]...))

mutable struct CSys
    val   :: Vector{Float64}
    ortho :: Matrix{Float64}
    function CSys(val::vNum, ortho = nothing::Union{Nothing, mNum})
        if ortho == nothing
            ortho = [[0.0, 0.0, 0.0] [1.0, 0.0, 0.0] [0.0, 1.0, 0.0] [0.0, 0.0, 1.0]]
        end
        x = new(Float64.(val), ortho)
    end
end

"""
Free memory in pointer. Called whenever a EG_() function has freeable arguments

Parameters
-------
object: the ego Object to delete
"""
free(obj) = Base.Libc.free(obj[])


"""
Wraps an egads CONTXT to trace memory

Members
----------
ctxt:
    an egads ego

delObj:
    If True then this ego takes ownership of obj and calls EG_close during garbage collection

name:
    optional name assignation (for debugging)
"""
mutable struct Ccontext
    ctxt   :: ego
    delObj :: Bool
    name   :: String
    function Ccontext(; obj=ego()::ego, delObj::Bool = false, name = "ctxt")
        if obj == C_NULL
            ptr = Ref{ego}()
            raiseStatus(ccall((:EG_open, C_egadslib), Cint, (Ptr{ego},), ptr))
            obj = ptr[] ; delObj = true
        else
            info = unsafe_load(obj)
            info.oclass != CONTXT && throw(ErrorMessage("struct Ccontext oclass needs to be CONTXT !!"))
        end
        x = new(obj,delObj,name)
        function f(x)
            if x.delObj
                z = ccall((:EG_close, C_egadslib), Cint, (ego,), x.ctxt)
                x.ctxt = ego() ; x.delObj = false
            end
        end
        finalizer(f, x)
    end
end


"""
Wraps an egads ego to trace memory. Called whenever an ego is created through create_Cego()

Members
----------
obj:
    an egads ego that is not context: CONTEXT has its own struct: Ccontext

ctxt:
    the Context associated with the ego

refs:
    list of ego's that are used internaly in obj. Used to ensure proper deletion order

delObj:
    If True then this ego takes ownership of obj and calls EG_deleteObject during garbage collection

name:
    optional name assignation (for debugging)
"""
mutable struct Cego
    obj    :: ego
    ctxt   :: Union{Ccontext,Nothing}
    refs   :: Union{Vector{Cego}, Vector{Nothing}}
    delObj :: Bool
    name   :: String
    function Cego(obj::ego, ctxt::Union{Ccontext, Nothing}, refs::Union{Vector{Cego},Vector{Nothing}},delObj::Bool, name::String)
        x = new(obj, ctxt, refs, delObj, name)
        function f(x)
            z = x.delObj ? deleteObject!(x) : 0
            x = create_Cego()
        end
        finalizer(f,x)
        return x
    end
end

function create_Cego(obj=ego()::ego; ctxt=nothing::Union{Ccontext,Nothing},
                     refs=nothing::Union{Nothing,Vector{Cego},Cego},delObj=true::Bool, name = "genC"::String)

    emptyVar(obj) && return Cego(ego(), nothing, [nothing], false, "void")
    if ctxt == nothing
        ptr = Ref{ego}()
        raiseStatus(ccall((:EG_getContext, C_egadslib), Cint, (ego, Ptr{ego}),obj, ptr))
        cobj = ptr[]
    else
        cobj = ctxt.ctxt
    end
    return Cego(obj, Ccontext(obj = cobj,delObj = false), vcat([refs]...), delObj, name)
end

function reset_Cego!(x::Union{Cego,Nothing})
    typeof(x) == Nothing && return
    x.obj = ego()
    x.refs = [nothing]
    x.delObj = false
    x.ctxt = nothing
end


"""
Deletes an EGADS object.

Parameters
-------
object: the ego Object to delete
"""
function deleteObject!(x::Cego)
    icode = ccall((:EG_deleteObject, C_egadslib), Cint, (ego,), x.obj)
    (icode < 0) && raiseStatus(icode)
    return icode
end


"""
Returns the version information for both EGADS and OpenCASCADE.

Returns
-------

imajor: major version number

iminor: minor version number

OCCrev: OpenCASCADE version string
"""
function revision( )
    major  = Ref{Cint}(0)
    minor  = Ref{Cint}(0)
    OCCrev = Ref{Cstring}()
    ccall((:EG_revision, C_egadslib), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cstring}), major, minor, OCCrev)
    return (imajor = major[], iminor = minor[], OCCrev = unsafe_string(OCCrev[]) )
end


"""
Loads a model object from disk and puts it in the Context

Parameters
----------
bitFlag: Options (additive):

        1 Don't split closed and periodic entities

        2 Split to maintain at least C1 in BSPLINEs

        4 Don't try maintaining Units on STEP read (always millimeters) 8 Try to merge Edges and Faces (with same geometry)

        16 Load unattached Edges as WireBodies (stp/step & igs/iges)

name:
    path of file to load (with extension - case insensitive):

    igs/iges    IGES file

    stp/step    STEP file

    brep        native OpenCASCADE file

    egads       native file format with persistent Attributes, splits ignored)

Returns
-------
model: the model object
"""
function loadModel(context::Ccontext, bitFlag::int, name::str)
    ptr = Ref{ego}()
    raiseStatus((ccall((:EG_loadModel, C_egadslib), Cint, (ego, Cint, Cstring, Ptr{ego}),
                    context.ctxt, bitFlag, Base.cconvert(Cstring,name), ptr)))
    return create_Cego(ptr[] ; ctxt =  context)
end


"""
Resets the Context's owning thread ID to the thread calling this function.

Parameters
----------
context:
    the Context Object to update
"""
updateThread!(context::Ccontext) = raiseStatus(ccall((:EG_updateThread, C_egadslib), Cint, (ego,), context.ctxt))


"""
Writes the geometry to disk. Only writes analytic geometric data for anything but EGADS output.
Will not overwrite an existing file of the same name unless overwrite = true.

Parameters
----------
object:
    the Body/Model Object to write

name:

    path of file to write, type based on extension (case insensitive):

    igs/iges   IGES file

    stp/step   STEP file

    brep       a native OpenCASCADE file

    egads      a native file format with persistent Attributes and
               the ability to write EBody and Tessellation data
overwrite:

    Overwrite an existing file
"""
function saveModel!(model::Cego, name::str, overwrite::Bool = false)
    overwrite && rm(name, force = true)
    raiseStatus(ccall((:EG_saveModel, C_egadslib), Cint, (ego, Cstring), model.obj, name))
end


"""
Makes a Transformation Object from a [4][3] translation/rotation matrix.
The rotation portion [1:3,:] must be "scaled" orthonormal (orthogonal with a
single scale factor).

Parameters
----------
context:
    the ego context

mat:
    the 12 values of the translation/rotation matrix.
    In matrix form [4][3] where [:,1:3] is a reference system


Returns
-------
xform: the transformation object
"""
function makeTransform(context::Ccontext, mat::mDbl)
    mat2 = vcat(mat...)
    ptr  = Ref{ego}()
    raiseStatus(ccall((:EG_makeTransform, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
                context.ctxt, mat2, ptr))
        return create_Cego(ptr[] ; ctxt = context)
end


"""
Returns the [4][3] transformation information.

Parameters
----------
xform:
    the Transformation Object.

"""
function getTransform(xform::Cego)
    mat = zeros(12)
    raiseStatus(ccall((:EG_getTransformation, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                xform.obj, mat))
    return reshape(mat, 4, 3)
end


"""
Creates a new EGADS Object by copying and reversing the input object.
Can be Geometry (flip the parameterization) or Topology (reverse the sense).
Not for Node, Body or Model. Surfaces reverse only the u parameter.

Parameters
----------
object:
    the Object to flip

Returns
-------
newObject: the resultant flipped ego
"""
function flipObject(object::Cego)
 ptr = Ref{ego}()
 raiseStatus(ccall((:EG_flipObject, C_egadslib), Cint, (ego, Ptr{ego}), object.obj, ptr))
 return create_Cego(ptr[] ; ctxt =  object.ctxt)
end


"""
Creates a new EGADS Object by copying and optionally transforming the input object.
If the Object is a Tessellation Object, then other can be a vector of
displacements that is 3 times the number of vertices of doubles in length.

Parameters
----------
object:
    the Object to copy

other:
    the transformation or context object (an ego) --
    Nothing for a strict copy
    a displacement vector for TESSELLATION Objects only
    (number of global indices by 3 doubles in length)

Returns
-------
newObject: the resultant new ego
"""
function copyObject(object::Cego, other::Union{Cego,Nothing} = nothing)
    ptr   = Ref{ego}()
    oform = emptyVar(other) ? C_NULL : typeof(other) == Cego ? other.obj : other[1:3]
    raiseStatus(ccall((:EG_copyObject, C_egadslib), Cint, (ego, Ptr{Cvoid}, Ptr{ego}),
                 object.obj, oform, ptr))
    return create_Cego(ptr[], ctxt = object.ctxt)
end


"""
Return information about the context

Returns
-------
oclass:
    object class type

mtype:
    object sub-type

topObj:
    is the top level BODY/MODEL that owns the object or context (if top)

prev:
    is the previous object in the threaded list (Nothing at CONTEXT)

next:
    is the next object in the list (Nothing is the end of the list)
"""
function getInfo(object::Union{Cego, Ccontext})
    oclass = Ref{Cint}(0)
    mtype  = Ref{Cint}(0)
    topObj = Ref{ego}()
    prev   = Ref{ego}()
    next   = Ref{ego}()
    obj    = typeof(object) == Cego ? object.obj : object.ctxt
    raiseStatus(ccall((:EG_getInfo, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}),
                    obj, oclass, mtype, topObj, prev, next))
    ctxt = typeof(object) == Cego ? object.ctxt : object
    return (oclass=oclass[], mtype=mtype[],
            topObj=create_Cego(topObj[];ctxt = ctxt, delObj = false),
             prev=create_Cego(prev[];   ctxt = ctxt, delObj = false),
             next=create_Cego(next[];   ctxt = ctxt, delObj = false))
end


"""
Sets the EGADS verbosity level, the default is 1.

Parameters
----------
outlevel:

    0 <= outlevel <= 3

    0-silent to 3-debug
"""
setOutLevel(context::Ccontext, outLevel::int) = ccall((:EG_setOutLevel, C_egadslib), Cint, (ego, Cint), context.ctxt, Int32(outLevel))


"""
Create a stream of data serializing the objects in the Model.

Returns
-------

stream:  the byte-stream
"""
function exportModel(model::Cego)
    nbyte  = Ref{Cint}(0)
    stream = Ref{Cstring}()
    raiseStatus(ccall((:EG_exportModel, C_egadslib), Cint,(ego, Ptr{Cint}, Ptr{Cstring}),
                model.obj, nbyte,stream))
    str    = unsafe_string(stream[])
    free(stream)
    return str
end



"""
Adds an attribute to the object. If an attribute exists with
the name it is overwritten with the new information.

Parameters
----------
object:
    the Object to attribute

name:
    The name of the attribute. Must not contain a space or other
    special characters

data:
    attribute data to add to the ego
    must be all intergers, all doubles, a string, or a csystem
"""
function attributeAdd!(object::Cego, name::str, data)
    ints  = C_NULL ; reals = C_NULL ; chars = C_NULL
    atype = nothing
    if typeof(data) == CSys
        atype = ATTRCSYS
        reals = data.val
        len   = length(data.val)
    else
        dataF = vcat([data]...)
        len = length(dataF)
        if typeof(dataF) <: vDbl
            atype = ATTRREAL
            reals = Cdouble.(dataF)
        elseif typeof(dataF) <: vInt
            atype = ATTRINT
            ints  = Cint.(dataF)
        elseif typeof(dataF) <: vStr
            atype = ATTRSTRING
            chars = length(dataF) == 1 ? data : dataF
        else
            throw(ErrorMessage(" attributeAdd Data = $data  is not all int, double, or str"))
        end
    end
    raiseStatus(ccall((:EG_attributeAdd, C_egadslib), Cint, (ego, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cstring),
                 object.obj, name, atype, len, ints, reals, chars) )
end


"""
Deletes an attribute from the object.

Parameters
----------
object:
    the Object

name:
    Name of the attribute to delete. If the aname is Nothing then
    all attributes are removed from this object
"""
function attributeDel!(object::Cego; name=nothing)
    cname = emptyVar(name) ? C_NULL : string(name)
    raiseStatus(ccall((:EG_attributeDel, C_egadslib), Cint, (ego, Cstring), object.obj, cname))
end


"""
Returns the number of attributes found with this object

Parameters
----------
object:
    the Object
"""
function attributeNum(object::Cego)
    ptr = Ref{Cint}(0)
    raiseStatus(ccall((:EG_attributeNum, C_egadslib), Cint, (ego, Ptr{Cint}),
                object.obj, ptr))
    return ptr[]
end


"""
Retrieves a specific attribute from the object

Parameters
----------
object:
    the Object to query

name
    the attribute name to query

Returns
-------
attrVal: the appropriate attribute value (of ints, reals or string)
"""
function attributeRet(object::Cego, name::str)
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    icode = ccall((:EG_attributeRet, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
             object.obj, name, atype, len, ints, reals, str)
    icode == EGADS_NOTFOUND && return
    raiseStatus(icode)

    if atype[] == ATTRINT
        return ptr_to_array(ints[],len[])
    elseif atype[] == ATTRREAL
        return ptr_to_array(reals[],len[])
    elseif atype[] == ATTRCSYS
            vr = ptr_to_array(reals[],len[] + 12)
            return CSys(vr[1:len[]], reshape(vr[len[]+1:end], 3,4))
    elseif atype[] == ATTRSTRING
        return unsafe_string(str[])
    end
    return
end


"""
Retrieves a specific attribute from the object

Parameters
----------
object:
    the Object to query

index:
    the index (1 to nAttr from attributeNum)

Returns
-------
attrVal: the appropriate attribute value (of ints, reals or string)
"""
function attributeGet(object::Cego, index::int)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    name  = Ref{Cstring}()
    raiseStatus(ccall((:EG_attributeGet, C_egadslib), Cint, (ego, Cint, Ptr{Cstring}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
                object.obj, Int32(index), name, atype, len, ints, reals, str))

    if atype[] == ATTRINT
        return (name = unsafe_string(name[]),attrVal = ptr_to_array(ints[],len[]))
    elseif atype[] == ATTRREAL
        return (name = unsafe_string(name[]),attrVal = ptr_to_array(reals[],len[]) )
    elseif atype[] == ATTRCSYS
        vr = ptr_to_array(reals[],len[] + 12)
        return (name = unsafe_string(name[]),attrVal = CSys(vr[1:len[]], reshape(vr[len[]+1:end], 3,4)))
    elseif atype[] == ATTRSTRING
        return (name = unsafe_string(name[]),attrVal = unsafe_string(str[]) )
    else
        return (name = nothing, attrVal = nothing)
    end
end


"""
Removes all attributes from the destination object, then copies the
attributes from the source

Parameters
----------
src:
    the source Object

dst:
    the Object to receive src’s attributes
"""
attributeDup!(src::Cego, dst::Cego) = raiseStatus(ccall((:EG_attributeDup, C_egadslib), Cint, (ego, ego), src.obj, dst.obj))


"""
Sets the attribution mode for the Context.

Parameters
----------
context:
    the Context Object

attrFlag:
    the mode flag: 0 - the default scheme, 1 - full attribution node
"""
setFullAttrs!(context::Ccontext, attrFlag::int) = raiseStatus(ccall((:EG_setFullAttrs, C_egadslib), Cint, (ego, Cint), context.ctxt, attrFlag))


"""
Add additional attributes with the same root name to a sequence.

Parameters
----------
object:
    the Object to attribute

name:
    The name of the attribute. Must not contain a space or other
    special characters

data:
    attribute data to add to the ego
    must be all intergers, all dbls, a string, or a csystem

Returns
-------
the new number of entries in the sequence
"""
function attributeAddSeq!(object::Cego, name::str, data)
    ints  = C_NULL ; reals = C_NULL ; chars = C_NULL
    atype = nothing
    if typeof(data) <: CSys
        atype = ATTRCSYS
        reals = data.val
        len = length(reals)
    else
        dataF = vcat([data]...)
        len = length(dataF)
        if typeof(dataF) <: vDbl
            atype = ATTRREAL
            reals = Cdouble.(dataF)
        elseif typeof(dataF) <: vInt
            atype = ATTRINT
            ints  = Cint.(dataF)
        elseif typeof(dataF) <: vStr
            atype = ATTRSTRING
            chars = length(dataF) == 1 ? data : dataF
        else
            throw(ErrorMessage(" attributeAddSeq! Data = $data  is not all int, double, or str"))
        end
    end
    s = ccall((:EG_attributeAddSeq, C_egadslib), Cint, (ego, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cstring),
               object.obj, name, atype, len, ints, reals, chars)
    s < 0 && raiseStatus(s)
    return s
end

"""
Returns the number of named sequenced attributes found on this object.

Parameters
----------
object:
    the Object

name:
    The name of the attribute. Must not contain a space or other
    special characters

Returns
-------
nSeq: the number of sequence attributes with the name
"""
function attributeNumSeq(object::Cego, name::str)
    ptr = Ref{Cint}(0)
    raiseStatus(ccall((:EG_attributeNumSeq, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}),
                object.obj, name, ptr))
    return ptr[]
end


"""
Retrieves a specific attribute from the object

Parameters
----------
object:
    the Object to query

name:
    the “root” name to query. Must not contain a space or other
    special characters

index:
    the sequence number (1 to nSeq)

Returns
-------
attrVal: appropriate attribute value (of ints, reals or string)
"""
function attributeRetSeq(object::Cego, name::str,index::int)
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    raiseStatus(ccall((:EG_attributeRetSeq, C_egadslib), Cint, (ego, Cstring, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
                object.obj, name, index, atype, len, ints, reals, str))
    if atype[] == ATTRINT
        return ptr_to_array(ints[],len[])
    elseif atype[] == ATTRREAL
        return ptr_to_array(reals[],len[])
    elseif atype[] == ATTRCSYS
            vr = ptr_to_array(reals[],len[]+12)
            return CSys(vr[1:len[]], reshape(vr[len[]+1:end], 3,4))
    elseif atype[] == ATTRSTRING
        return unsafe_string(str[])
    end
end

"""
Creates a geometric object:

Parameters
----------
context:
    the Context Object

oclass:
    PCURVE, CURVE or SURFACE

mtype:
    For PCURVE/CURVE:
        LINE, CIRCLE, ELLIPSE, PARABOLA, HYPERBOLA, TRIMMED,
        BEZIER, BSPLINE, OFFSET

    For SURFACE:
        PLANE, SPHERICAL, CYLINDRICAL, REVOLUTION, TORIODAL,
        TRIMMED, BEZIER, BSPLINE, OFFSET, CONICAL, EXTRUSION

reals:
    is the pointer to a block of double precision reals. The
    content and length depends on the oclass/mtype

Optional
--------
ints:
    is a pointer to the block of integer information. Required for
    either BEZIER or BSPLINE

rGeom:
    The reference geometry object (if Nothing use Nothing)

Returns
-------
nGeom: the ego to the new Geometry Object
"""
function makeGeometry(context::Ccontext, oclass::int, mtype::int, reals::Union{Real,vNum, mNum};
                      rGeom=nothing::Union{Nothing,Cego},
                      ints=nothing::Union{Nothing,vInt}, name = "genC"::String)

    reals = Cdouble.(vcat([reals]...))
    if (oclass != CURVE && oclass != PCURVE && oclass != SURFACE)
        throw(ErrorMessage("makeGeometry only accepts CURVE, PCURVE or SURFACE"))
    end
    geo  = rGeom == nothing ? C_NULL : rGeom.obj
    vint = ints  == nothing ? C_NULL : Cint.(ints)

    ptr   = Ref{ego}()
    lr    = length(reals)
    if mtype == BSPLINE
        flag   = ints[1]
        nCP    = oclass == SURFACE ? ints[3] * ints[6] : ints[3]
        nKnots = oclass == SURFACE ? ints[4] + ints[6] : ints[4]
        fact   = oclass == PCURVE  ? 2 : 3
        (flag != 2 && lr != nKnots + fact * nCP) &&
            throw(ErrorMessage(" makeGeometry with periodic BSPLINE: real data for oclass = $oclass needs to be nCP * $fact + nKnots = $nCP * $fact +$nKnots. I have and length(reals) = $lr"))
        (flag == 2 && lr != nKnots + fact * nCP + nCP) &&
            throw(ErrorMessage(" makeGeometry with rational  BSPLINE: real data for oclass = $oclass needs to be nCP * $fact + nKnots + $nCP = $nCP * $fact +$nKnots +$nCP. I have and length(reals) = $lr"))
    elseif mtype == BEZIER
        flag = ints[1]
        nCP  = oclass == SURFACE ? ints[3] * ints[5] : ints[3]
        fact = oclass == PCURVE ? 2 : 3
        (flag != 2 && lr != fact * nCP) &&
            throw(ErrorMessage(" makeGeometry with periodic Bezier: real data for oclass = $oclass needs to be nCP * $fact = $nCP * $fact. I have and length(reals) = $lr"))
        (flag == 2 && lr != fact * nCP + nCP) &&
            throw(ErrorMessage(" makeGeometry with rational Bezier: real data for oclass = $oclass needs to be nCP * $fact + $nCP = $nCP * $fact + $nCP. I have and length(reals) = $lr"))
    end

    raiseStatus(ccall((:EG_makeGeometry, C_egadslib), Cint, (ego, Cint, Cint, ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}),
                context.ctxt, oclass, mtype, geo, vint, reals, ptr))

    return create_Cego(ptr[] ; ctxt =  context, refs = rGeom, name = name)
end

"""
Returns information about the geometric object.
Notes: ints is returned for either mtype = BEZIER or BSPLINE.

Parameters
----------
object:
    the Geometry Object

Returns
-------
oclass:
    the returned Object Class: PCURVE, CURVE or SURFACE

mtype:
    the returned Member Type (depends on oclass)

rGeom:
    the returned reference Geometry Object (Nothing if Nothing)

reals:
    the returned pointer to real data used to describe the geometry

ints:
    the returned pointer to integer information (Nothing if Nothing)
"""
function getGeometry(object::Cego)
    oclass  = Ref{Cint}(0)
    mtype   = Ref{Cint}(0)
    refGeom = Ref{ego}()
    rvec    = Ref{Ptr{Cdouble}}()
    ivec    = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getGeometry, C_egadslib), Cint,(ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}),
              object.obj, oclass, mtype, refGeom, ivec, rvec))

    nivec = Ref{Cint}(0)
    nrvec = Ref{Cint}(0)
    ccall((:EG_getGeometryLen, C_egadslib), Cvoid, (ego, Ptr{Cint}, Ptr{Cint}),object.obj, nivec, nrvec)
    li = nivec[]
    lr = nrvec[]

    j_ivec = ((ivec[] != C_NULL && li > 0) ? ptr_to_array(ivec[], li) : nothing)

    j_rvec = ((rvec[] != C_NULL && lr > 0) ? ptr_to_array(rvec[], lr) : nothing)
    free(ivec) ; free(rvec)

    if mtype[] == BSPLINE
       if oclass[] == PCURVE || oclass[] == CURVE
           flag   = j_ivec[1]
           nCP    = j_ivec[3]
           nKnots = j_ivec[4]
           fact   = (oclass[] == PCURVE ? 2 : 3)
           reals  = j_rvec[1:nKnots+fact*nCP]
           (flag == 2) && (append!(reals,j_rvec[lr-nCP+1:lr]))
       else
           flag    = j_ivec[1]
           nKnots  = j_ivec[4] + j_ivec[7]
           nCP     = j_ivec[3] * j_ivec[6]
           reals   = j_rvec[1:nKnots+3*nCP]
           (flag == 2) && (append!(reals,j_rvec[lr-nCP+1:lr]))
       end
    elseif mtype[] == BEZIER
       fact  = (oclass[] == PCURVE ? 2 : 3)
       nrvec = Int(nrvec[]/fact)
       reals = j_rvec[1:fact * nrvec]
    else
       reals = j_rvec
    end

    return (oclass=oclass[], mtype=mtype[],rGeom=refGeom[], ints=j_ivec, reals = length(reals) == 1 ? reals[1] : reals)
end


"""
This function produces a BSpline Surface that is not fit or approximated
in any way, and is true to the input curves.

Parameters
----------
curves:
    a pointer to a vector of egos containing non-periodic,
    non-rational BSPLINE curves properly positioned and ordered

degree:
    degree of the BSpline used in the skinning direction

Returns
-------
bspline: the new BSpline Surface Object
"""
function skinning(curves::Vector{Cego}, degree::int)
    ptr     = Ref{ego}()
    raiseStatus(ccall((:EG_skinning, C_egadslib), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}),
                length(curves), getfield.(curves,:obj), degree, ptr))
    return create_Cego(ptr[]; ctxt =  curves[1].ctxt)
end


"""
Computes and returns the resultant geometry object created by
approximating the data by a BSpline (OCC or EGADS method).

Parameters
----------
context:
    the Context Object

sizes:
    a vector of 2 integers that specifies the size and dimensionality of
    the data. If the second is zero, then a CURVE is fit and the first
    integer is the length pf the number of points to fit. If the second
    integer is nonzero then the input data reflects a 2D map

xyzs:
    the data to fit -> matrix[3][sizes[1]]

Optional:
-----------
mDeg: default 0

    the maximum degree used by OCC [3-8], or cubic by EGADS [0-2]

    0 - fixes the bounds and uses natural end conditions

    1 - fixes the bounds and maintains the slope input at the bounds

    2 - fixes the bounds & quadratically maintains the slope at 2 nd order

tol: default 1.e-8

    is the tolerance to use for the BSpline approximation procedure,
    zero for a SURFACE fit (OCC).


Returns
-------
bspline: the approximated (or fit) BSpline resultant ego
"""
function approximate(context::Ccontext, sizes::vInt, xyzs::mDbl ;
                     deg=0::Int, tol=1.e-08::dbl)
    ptr      = Ref{ego}()
    raiseStatus(ccall((:EG_approximate, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}),
                context.ctxt, Cint(deg), Cdouble(tol), sizes, vcat(xyzs...), ptr))
    return create_Cego(ptr[] ; ctxt =  context)
end


"""
Computes and returns the resultant geometry object created by approximating
the triangulation by a BSpline surface.

Parameters
----------
context:
    the Context Object

xyzs:
    the coordinates to fit (list of [x,y,z] triads)

tris:
    triangle indices (1 bias) (list of 3-tuples)

Optional
--------
tric:
    neighbor triangle indices (1 bias) -- 0 or (-) at bounds
    Nothing -- will compute (len(tris) in length, if not Nothing)

tol:
    the tolerance to use for the BSpline approximation procedure

Returns
-------
bspline: The approximated (or fit) BSpline resultant ego
"""
function fitTriangles(context::Ccontext, xyzs::Union{vDbl,mDbl}, tris::Union{vInt,mInt} ;
                      tric=nothing::Union{Nothing,vInt},
                      tol=1.e-08::dbl)

    vint = tric == nothing ? C_NULL : tric

    flat_xyz  = vcat(xyzs...)
    flat_tris = vcat(tris...)
    ntris     = Cint(length(flat_tris)/ 3)
    nxyz      = Cint(length(flat_xyz) / 3)

    if vint != C_NULL
        ntric = length(vint)
        Int(ntric  / 3 ) != ntris &&
            (throw(ErrorMessage("fitTriangles tric needs to be integer vector of length $ntris != $ntric")))
    end
    ptr       = Ref{ego}()
    raiseStatus(ccall((:EG_fitTriangles, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{ego}),
                context.ctxt, nxyz,flat_xyz,ntris, flat_tris, vint, tol, ptr))
    return create_Cego(ptr[] ; ctxt =  context)
end


"""
Computes from the input Surface and returns the isocline curve.

Parameters
----------
surface:
    the Surface Object

iUV:
    the type of isocline: UISO (0) constant U or VISO (1) constant V

value:
    the value used for the isocline

Returns
-------
curve: the returned resultant curve object
"""
function isoCline(surface::Cego, iUV::int, value::dbl)
    c = Ref{ego}()
    raiseStatus(ccall((:EG_isoCline, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{ego}),
                surface.obj, iUV, value, c))
    return create_Cego(c[] ; ctxt =  surface.ctxt)
end


"""
Returns the valid range of the object: may be one of PCURVE,
CURVE, SURFACE, EDGE or FACE

Parameters
----------
object:
    the input Object (PCURVE, CURVE, EDGE, EEDGE, SURFACE, FACE or EFACE)


Returns
-------
range:
    for PCURVE, CURVE or EDGE returns 2 values, t-start and t-end

    for SURFACE or FACE returns [2][2] values: ran[1,:] = [u-min, u-max] ; ran[2,:] = [v-min and v-max]

periodic:
    0 for non-periodic, 1 for periodic in t or u 2 for periodic in

    v (or-able)
"""
function getRange(object::Cego)
    geom = getInfo(object)
    r = ( (geom.oclass == PCURVE || geom.oclass == CURVE || geom.oclass == EDGE
        || geom.oclass == EEDGE) ? zeros(2) : zeros(4))
    p = Ref{Cint}(0)
    raiseStatus(ccall((:EG_getRange, C_egadslib), Cint,(ego, Ptr{Cdouble}, Ptr{Cint}),
                object.obj, r,p))
    ran = length(r) == 2 ? r : [r[[1,3]] r[[2,4]]]
    return (range = ran, periodic = p[])
end


"""
Get the arc-length of the  dobject

Parameters
----------
object:
    the input Object (PCURVE, CURVE or EDGE)

t1:
    The starting t value

t2:
    The end t value

Returns
-------
alen = the arc-length
"""
function arcLength(object::Cego, t1::Real, t2::Real)
    alen = Ref{Cdouble}(0.0)
    raiseStatus(ccall((:EG_arcLength, C_egadslib), Cint, (ego, Cdouble, Cdouble, Ptr{Cdouble}),
                object.obj, t1, t2, alen))
    return alen[]
end


"""
Returns the result of evaluating the object at a parameter
point. May be used for PCURVE, CURVE, SURFACE, EDGE or FACE.

Parameters
----------
object:
    The input Object

params:
    The parametric location

    For NODE: Nothing

    For PCURVE, CURVE, EDGE, EEDGE: The t-location

    For SURFACE, FACE, EFACE: The (u, v) coordiantes


Returns
-------
    For NODE: X

    For PCURVE, CURVE, EDGE: X, dX and d2X

    For SURFACE, FACE: X, dX = [X_u,X_v], ddX = [X_uu, X_uv, X_vv]

"""
function evaluate(object::Cego, params::Union{Real,vNum, Nothing}=nothing)

    info = getInfo(object)
    if info.oclass == NODE
        xyz = zeros(3)
    elseif info.oclass == PCURVE
        xyz = zeros(6)
    elseif info.oclass == EDGE || info.oclass == EEDGE ||
           info.oclass == CURVE
           xyz = zeros(9)
    elseif info.oclass == FACE || info.oclass == EFACE ||
           info. oclass == SURFACE
           xyz = zeros(18)
    else
        @info " evaluate object class $info.oclass  Not accepted geometry "
        raiseStatus(EGADS_GEOMERR)
    end
    param = params == nothing ? C_NULL : vcat([params]...)
    raiseStatus(ccall((:EG_evaluate, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object.obj, param, xyz))
    # MAYBE RESHAPE INTO  3 x types
    if info.oclass == NODE
        return xyz
    elseif info.oclass == PCURVE
        return (X = xyz[1:2], dX = xyz[3:4], ddX = xyz[5:6])
    elseif (info.oclass == EDGE || info.oclass == EEDGE ||
            info.oclass == CURVE)
        return (X = xyz[1:3], dX = xyz[4:6], ddX = xyz[7:9])
    else
        return (X = xyz[1:3], dX = [xyz[4:6],xyz[7:9]],
                ddX = [xyz[10:12],xyz[13:15],xyz[16:18]])
    end
end


"""
Returns the result of inverse evaluation on the object.
For topology the result is limited to inside the EDGE/FACE valid bounds.

Parameters
----------
object:
    the input Object

pos:
    [u,v] for a PCURVE and [x,y,z] for all others

Returns
-------
parms:
    the returned parameter(s) found for the nearest position on the object:

    for PCURVE, CURVE, EDGE, or EEDGE the one value is t

    for SURFACE, FACE, or EFACE the 2 values are u then v

result:
    the closest position found is returned:

    [u,v] for a PCURVE (2) and [x,y,z] for all others (3)

Note: When using this with a FACE the timing is significantly slower
      than making the call with the FACE's SURFACE (due to the clipping).
      If you don't need the limiting call EG_invEvaluate with the
      underlying SURFACE.
"""
function invEvaluate(object::Cego, pos::vNum)
    r = zeros(3)
    p = zeros(2)
    raiseStatus(ccall((:EG_invEvaluate, C_egadslib), Cint,(ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                object.obj, Cdouble.(pos), p,r))
    geoPtr = getInfo(object)

    if geoPtr.oclass == PCURVE
        return (params = p[1], result = [r[1], r[2]])
    elseif (geoPtr.oclass == EDGE ||
            geoPtr.oclass == CURVE)
        return (params = p[1], result = r)
    elseif (geoPtr.oclass == FACE ||
            geoPtr.oclass == SURFACE)
        return (params=p, result=r)
    end
    return (params=nothing, result=nothing)
end


"""
Returns the curvature and principle directions/tangents

Parameters
----------
object:
    the input Object

params:
    The parametric location

    For PCURVE, CURVE, EDGE, EEDGE: The t-location

    For SURFACE, FACE, EFACE: The (u, v) coordiantes

Returns
-------
For PCURVE: curvature, direction [dir.x,dir.y]

For CURVE, EDGE: 1 curvature direction -> [dir.x ,dir.y, dir.z]

For SURFACE, FACE: 2 curvature directions -> matrix[2][3] :  dir[1,:] = [dir1.x, dir1.y, dir1.z]
"""
function curvature(object::Cego,params::Union{Real,vNum})
    param  = vcat(params...) #allows feeding pars as a double (not vector)
    geoPtr = unsafe_load(object.obj)
    if geoPtr.oclass == PCURVE
        res = zeros(3)
    elseif (geoPtr.oclass == EDGE  ||
            geoPtr.oclass == EEDGE ||
            geoPtr.oclass == CURVE)
        res = zeros(4)
    elseif (geoPtr.oclass == FACE ||
            geoPtr.oclass == EFACE ||
            geoPtr.oclass == SURFACE)
        res = zeros(8)
    else
        @info "curvature object oclass $geoPtr.oclass not recognised "
        raiseStatus(EGADS_GEOMERR)
    end
    raiseStatus(ccall((:EG_curvature, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object.obj, param, res))
    if length(res) != 8
        return (curvature = res[1], direction = res[2:end])
    else
        return (curvature = res[[1,5]], direction = [res[2:4] res[6:8]]' )
    end
end



"""
Computes and returns the other curve that matches the input curve.
If the input curve is a PCURVE, the output is a 3D CURVE (and vice versa).

Parameters
----------
object:
    the input Object (SURFACE or FACE)

curve:
    the input PCURVE or CURVE/EDGE object

tol:
    is the tolerance to use when fitting the output curve

Returns
-------
ocurve: the approximated resultant curve object
"""
function otherCurve(surface::Cego, curve::Cego , tol::dbl=0.0)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_otherCurve, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{ego}),
                surface.obj, curve.obj, tol,ptr))
    return create_Cego(ptr[] ; ctxt =  surface.ctxt)
end


"""
Compares two objects for geometric equivalence.

Parameters
----------
object1:
    the first input Object (NODE, CURVE, EDGE, SURFACE or FACE)

object2:
    an object to compare with
"""
isSame(object1::Cego, object2::Cego) = return EGADS_SUCCESS == ccall((:EG_isSame, C_egadslib), Cint, (ego, ego),object1.obj,object2.obj)


"""
Computes and returns the BSpline representation of the geometric object

Parameters
----------
object:
    the input Object (PCURVE, CURVE, EDGE, SURFACE or FACE)

Returns
-------
bspline: the approximated resultant BSPLINE object
can be the input ego when a BSPLINE and the same limits
"""
function convertToBSpline(object::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_convertToBSpline, C_egadslib), Cint, (ego, Ptr{ego}), object.obj, ptr))
    return create_Cego(ptr[] ; ctxt =  object.ctxt, delObj = (object.obj != ptr[] ))
end


"""
Computes and returns the BSpline representation of the geometric object

Required when converting Geometry Objects with infinite range.

Parameters
----------
object:
    the input Object (PCURVE, CURVE or SURFACE)

range:
    t range (2) or [u, v] box (4) to limit the conversion

Returns
-------
bspline: the returned approximated resultant BSPLINE object
can be the input ego when a BSPLINE and the same limits
"""
function convertToBSplineRange(geom::Cego, range::Union{vNum,mNum})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_convertToBSplineRange, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
                geom.obj, range, ptr))
    return create_Cego(ptr[] ; ctxt =  geom.ctxt)
end


"""
Creates and returns a topological object:

Parameters
----------
context:
    the Context Object

oclass:
    either NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

Optional
--------
mtype:
    for EDGE is TWONODE, ONENODE or DEGENERATE

    for LOOP is OPEN or CLOSED

    for FACE is either SFORWARD or SREVERSE

    for SHELL is OPEN or CLOSED

    BODY is either WIREBODY, FACEBODY, SHEETBODY or SOLIDBODY

geom:
    reference geometry object required for EDGEs and FACEs
    (optional for LOOP)

reals:
    may be Nothing except for:

    NODE which contains the [x,y,z] location

    EDGE is [t-min t-max] (the parametric bounds)

    FACE is [u-min, u-max, v-min, v-max] (the parametric bounds)

children:
    chldrn a list of children objects (nchild in length)
    if LOOP and has reference SURFACE, then 2*nchild in length
    (PCURVES follow)

senses:
    a list of integer senses for the children (required for FACES
    & LOOPs only)

Returns
-------
topo: resultant new topological ego
"""
function makeTopology!(context::Ccontext, oclass::int;
                       mtype=0::Int, geom=nothing::Union{Nothing, Cego},
                       reals=nothing::Union{Nothing,vDbl,mDbl},
                       senses=nothing::Union{Nothing,vInt},
                       name="genCtopo"::String,
                       children=nothing::Union{Nothing, Cego, Vector{Cego}})
    geo = geom   == nothing ? C_NULL : geom.obj
    rea = reals  == nothing ? C_NULL : Cdouble.(vcat(reals...))
    sen = senses == nothing ? C_NULL : Cint.(vcat(senses...))
    childObj = C_NULL ; nchild = 0
    childC   = vcat([children]...)

    if childC != [nothing]
        nchild = length(childC)
        childObj = getfield.(childC,:obj)
    end

    mtype = Cint(mtype)
    topo_ptr = Ref{ego}()
    lR = length(rea)
    if (oclass == NODE && lR != 3)
        throw(ErrorMessage("NODE must be have reals of length 3 != $lR"))
    elseif oclass == EDGE
        (mtype != TWONODE && mtype != ONENODE && mtype != CLOSED) &&
            throw(ErrorMessage("EDGE must be TWONODE, ONENODE or CLOSED"))
        (lR != 2) &&
            throw(ErrorMessage("EDGE must be have data of length 2 != $lR"))
    elseif oclass == LOOP
        (mtype != OPEN && mtype != CLOSED) &&
            throw(ErrorMessage("LOOP must be OPEN or CLOSED"))
        (sen == C_NULL) &&
            (throw(ErrorMessage("LOOP needs senses !! ")))
    elseif (oclass == FACE  && mtype != SFORWARD  && mtype != SREVERSE)
            throw(ErrorMessage("FACE must be SFORWARD or SREVERSE"))
    elseif (oclass == SHELL && mtype != OPEN      && mtype != CLOSED)
            throw(ErrorMessage("SHELL must be OPEN or CLOSED"))
    elseif (oclass == BODY  && mtype != WIREBODY  && mtype != FACEBODY   &&
                               mtype != SHEETBODY && mtype != SOLIDBODY)
            throw(ErrorMessage("BODY must be WIREBODY, FACEBODY,SHEETBODY or SOLIDBODY"))
    end

    # Store references to children so they will be deleted after the new ego
    refs = nchild == 0 ? nothing : childC

    # Account for effective topology
    if oclass == MODEL
        mtype  = nchild
        nchild = 0
        for i in childC
            info = getInfo(i)
            info.oclass == BODY && (nchild += 1)
        end
    end

    # children contains both edges and p-curves
    (oclass == LOOP && geo != C_NULL) && (nchild = Int(nchild/2))
    nsenses = sen == C_NULL ? 0 : length(sen)
    (nchild > 0 && nsenses > 0 && nchild != nsenses) &&
        throw(ErrorMessage(" Number of children is inconsistent with senses $nchild != $nsenses"))


     raiseStatus(ccall((:EG_makeTopology, C_egadslib), Cint,
     (ego, ego, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{ego}, Ptr{Cint}, Ptr{ego}),
      context.ctxt, geo, oclass, mtype, rea, nchild,childObj, sen, topo_ptr))

    (geo != C_NULL) && push!(refs,geom)

    if oclass == BODY
        refs = nothing
    elseif oclass == MODEL
        setfield!.(childC, :delObj, false)
    end

    return create_Cego(topo_ptr[] ; ctxt =  context, name = name, refs = refs)
end



"""
Returns information about the topological object.

Parameters
----------
topo:
    the Topology or Effective Topology Object to query

Returns
-------
geom:
    The reference geometry object (if Nothing this is returned as Nothing)

oclass:
    is NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

mtype:
    for EDGE is TWONODE, ONENODE or DEGENERATE

    for LOOP is OPEN or CLOSED

    for FACE is either SFORWARD or SREVERSE

    for SHELL is OPEN or CLOSED

    BODY is either WIREBODY, FACEBODY, SHEETBODY or SOLIDBODY

reals:
    will retrieve at most 4 doubles:

    for NODE this contains the [x,y,z] location

    EDGE is the t-min and t-max (the parametric bounds)

    FACE returns the [u,v] box (the limits first for u then for v)

    number of children (lesser) topological objects

children:
    the returned pointer to a vector of children objects

    if a LOOP with a reference SURFACE, then 2*nchild in length (PCurves follow)

    if a MODEL -- nchild is the number of Body Objects, mtype the total

senses:
    is the returned pointer to a block of integer senses for

    the children.
"""
function getTopology(topo::Cego)
    geom   = Ref{ego}()
    oclass = Ref{Cint}(0)
    mtype  = Ref{Cint}(0)
    lims   = zeros(4)
    nChild = Ref{Cint}(0)
    sen    = Ref{Ptr{Cint}}(0)
    child  = Ref{Ptr{ego}}()

    raiseStatus(ccall((:EG_getTopology, C_egadslib), Cint, (ego, Ptr{ego},
                 Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint},Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}),
                 topo.obj, geom, oclass, mtype,lims, nChild, child, sen))

    nC = (oclass[] == LOOP  && geom[] != C_NULL) ? nChild[] * 2 :
         (oclass[] == MODEL && mtype[] > 0     ) ? mtype[] : nChild[]

    aux = ptr_to_array(child[], nC)

    children = [create_Cego(unsafe_load(child[], j) ; ctxt =  topo.ctxt, delObj = false) for j = 1:nC]

    senses   = ptr_to_array(sen[], nChild[])

    reals =  oclass[] == NODE ? lims[1:3] : oclass[] == EDGE ? lims[1:2] :
            (oclass[] == FACE || oclass[] == EFACE) ? [lims[[1,3]] lims[[2,4]]] : nothing

    return (geom = create_Cego(geom[] ; ctxt =  topo.ctxt, delObj = false),
            oclass = oclass[], mtype=mtype[], reals=reals, senses = senses, children = children)
end


"""
Creates a simple FACE from a LOOP or a SURFACE. Also can be
used to hollow a single LOOPed existing FACE. This function
creates any required NODEs, EDGEs and LOOPs.

Parameters
----------
object:
    Either a LOOP (for a planar cap), a SURFACE with [u,v] bounds,
    or a FACE to be hollowed out

mtype:
    Is either SFORWARD or SREVERSE
    For LOOPs you may want to look at the orientation using
    getArea, ignored when the input object is a FACE

rdata:
    May be Nothing for LOOPs, but must be the limits for a SURFACE (4
    values), the hollow/offset distance and fillet radius (zero is
    for no fillets) for a FACE input object (2 values)

Returns
-------
face: the resultant returned topological FACE object (a return of
      EGADS_OUTSIDE is the indication that offset distance was too
      large to produce any cutouts, and this result is the input
      object)
"""
function makeFace(object::Cego, mtype::int ; rdata::Union{mNum,vNum,Nothing} = nothing)
    limits = typeof(rdata) <: mNum ? Cdouble.(vcat(rdata'...)) :
             typeof(rdata) <: vNum ? Cdouble.(rdata) : C_NULL

    ptr    = Ref{ego}()
    icode  = ccall((:EG_makeFace, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}),
                    object.obj, mtype, limits, ptr)

    icode < 0 && raiseStatus(icode)

    return create_Cego(ptr[]; ctxt =  object.ctxt, refs = object)
end


"""
Creates a LOOP from a list of EDGE Objects, where the EDGEs do
not have to be topologically connected. The tolerance is used
to build the NODEs for the LOOP. The orientation is set by the
first non-NULL entry in the list, which is taken in the
positive sense. This is designed to be executed until all list
entries are exhausted.

Parameters
----------
edges:
    list of EDGEs, of which some may be NULL
    !! The edges used in loop are converted to Cego()

Optional
--------
geom:
    SURFACE Object for non-planar LOOPs to be used to bound
    FACEs (can be Nothing)

toler:
    tolerance used for the operation (0.0 - use EDGE tolerances)

Returns
-------
loop:
    the resultant LOOP Object

nedges:
    the number of non- entries in edges when returned
    or error code
"""
function makeLoop!(edges::Union{Cego, Vector{Cego}};
                   geom = nothing::Union{Nothing,Cego},
                   toler = 0.0::dbl, name ="genLoop"::String)

    geo  = geom == nothing ? C_NULL : geom.obj
    nedge    = length(edges)
    loop_ptr = Ref{ego}()

    # Create a copy of the edge pointers
    piv  = findall(x -> !emptyVar(x), edges )
    # Create a copy of the edge pointers
    ePtr = [edges[piv[i]].obj for i = 1:length(piv) ]

    ctxt = edges[piv[1]].ctxt

    # This function makes ePtr .== C_NULL !!
    s    = ccall((:EG_makeLoop, C_egadslib), Cint,(Cint, Ptr{ego}, ego, Cdouble, Ptr{ego}),
                nedge, ePtr, geo, toler, loop_ptr)

    piv        = findall(x -> emptyVar(x), ePtr)

    eRefs      = deepcopy.(edges[piv])
    edges[piv].= create_Cego.()

    loop = create_Cego(loop_ptr[]; ctxt = ctxt, name = name, refs = eRefs)
    s < EGADS_SUCCESS && raiseStatus(s)
    return loop
end



"""
Creates a simple SOLIDBODY. Can be either a box, cylinder,
sphere, cone, or torus.

Parameters
----------
context:
    the Context Object

stype:
    BOX, SPHERE, CONE, CYLINDER or TORUS

data:
    Depends on stype

    For box [x,y,z] then [dx,dy,dz] for size of box

    For sphere [x,y,z] of center then radius

    For cone apex [x,y,z], base center [x,y,z], then radius

    For cylinder 2 axis points and the radius

    For torus [x,y,z] of center, direction of rotation, then major

    radius and minor radius

returns:
    the resultant topological BODY object
"""
function makeSolidBody(context::Ccontext, stype::int, data::vNum)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeSolidBody, C_egadslib), Cint, (ego, Cint,Ptr{Cdouble}, Ptr{ego}),
                context.ctxt, stype, Cdouble.(data), ptr))
    return create_Cego(ptr[]; ctxt =  context)
end


"""
Returns topologically connected objects:

Parameters
----------
body:
    the Body/EBody Object

oclass:
    is NODE, EGDE, LOOP, FACE or SHELL -- must not be the same class as ref

ref:
    reference topological object or NULL. Sets the context for the
    returned objects (i.e. all objects of a class [oclass] in the
    tree looking towards that class from ref).  starts from
    the BODY (for example all NODEs in the BODY)

Returns
-------
topos: the returned number of requested topological objects is a
       returned pointer to the block of objects
"""
function getBodyTopos(body::Cego, oclass::int; ref::Union{Cego,Nothing} = nothing)

    ntopo = Ref{Cint}(0)
    topos = Ref{Ptr{ego}}()
    ref   = typeof(ref) == Cego ? ref.obj : C_NULL
    raiseStatus( ccall((:EG_getBodyTopos, C_egadslib), Cint,
                (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{ego}}),
                 body.obj, ref, oclass, ntopo, topos))

    topo_arr = ptr_to_array(topos[], ntopo[])
    free(topos)
    return [create_Cego(topo_arr[j]; ctxt =  body.ctxt, delObj = false) for j = 1:ntopo[] ]
end


"""
Return the (1-based) index of the topological object in the BODY/EBODY

Parameters
----------
body:
    the Body/EBody Object

object:
    the Topology Object in the Body/EBody

Returns
-------
the index
"""
indexBodyTopo(body::Cego,object::Cego) = ccall((:EG_indexBodyTopo, C_egadslib), Cint, (ego, ego), body.obj, object.obj)


"""
Returns the topological object (based on index) in the BODY/EBODY

Parameters
----------
body:
    the Body/EBody Object

oclass:
    is NODE, EDGE, LOOP, FACE or SHELL

index:
    is the index (bias 1) of the entity requested

Returns
-------
obj: the topological object in the Body
"""
function objectBodyTopo(body::Cego, oclass::int, index::int)
    ptr = Ref{ego}()
    loadBody = unsafe_load(body.obj)
    raiseStatus(ccall((:EG_objectBodyTopo, C_egadslib), Cint, (ego, Cint, Cint, Ptr{ego}),
                body.obj, Cint(oclass), Cint(index), ptr))
    return create_Cego(ptr[] ; ctxt =  body.ctxt, delObj = false)
end


"""
Computes the surface area from a LOOP, a surface or a FACE.
When a LOOP is used a planar surface is fit
and the resultant area can be negative if the orientation of
the fit is opposite of the LOOP.

Parameters
----------
object:
    either a Loop (for a planar cap), a surface with [u; v] bounds or a Face

limits:
    may be Nothing except must contain the limits for a surface (4 or 2 x 2)

Returns
-------
area: the surface area
"""
function getArea(object; limits = nothing)
    info = getInfo(object)
    area = Ref{Cdouble}(0.0)
    lims = emptyVar(limits) ? C_NULL : vcat(limits'...)
    info = getInfo(object)
    if info.oclass == SURFACE && (lims == C_NULL || length(lims) != 4)
        throw(ErrorMessage("getArea for surfaces needs 4 limits (u = limits[1,:] v = limits[2,:] ) "))
    end
    raiseStatus(ccall((:EG_getArea, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object.obj, lims, area))
    return area[]
end


"""
Computes the Cartesian bounding box around the object:

Parameters
----------
object:
    any topological object

Returns
-------
bbox: [(x_min,y_min,z_min),
       (x_max,y_max,z_max)]
"""
function getBoundingBox(topo::Cego)
    bbox = zeros(6)
    raiseStatus(ccall((:EG_getBoundingBox, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo.obj, bbox))
    return reshape(bbox, 2,3)
end


"""
Computes and returns the physical and inertial properties of a topological object

Parameters
----------
object:
    can be EDGE, LOOP, FACE, SHELL, BODY or Effective Topology counterpart

Returns
-------
volume,

area (length for EDGE, LOOP or WIREBODY),

CG: center of gravity vector[3]

IM: inertia matrix (at CG) [9]
"""
function getMassProperties(object::Cego)
    props = zeros(14)
    raiseStatus(ccall((:EG_getMassProperties, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                object.obj, props))
    return (volume = props[1], area = props[2], CG = props[3:5], IM = props[6:14])
end


"""
Compares two topological objects for equivalence.

Parameters
----------
topo1:
    the first input Topology Object

topo2:
    the second input Topology Object to compare with
"""
isEquivalent(topo1::Cego, topo2::Cego) = return EGADS_SUCCESS == ccall((:EG_isEquivalent, C_egadslib), Cint, (ego, ego), topo1.obj, topo2.obj)


"""
Computes whether the point is on or contained within the object.
Works with EDGEs and FACEs by projection. SHELLs must be CLOSED.

Parameters
----------
object:
    the object, can be EDGE, FACE, SHELL or SOLIDBODY

xyz:
    the coordinate location to check

Returns
-------
True if xyz is in the topology, False otherwise
"""
function inTopology(object::Cego, xyz::vNum)
    s = ccall((:EG_inTopology, C_egadslib), Cint, (ego, Ptr{Cdouble}), object.obj, xyz)
    (s < EGADS_SUCCESS) && raiseStatus(s)
    return s == EGADS_SUCCESS
end


"""
Computes the result of the [u,v] location in the valid part of the FACE.

Parameters
----------
face:
    the Face Object

uv:
    the parametric location to check

Returns
-------
true if uv is in the FACE, false otherwise
"""
function inFace(face::Cego, uv::vNum)
    s = ccall((:EG_inFace, C_egadslib), Cint, (ego, Ptr{Cdouble}), face.obj, uv)
    (s < EGADS_SUCCESS) && (raiseStatus(s) )
    return s == EGADS_SUCCESS
end


"""
Computes on the EDGE/PCURVE to get the appropriate [u,v] on the FACE.

Parameters
----------
face:
    the Face Object

edge:
    the EDGE object

t:
    the parametric value to use for the evaluation

sense:
    can be 0, but must be specified (+/-1) if the EDGE is found
    in the FACE twice that denotes the position in the LOOP to use.
    EGADS_TOPOERR is returned when sense==0 and an EDGE is found twice.

Returns
-------
uv: the resulting [u,v] evaluated at t
"""
function getEdgeUV(face::Cego, edge::Cego, sense::int, t::Real)
    uv = zeros(2)
    raiseStatus(ccall((:EG_getEdgeUV, C_egadslib), Cint, (ego, ego, Cint, Cdouble, Ptr{Cdouble}),
                face.obj, edge.obj, Cint.(sense), t, uv))
    return uv
end


"""
Creates a MODEL from a collection of Objects. The Objects can
be either BODYs (not WIREBODY), SHELLs and/or FACEs. After the
sewing operation, any unconnected Objects are returned as
BODYs.

Parameters
----------
objects:
    List of Objects to sew together

Optional
--------
toler:
    Tolerance used for the operation (0.0 - use Face tolerances)

flag:
    Indicates whether to produce manifold/non-manifold geometry

Returns
-------
model: The resultant MODEL object
"""
function sewFaces(objects::Vector{Cego}; toler=0.0::dbl, flag = 1::Int)

    ptr   = Ref{ego}()
    raiseStatus(ccall((:EG_sewFaces, C_egadslib), Cint, (Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}),
                 length(objects),getfield.(objects,:obj), toler, flag, ptr))
    return create_Cego(ptr[] ; ctxt =  objects[1].ctxt, refs = objects)
end


"""
Creates a new SHEETBODY or SOLIDBODY from an input SHEETBODY or SOLIDBODY
and a list of FACEs to modify. The FACEs are input in pairs where
the first must be an Object in the BODY and the second either
a new FACE or Nothing. The Nothing replacement flags removal of the FACE in the BODY.

Parameters
----------
body:
    the Body Object to adjust (either a SHEETBODY or a SOLIDBODY)

faces:
    list of FACE tuple pairs, where the first must be a FACE in the BODY
    and second is either the FACE to use as a replacement or a Nothing which
    indicates that the FACE is to be removed from the BODY

Returns
-------
result: the resultant BODY object, either a SHEETBODY or a SOLIDBODY
        (where the input was a SOLIDBODY and all FACEs are replaced
        in a way that the LOOPs match up)
"""
function replaceFaces(body::Cego, objects)

    obj_flat = vcat(objects...)
    nfaces   = Int32(length(obj_flat) / 2)
    (length(obj_flat) != 2 * nfaces) &&
        (throw(ErrorMessage("replaceFaces object is not a pair ")))
    pfaces = Array{ego}(undef, 2 * nfaces)
    for i = 1: 2 * nfaces
        pfaces[i] = emptyVar(obj_flat[i]) ? C_NULL : obj_flat[i].obj
    end
    res    = Ref{ego}()
    info   = getInfo(body)
    raiseStatus(ccall((:EG_replaceFaces, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
                       body.obj,nfaces, pfaces,res))
    return create_Cego(res[] ; ctxt =  body.ctxt)
end


"""
Get the body that this object is contained within

Parameters
----------
object:
    the input topology Object (can be an EBody to retrieve the referenced Body)

Returns
-------
body: The BODY ego
"""
function getBody(object::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_getBody, C_egadslib), Cint, (ego, Ptr{ego}),object.obj,ptr))
    return create_Cego(ptr[] ; ctxt =  object.ctxt, delObj = false)
end


"""
Parameters
----------
topo:
    the Topology Object

Returns
-------
toler: the maximum tolerance defined for the object's hierarchy
"""
function tolerance(topo::Cego)
    tol = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_tolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo.obj, tol))
    return tol[]
end


"""
Parameters
----------
topo:
    the Topology Object

Returns
-------
toler: the internal tolerance defined for the object
"""
function getTolerance(topo::Cego)
    tol = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_getTolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo.obj, tol))
    return tol[]
end


"""
Examines the EDGEs in one BODY against all of the EDGEs in another.
If the number of NODEs, the NODE locations, the EDGE bounding boxes
and the EDGE arc lengths match it is assumed that the EDGEs match.
A list of pairs of indices are returned.

Parameters
----------
body1:
    the first input Body Object

body2:
    the second Body Object

toler:
    the tolerance used (can be zero to use entity tolerances)

Returns
-------
matches: a list of the tuples of matching indices
"""
function matchBodyEdges(body1::Cego, body2::Cego ; toler::dbl=0.0)
    nMatch  = Ref{Cint}(0)
    matches = Ref{Ptr{Cint}}()
    v1 = unsafe_load(body1.obj)
    v2 = unsafe_load(body2.obj)

    (v1.oclass != BODY || v2.oclass != BODY) &&
        (throw(ErrorMessage("matchBodyFaces Both objects must have object class BODY")))
    raiseStatus(ccall((:EG_matchBodyEdges, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}),
                body1.obj, body2.obj, toler, nMatch, matches))
    arr  = ptr_to_array(matches[],nMatch[],2)
    free(matches)
    return arr
end


"""
Examines the FACEs in one BODY against all of the FACEs in
another. If the number of LOOPs, number of NODEs, the NODE
locations, the number of EDGEs and the EDGE bounding boxes as
well as the EDGE arc lengths match it is assumed that the
FACEs match. A list of pairs of indices are returned.

Parameters
----------
body1:
    the first input Body Object

body2:
    the second Body Object

toler:
    the tolerance used (can be zero to use entity tolerances)

Returns
-------

matches: a list of the tuples of matching indices

Note: This is useful for the situation where there are
      glancing FACEs and a UNION operation fails (or would
      fail). Simply find the matching FACEs and do not include them
      in a call to EG_sewFaces.
"""
function matchBodyFaces(body1::Cego, body2::Cego; toler::dbl=0.0)
    nMatch  = Ref{Cint}(0)
    matches = Ref{Ptr{Cint}}()
    v1 = unsafe_load(body1.obj)
    v2 = unsafe_load(body2.obj)

    (v1.oclass != BODY || v2.oclass != BODY) &&
        (throw(ErrorMessage("matchBodyFaces Both objects must have object class BODY")))

    raiseStatus(ccall((:EG_matchBodyFaces, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}),
                body1.obj, body2.obj, toler, nMatch, matches))
    arr = ptr_to_array(matches[],nMatch[],2)
    free(matches)
    return arr
end


"""
Checks for topological equivalence between the the BODY self and the BODY dst.
If necessary, produces a mapping (indices in src which map to dst) and
places these as attributes on the resultant BODY mapped (named .nMap, .eMap and .fMap).
Also may modify BSplines associated with FACEs.

Parameters
----------
src:
    the source Body Object (not a WIREBODY)

dst:
    destination body object

fAttr:
    the FACE attribute string used to map FACEs

Returns
-------
mapped: the mapped resultant BODY object copied from dst
        If Nothing and s == EGADS_SUCCESS, dst is equivalent
        and can be used directly in EG_mapTessBody

Note: It is the responsibility of the caller to have
      uniquely attributed all FACEs in both src and dst
      to aid in the mapping for all but FACEBODYs.
"""
function mapBody(src::Cego, dst::Cego, fAttr::str)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_mapBody, C_egadslib), Cint, (ego, ego, Cstring, Ptr{ego}),
                src.obj, dst.obj, fAttr, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Computes the Winding Angle along an Edge

The Winding Angle is measured from one Face ``winding''
around to the other based on the normals.
An Edge with a single Face always returns 180.0.

Parameters
----------
edge:
    the input Edge (not DEGENERATE)
    returns an error if there are 3 or more Faces attached to the Edge

t:
    the t-value along the Edge used to compute the Winding Angle

Returns
-------
angle: Winding Angle in degrees (0.0 -- 360.0)
"""
function getWindingAngle(edge::Cego, t::Real)
    ptr = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_getWindingAngle, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}),
                edge.obj, t, ptr))
    return ptr[]
end


"""
Creates a discretization object from a geometry-based Object.

Parameters
----------
geom:
    the input Object, may be a CURVE or SURFACE

limits:
    the bounds of the tessellation (like range)

sizes:
    a set of 2 integers that specifies the size and dimensionality of the
    data. The second is assumed zero for a CURVE and in this case
    the first integer is the length of the number of evenly spaced (in t)
    points created. The second integer must be nonzero for SURFACEs and
    this then specifies the density of the [u,v] map of coordinates
    produced (again evenly spaced in the parametric space).
    If a value of sizes is negative, then the fill is reversed for that coordinate.

Returns
-------
tess: the returned resultant Tessellation of the geometric entity
"""
function makeTessGeom(geom::Cego, limits::vNum, sizes::vInt)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeTessGeom, C_egadslib), Cint,
                (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{ego}),
                geom.obj, Float64.(vcat(limits...)), Cint.(sizes), ptr))
    return create_Cego(ptr[] ; ctxt =  geom.ctxt, refs = geom)
end


"""
Retrieves the data associated with the discretization of a geometry-based Object.

Parameters
----------
tess:
    the input geometric Tessellation Object

Returns
-------
sizes: a returned set of 2 integers that specifies the size and dimensionality
       of the data. If the second is zero, then it is from a CURVE and
       the first integer is the length of the number of [x,y,z] triads.
       If the second integer is nonzero then the input data reflects a 2D map of coordinates.

xyz: the returned pointer to the suite of coordinate data -> matrix[3][sizes]
"""
function getTessGeom(tess::Cego)
    sizes = zeros(Cint,2)
    xyz   = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_getTessGeom, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cdouble}}),
                tess.obj, sizes, xyz))
    return (sizes = sizes, xyzs = ptr_to_array(xyz[] , sizes[1] * sizes[2] , 3) )
end


"""
Creates a discretization object from a Topological BODY Object.

Parameters
----------
body:
    the input Body or closed EBody Object, may be any Body type

parms:
    a set of 3 parameters that drive the EDGE discretization and the FACE
    triangulation. The first is the maximum length of an EDGE segment
    or triangle side (in physical space). A zero is flag that allows for
    any length. The second is a curvature-based value that looks locally
    at the deviation between the centroid of the discrete object and the
    underlying geometry. Any deviation larger than the input value will
    cause the tessellation to be enhanced in those regions. The third is
    the maximum interior dihedral angle (in degrees) between triangle
    facets (or Edge segment tangents for a WIREBODY tessellation), note
    that a zero ignores this phase.

Returns
-------
tess: the resulting egads tessellation object
"""
function makeTessBody(body::Cego, parms::vNum)
    ptr   = Ref{ego}()
    raiseStatus(ccall((:EG_makeTessBody, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
            body.obj, Float64.(parms),ptr))
    return create_Cego(ptr[] ; ctxt =  body.ctxt, refs = body)
end


"""
Takes a triangulation as input and outputs a full quadrilateralization
of the Body. The algorithm uses the bounds of each Face (the discretization
of the Loops) and drives the interior towards regularization (4 quad sides
attached to a vertex) without regard to spacing but maintaining a valid mesh.
This is the recommended quad approach in that it is robust and does not require
manual intervention like makeQuads (plus retrieving the quads is much simpler
and does nor require invoking getQuads and getPatch).

Parameters
----------
tess:
    the input Tessellation Object (not for WIREBODY)

Returns
-------
qtess: a triangle-based Tessellation Object, but with pairs of triangles
       (as retrieved by calls to getTessFace) representing each quadrilateral.
       This is marked by the following attributes:
       ".tessType" (ATTRSTRING) is set to "Quad"
       ".mixed" with type ATTRINT and the length is the number of Faces in the Body,
        where the values are the number of quads (triangle pairs) per Face
"""
function quadTess(tess::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_quadTess, C_egadslib), Cint, (ego, Ptr{ego}), tess.obj, ptr))
    return create_Cego(ptr[] ; ctxt =  tess.ctxt, refs = tess)
end


"""
Redoes the discretization for specified objects from within a BODY TESSELLATION.

Parameters
----------
tess:
    the Tessellation Object to modify

facedg:
    list of FACE and/or EDGE objects from within the BODY used to
    create the TESSELLATION object. First all specified Edges are
    rediscretized. Then any listed Face and the Faces touched by
    the retessellated Edges are retriangulated.
    Note that Quad Patches associated with Faces whose Edges were
    redone will be removed.

parms:
    a set of 3 parameters that drive the EDGE discretization and the FACE
    triangulation. The first is the maximum length of an EDGE segment
    or triangle side (in physical space). A zero is flag that allows for
    any length. The second is a curvature-based value that looks locally
    at the deviation between the centroid of the discrete object and the
    underlying geometry. Any deviation larger than the input value will
    cause the tessellation to be enhanced in those regions. The third is
    the maximum interior dihedral angle (in degrees) between triangle
    facets (or Edge segment tangents for a WIREBODY tessellation), note
    that a zero ignores this phase.
"""
remakeTess!(tess::Cego, facedg::Vector{Cego}, parms) = raiseStatus(ccall((:EG_remakeTess, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{Cdouble}),
                                                tess.obj, length(facedg), getfield.(facedg,:obj), Float64.(parms)))


"""
Maps the input discretization object to another BODY Object.
The topologies of the BODY that created the input tessellation
must match the topology of the body argument (the use of EG_mapBody
can be used to assist).

Parameters
----------
tess:
    the Body Tessellation Object used to create the tessellation on body

body:
    the BODY object (with a matching Topology) used to map the
    tessellation.

Returns
-------
newTess: the resultant TESSELLATION object. The triangulation is
         simply copied but the uv and xyz positions reflect
         the input body (above).

Note: Invoking EG_moveEdgeVert, EG_deleteEdgeVert and/or EG_insertEdgeVerts
      in the source tessellation before calling this routine invalidates
      the ability of EG_mapTessBody to perform its function.
"""
function mapTessBody(tess::Cego, body::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_mapTessBody, C_egadslib), Cint, (ego, ego, Ptr{ego}),
                tess.obj, body.obj, ptr))
    return create_Cego(ptr[] ; ctxt =  tess.ctxt, refs = body)
end


"""
Provides the triangle and the vertex weights for each of the input requests
or the evaluated positions in a mapped tessellation.

Parameters
----------
tess:
    the input Body Tessellation Object

ifaces:
    the face indices for each request - minus index refers to the use of a
    mapped Face index from EG_mapBody and EG_mapTessBody

uv:
    the UV positions in the face for each request (2*len(ifaces) in length)

mapped:
    perform mapped evaluations


Returns
-------
weight: the vertex weights in the triangle that refer to the requested position
        (any negative weight indicates that the point was extrapolated)
        -or-
        the evaluated position based on the input uvs (when mapped is True) -> matrix[3][npt]

tris: the resultant 1-bias triangle index (if not mapped)
      if input as NULL then this function will perform mapped evaluations
"""
function locateTessBody(tess::Cego, ifaces::vInt, uv::Union{vNum, mNum}, mapped = false)

    npt     = length(ifaces)
    weights = zeros(3 * npt)
    itri    = mapped ? zeros(Cint, npt) : C_NULL
    idx     = collect(Cint, vcat(ifaces...))
    raiseStatus(ccall((:EG_locateTessBody, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
                tess.obj, npt, idx, vcat(uv...), itri, weights))
    return mapped ? (weight = reshape(weights,3,npt), tris = itri) : reshape(weights,3,npt)
end


"""
Retrieves the data associated with the discretization of an EDGE
from a Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

eIndex:
    the EDGE index (1 bias). The EDGE Objects and number of EDGEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
    A minus refers to the use of a mapped (+) Edge index from applying
    the functions mapBody and mapTessBody


Returns
-------
xyzs: the set of coordinate data for each vertex -> matrix[3][len]

ts: the vector to the parameter values associated with each vertex -> matrix[2][len]
"""
function getTessEdge(tess::Cego, eIndex::int)
    len = Ref{Cint}(0)
    xyz = Ref{Ptr{Cdouble}}()
    t   = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_getTessEdge, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}),
                tess.obj, Cint(eIndex), len, xyz, t))
    t   = ptr_to_array(t[],len[])
    return (xyzs= ptr_to_array(xyz[],len[],3), ts=t)
end


"""
Retrieves the data associated with the discretization of a FACE from a
Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

fIndex:
    the FACE index (1 bias). The FACE Objects and number of FACEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
    A minus refers to the use of a mapped (+) FACE index (if it exists)

Returns
-------
xyz: set of coordinate data for each vertex [3][len]

uv: parameter values associated with each vertex [2][len]

ptype: vertex type (-1 - internal, 0 - NODE, >0 EDGE)

pindx: vertex index (-1 internal)

tris: returned matrix to triangle vertex indices -> matrix[3][ntri]
      orientation consistent with FACE's mtype

tric: returned matrix to triangle neighbor information looking at opposing side -> matrix[3][ntri]:
      negative is Edge index for an external side
"""
function getTessFace(tess::Cego, fIndex::int)
    len    = Ref{Cint}(0)
    ntri   = Ref{Cint}(0)

    xyz    = Ref{Ptr{Cdouble}}()
    uv     = Ref{Ptr{Cdouble}}()
    ptype  = Ref{Ptr{Cint}}()
    pindex = Ref{Ptr{Cint}}()
    tric   = Ref{Ptr{Cint}}()
    tris   = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessFace, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}),
                tess.obj, Cint(fIndex), len, xyz, uv, ptype, pindex, ntri, tris, tric))
    xyz  = ptr_to_array(xyz[] , len[],3)
    uv   = ptr_to_array(uv[]  , len[],2)
    tris = ptr_to_array(tris[],ntri[],3)
    tric = ptr_to_array(tric[],ntri[],3)
    return (xyz=xyz, uv=uv, ptype=ptype[], pindx=pindex[], tris=tris, tric=tric)
end


"""
Retrieves the data for the LOOPs associated with the discretization of a FACE
from a Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

fIndex:
    the FACE index (1 bias). The FACE Objects and number of FACEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.

Returns
-------
lIndices: the returned pointer to a vector of the last index
          for each LOOP (nloop in length).

Notes: (1) all boundary vertices are listed first for any FACE tessellation,

       (2) outer LOOP is ordered in the counter-clockwise direction, and

       (3) inner LOOP(s) are ordered in the clockwise direction.
"""
function getTessLoops(tess::Cego, fIndex::int)
    ptr1 = Ref{Cint}(0)
    ptr2 = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessLoops, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}),
                tess.obj, fIndex,ptr1,ptr2))
    return ptr_to_array(ptr2[], ptr1[])
end


"""
Creates Quadrilateral Patches for the indicated FACE and updates
the Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

qparms:
    a set of 3 parameters that drive the Quadrilateral patching for the
    FACE. Any may be set to zero to indicate the use of the default value:

    parms[1] EDGE matching tolerance expressed as the deviation
             from an aligned dot product [default: 0.05]

    parms[2] Maximum quad side ratio point count to allow
            [default: 3.0]

    parms[3] Number of smoothing loops [default: 0.0]

fIndex:
    the FACE index (1 bias)
"""
makeQuads!(tess::Cego, qparms::vNum, fIndex::int) = raiseStatus(ccall((:EG_makeQuads, C_egadslib), Cint, (ego, Ptr{Cdouble}, Cint),
                                                        tess.obj,Float32.(qparms), Cint(fIndex)))

"""
Returns a list of Face indices found in the Body-based Tessellation Object
that has been successfully Quadded via makeQuads.

Parameters
----------
tess:
    the input Body Tessellation Object

Returns
-------
fList: the returned pointer the Face indices (1 bias) - nface in length
       The Face Objects themselves can be retrieved via getBodyTopos
"""
function getTessQuads(tess::Cego)
    nface = Ref{Cint}(0)
    fList = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessQuads, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cint}}),
                tess.obj, nface,fList))
    arr   = ptr_to_array(fList[], nface[])
    free(fList)
    return arr
end


"""
Retrieves the data associated with the Quad-patching of a FACE
from a Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

fIndex:
    the FACE index (1 bias). The FACE Objects and number of FACEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo

Returns
-------
xyz: vertex coordinate data -> matrix[3][len]

uv: vertex asscaited parameter values -> matrix[2][len]

ptype: vertex type (-1 - internal, 0 - NODE, >0 EDGE)

pindx: vertex index (-1 internal)

npatch: number of patches
"""
function getQuads(tess::Cego, fIndex::int)
    len    = Ref{Cint}(0)
    xyz    = Ref{Ptr{Float64}}()
    uv     = Ref{Ptr{Float64}}()
    ptype  = Ref{Ptr{Cint}}()
    pindex = Ref{Ptr{Cint}}()
    npatch = Ref{Cint}(0)
    raiseStatus(ccall((:EG_getQuads, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}),
                 tess.obj, fIndex, len, xyz, uv, ptype, pindex, npatch))
    ptype  = ptr_to_array(ptype[],len[])
    pindex = ptr_to_array(pindex[],len[])
    xyz    = ptr_to_array(xyz[],len[],3)
    uv     = ptr_to_array(uv[],len[] ,2)
    return (xyz = xyz, uv = uv, ptype = ptype, pindx = pindex, npatch = npatch[])
end


"""
Retrieves the data associated with the Patch of a FACE
from a Body-based Tessellation Object.

Parameters
----------
tess:
    the input Body Tessellation Object

fIndex:
    the FACE index (1 bias). The FACE Objects and number of FACEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo

pIndex:
    the patch index (1-npatch from EG_getQuads)

Returns
-------
n1: patch size in the first direction (indexed by i)

n2: patch size in the second direction (indexed by j)

pvindex: list of n1*n2 indices that define the patch

pbounds: list of the neighbor bounding information for the patch
         (2*(n1-1)+2*(n2-1) in length). The first represents the segments at
         the base (j at base and increasing in i), the next is at the right (with i
         at max and j increasing). The third is the top (with j at max and i decreasing)
         and finally the left (i at min and j decreasing)
"""
function getPatch(tess::Cego, fIndex::int, pIndex::int)
    nu     = Ref{Cint}(0)
    nv     = Ref{Cint}(0)
    ipts   = Ref{Ptr{Cint}}(0)
    bounds = Ref{Ptr{Cint}}(0)
    raiseStatus(ccall((:EG_getPatch, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}),
                tess.obj, Cint(fIndex), Cint(pIndex), nu, nv, ipts, bounds))
    arr1   = ptr_to_array(  ipts[], nu[] * nv[])
    arr2   = ptr_to_array(bounds[], 2 * (nu[] - 1) + 2 * (nv[] -1))
    return (n1 = nu[], n2 = nv[], pvindex = arr1, pbounds = arr2)
end


"""
Moves the position of an EDGE vertex in a Body-based Tessellation Object.
Will invalidate the Quad patches on any FACEs touching the EDGE.

Parameters
----------
tess:
    the Body Tessellation Object (not on WIREBODIES)

eIndex:
    the EDGE index

vIndex:
    the Vertex index in the EDGE (2 - nVert-1)

t:
    the new parameter value on the EDGE for the point
"""
moveEdgeVert!(tess::Cego, eIndex::int, vIndex::int, t::Real) = raiseStatus(ccall((:EG_moveEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cdouble),
                                               tess.obj, Cint(eIndex), Cint(vIndex),Cdouble(t)))

"""
Deletes an EDGE vertex from a Body-based Tessellation Object.
Will invalidate the Quad patches on any FACEs touching the EDGE.

Parameters
----------
tess:
    the Body Tessellation Object (not on WIREBODIES)

eIndex:
    the EDGE index

vIndex:
    the Vertex index in the EDGE (2 - nVert-1)

dir:
    the direction to collapse any triangles (either -1 or 1)
"""
deleteEdgeVert!(tess::Cego, eIndex::int, vIndex::int, dir::int) = raiseStatus(ccall((:EG_deleteEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cint),
                                               tess.obj, Cint(eIndex), Cint(vIndex), Cint(dir)))

"""
Inserts vertices into the EDGE discretization of a Body Tessellation Object.
This will invalidate the Quad patches on any FACEs touching the EDGE.

Parameters
----------
tess:
    the Body Tessellation Object (not on WIREBODIES)

eIndex:
    the EDGE index

vIndex:
    the Vertex index in the EDGE (2 - nVert-1)

ts:
    the t values for the new points. Must be monotonically increasing and
    be greater than the t of vIndex and less than the t of vIndex+1.
"""
insertEdgeVerts!(tess::Cego, eIndex::int, vIndex::int, ts::Union{Real,vNum}) = raiseStatus(ccall((:EG_insertEdgeVerts, C_egadslib), Cint, (ego, Cint, Cint, Cint, Ptr{Cdouble}),
                                               tess.obj, eIndex, vIndex, length(vcat(ts)), collect(Cdouble,vcat([ts]...))))


"""
Opens an existing Tessellation Object for replacing EDGE/FACE discretizations.

Parameters
----------
tess:
    the Tessellation Object to open
"""
openTessBody(tess::Cego) = raiseStatus(ccall((:EG_openTessBody, C_egadslib), Cint, (ego,), tess.obj))


"""
Creates an empty (open) discretization object for a Topological BODY Object.

Parameters
----------
    body: body the input object, may be any Body type

Returns
-------
tess: resultant empty TESSELLATION object where each EDGE in the BODY must be
      filled via a call to EG_setTessEdge and each FACE must be filled with
      invocations of EG_setTessFace. The TESSSELLATION object is considered
      open until all EDGEs have been set (for a WIREBODY) or all FACEs have been
      set (for other Body types)
"""
function initTessBody(body::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_initTessBody, C_egadslib), Cint, (ego, Ptr{ego}), body.obj, ptr))
    return create_Cego(ptr[] ; ctxt =  body.ctxt, refs = body)
end


"""
Returns the status of a TESSELLATION Object.

Note: Placing the attribute ".mixed" on tess before invoking this
      function allows for tri/quad (2 tris) tessellations
      The type must be ATTRINT and the length is the number of FACEs,
      where the values are the number of quads

Parameters
----------
tess:
    the Tessellation Object to query

Returns
-------
body: the returned associated BODY Object

stat: the state of the tessellation: -1 - closed but warned, 0 - open, 1 - OK, 2 - displaced

npts: the number of global points in the tessellation (0 -- open)

icode: EGADS_SUCCESS -- complete, EGADS_OUTSIDE -- still open
"""
function statusTessBody(tess::Cego)
    ptr = Ref{ego}()
    np  = Ref{Cint}(0)
    st  = Ref{Cint}(0)
    raiseStatus(ccall((:EG_statusTessBody, C_egadslib), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}),
                tess.obj, ptr, st, np))
    return ( body = create_Cego(ptr[] ; ctxt =  tess.ctxt, delObj = false), stat = st[], npts = np[])
end


"""
Sets the data associated with the discretization of an EDGE
for an open Body-based Tessellation Object.

Notes: (1) all vertices must be specified in increasing t

       (2) the coordinates for the first and last vertex MUST match the
           appropriate NODE's coordinates

       (3) problems are reported to Standard Out regardless of the OutLevel

Parameters
----------
tess:
    the open Tessellation Object

eIndex:
    the EDGE index. The EDGE Objects and number of EDGEs
    can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
    If this EDGE already has assigned data, it is overwritten.

xyzs:
    the matrix [3][nPts] to the set of coordinate data

ts:
    the pointer to the parameter values associated with each vertex

"""
setTessEdge!(tess::Cego, eIndex::int, xyz::mNum, t::vNum) = raiseStatus(ccall((:EG_setTessEdge, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
                                      tess.obj, Cint(eIndex), Cint(length(t)), Float64.(vcat(xyz...)), Float64.(t)))


"""
Sets the data associated with the discretization of a FACE for an open Body-based Tessellation Object.

Parameters
----------
tess:
    the open Tessellation Object

fIndex:
    the FACE index (1 bias). The FACE Objects and number of FACEs can be
    retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
    If this FACE already has assigned data, it is overwritten.

xyz:
    the matrix to the set of coordinate data for all vertices [3][len]

uv:
    the matrix to the vertex parameter values [2][len]

tris:
    the matrix to triangle vertex indices [3][ntri]
"""
function setTessFace!(tess::Cego, fIndex::int, xyz::mNum, uv::mNum, tris::mInt)
    len1 = Int32(length(vec(uv))   / 2)
    len2 = Int32(length(vec(tris)) / 3)
    return raiseStatus(ccall((:EG_setTessFace, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}),
                        tess.obj, Cint(fIndex), len1, Float64.(vcat(xyz...)),Float64.(vcat(uv...)), len2, Cint.(vcat(tris...))))
end


"""
Completes the discretization for specified objects
for the input TESSELLATION object.

Note: an open TESSELLATION object is created by EG_initTessBody and
      can be partially filled via EG_setTessEdge and/or EG_setTessFace
      before this function is invoked.

Parameters
----------
tess:
    the Open Tessellation Object to close

parms:
    a set of 3 parameters that drive the EDGE discretization and the
    FACE triangulation. The first is the maximum length of an EDGE
    segment or triangle side (in physical space). A zero is flag that
    allows for any length. The second is a curvature-based value that
    looks locally at the deviation between the centroid of the discrete
    object and the underlying geometry. Any deviation larger than the
    input value will cause the tessellation to be enhanced in those
    regions. The third is the maximum interior dihedral angle (in degrees)
    between triangle facets (or Edge segment tangents for a WIREBODY
    tessellation), note that a zero ignores this phase.
"""
finishTess!(tess::Cego, params::vNum) = raiseStatus(ccall((:EG_finishTess, C_egadslib), Cint,
                                         (ego, Ptr{Cdouble}), tess.obj, dbl.(params)))


"""
Perform Local to Global index lookup. Tessellation Object must be closed.

Parameters
----------
tess:
    the closed Tessellation Object

ind:
    the topological index (1 bias) - 0 Node, (-) Edge, (+) Face

locl:
    the local (or Node) index

Returns
-------
gbl: the returned global index
"""
function localToGlobal(tess::Cego, ind::int, locl::int)
    ptr  = Ref{Cint}(0)
    raiseStatus(ccall((:EG_localToGlobal, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}),
                tess.obj, Cint(ind), Cint(locl), ptr));
    return ptr[]
end


"""
Returns the point type and index (like from EG_getTessFace)
with coordinates.

Parameters
----------
tess:
    the closed Tessellation Object

gbl:
    the global index

Returns
-------
ptype: the point type (-) Face local index, (0) Node, (+) Edge local index

pindex: the point topological index (1 bias)

xyz: the coordinates at this global index (can be NULL for no return)
"""
function getGlobal(tess::Cego, gbl::int)
    ptr1 = Ref{Cint}(0)
    ptr2 = Ref{Cint}(0)
    xyz  = zeros(3)
    raiseStatus(ccall((:EG_getGlobal, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
                tess.obj, gbl, ptr1, ptr2, xyz))
    return (ptype = ptr1[], pindex = ptr2[], xyz = xyz)
end


"""
Computes and returns the physical and inertial properties of a Tessellation Object.

Parameters
----------
tess:
    the Body Tessellation Object used to compute the Mass Properties
Returns
-------
volume,

area (length for EDGE, LOOP or WIREBODY),

CG: center of gravity (3),

IM: inertia matrix at CG (9)
"""
function tessMassProps(tess::Cego)
    props = zeros(14)
    raiseStatus(ccall((:EG_tessMassProps, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                tess.obj, props))
    return (volume = props[1], area = props[2], CG = props[3:5], IM = props[6:14])
end


"""
Performs the Boolean Operations on the source BODY Object(s).
The tool Object(s) are applied depending on the operation.
This supports Intersection, Splitter, Subtraction and Union.

Parameters
----------
src:
    the source in the form of a Model or Body Object

tool:
    the tool BODY object(s) in the form of a MODEL or BODY

oper:
    1-SUBTRACTION, 2-INTERSECTION, 3-FUSION, 4-SPLITTER

tol:
    the tolerance applied to perform the operation.
    If set to zero, the minimum tolerance is applied

Returns
-------
model: the resultant MODEL object (this is because there may be
       multiple bodies from either the subtraction or intersection operation).

Note: The BODY Object(s) for src and tool can be of any type but
      the results depend on the operation.
      Only works for OpenCASCADE 7.3.1 and higher.
"""
function generalBoolean(src::Cego, tool::Cego, oper::int ; tol::dbl=0.0)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_generalBoolean, C_egadslib), Cint, (ego, ego,Cint, Cdouble, Ptr{ego}),
                src.obj, tool.obj, Cint(oper), tol, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Fuses (unions) two SHEETBODYs resulting in a single SHEETBODY.

Parameters
----------
src:
    the source Sheet Body Object

tool:
   the tool SHEETBODY object

Returns
-------
sheet: the resultant SHEETBODY object
"""
function fuseSheets(src::Cego, tool::Cego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_fuseSheets, C_egadslib), Cint, (ego, ego, Ptr{ego}),
                src.obj, tool.obj, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Intersects the source BODY Object (that has the type
SOLIDBODY, SHEETBODY or FACEBODY) with a surface or
surfaces. The tool object contains the intersecting geometry
in the form of a FACEBODY, SHEETBODY, SOLIDBODY or a single
FACE.

Parameters
----------
src:
    the source Body Object

tool:
    the FACE/FACEBODY/SHEETBODY/SOLIDBODY tool object


Returns
-------
pairs: List of the FACE/EDGE object pairs -> matrix[2][nEdges]

model: The resultant MODEL object which contains the set of WIREBODY
       BODY objects (this is because there may be multiple LOOPS as a
       result of the operation).  The EDGE objects contained within
       the LOOPS have the attributes of the FACE in src responsible
       for that EDGE.

Note: The EDGE objects contained within the LOOPS have the attributes
      of the FACE in src responsible for that EDGE.
"""
function intersection(src::Cego, tool::Cego)
    nEdge = Ref{Cint}(0)
    pairs = Ref{Ptr{ego}}()
    model = Ref{ego}()
    raiseStatus(ccall((:EG_intersection, C_egadslib), Cint, (ego, ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{ego}),
                src.obj, tool.obj, nEdge, pairs, model))

    if nEdge[] == 0
        free(pairs)
        return (pairs = create_Cego(nothing), model = create_Cego(nothing))
    end
    arr = ptr_to_array(pairs[], nEdge[] * 2)
    free(pairs)
    return (pairs = reshape([create_Cego(arr[j]; ctxt =  src.ctxt, delObj = false) for j =1:length(arr)],2,Int64(nEdge[])),
            model = create_Cego(model[]; ctxt =  src.ctxt) )
end


"""
Imprints EDGE/LOOPs on the source BODY Object (that has the type
SOLIDBODY, SHEETBODY or FACEBODY). The EDGE/LOOPs are
paired with the FACEs in the source that will be scribed with the
EDGE/LOOP.

Parameters
---------
src:
    the source Body Object

pairs:
    list of FACE/EDGE and/or FACE/LOOP object pairs to scribe
    -> matrix[2][nObj] in len -- can be the output from intersect result

Returns
-------
result: the resultant BODY object (with the same type as the input source
        object, though the splitting of FACEBODY objects results in a
        SHEETBODY)
"""
function imprintBody(src::Cego, pairs::Matrix{Cego})
    _,ps = size(pairs)
    fpairs = vcat( getfield.(pairs,:obj)...)
    ptr    = Ref{ego}()
    raiseStatus(ccall((:EG_imprintBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
               src.obj,ps, fpairs, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Fillets the EDGEs on the source BODY Object (that has the type
SOLIDBODY or SHEETBODY).

Parameters
----------
src:
    the source Body Object

edges:
    list of EDGE objects to fillet

radius:
    the radius of the fillets created

Returns
-------
result: the resultant BODY object (with the same type as the input
        source object)

maps: list of Face mappings (in the result) which includes
      operations and an index to src where
      the Face originated -> matrix[2][nFaces]
"""
function filletBody(src::Cego, edges::Vector{Cego}, radius::Real)
    nedge = length(edges)
    res   = Ref{ego}()
    map   = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_filletBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}),
                src.obj, nedge, getfield.(edges,:obj),Cdouble(radius), res,map))
    Cres  = create_Cego(res[] ; ctxt =  src.ctxt)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(map[], length(bodyT), 2)
    free(map)
    return ( result = Cres, maps = arr)
end


"""
Chamfers the EDGEs on the source BODY Object (that has the type SOLIDBODY or SHEETBODY).

Parameters
----------
src:
    the source Body Object

edges:
    list of EDGE objects to chamfer - same len as faces

faces:
    list of FACE objects to measure dis1 from - same len as edges

dis1:
    the distance from the FACE object to chamfer

dis2:
    the distance from the other FACE to chamfer


Returns
-------
result: the resultant BODY object (with the same type as the input source object)

maps: list of Face mappings (in the result) which includes operations and an
      index to src where the Face originated -> matrix[2][nFaces]
"""
function chamferBody(src::Cego, edges::Vector{Cego}, faces::Vector{Cego}, dis1::dbl, dis2::dbl)
    n     = length(edges)
    n    != length(faces) && (throw(ErrorMessage("edges and faces must be same length")))
    fmap  = Ref{Ptr{Cint}}()
    res   = Ref{ego}()
    raiseStatus(ccall((:EG_chamferBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}, Cdouble, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}),
                src.obj, n, getproperty.(edges,:obj),getproperty.(faces,:obj), dis1, dis2, res, fmap))
    Cres  = create_Cego(res[]; ctxt =  src.ctxt)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(fmap[], length(bodyT),2)
    free(fmap)
    return (result = Cres, maps = arr)
end


"""
A hollowed solid is built from an initial SOLIDBODY Object and a set of
FACEs that initially bound the solid. These FACEs are removed and the
remaining FACEs become the walls of the hollowed solid with the
specified thickness. If there are no FACEs specified then the Body
is offset by the specified distance (which can be negative).

Parameters
----------
src:
    the source Body Object (SOLIDBODY, SHEETBODY or FACEBODY)

faces:
    face or list of faces to remove - (nothing performs an Offset)

off:
    the wall thickness (offset) of the hollowed result

join:
    0 - fillet-like corners, 1 - expanded corners

Returns
-------
result: the resultant BODY object

maps: list of Face mappings (in the result) which includes operations and an
      index to src where the Face originated -> matrix[2][nFaces]
"""
function hollowBody(src::Cego, faces::Union{Cego,Vector{Cego}, Nothing}, offset::Real, join::int)
    Cface = emptyVar(faces) ? C_NULL : getfield.(vcat([faces]...),:obj)
    nface = Cface == C_NULL ? 0 : length(Cface)
    fmap  = Ref{Ptr{Cint}}()
    res   = Ref{ego}()
    raiseStatus(ccall((:EG_hollowBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}, Ptr{Ptr{Cint}}),
                src.obj, nface, Cface, Cdouble(offset), Cint(join), res, fmap))
    Cres  = create_Cego(res[]; ctxt =  src.ctxt)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(fmap[], length(bodyT), 2)
    free(fmap)
    return (result = Cres, maps = arr)
end


"""
Rotates the source Object about the axis through the angle
specified. If the Object is either a LOOP or WIREBODY the
result is a SHEETBODY. If the source is either a FACE or
FACEBODY then the returned Object is a SOLIDBODY.

Parameters
----------
src:
    the source Body Object (not SHEETBODY or SOLIDBODY)

angle:
    the angle to rotate the object through [0-360 Degrees]

axis:
      [2][3] or [6]: a point (on the axis)[3] and a direction [3]

Returns
-------
result: the resultant BODY object (type is one greater than the input
        source object)
"""
function rotate(src::Cego, angle::Real, axis::Union{vNum, mNum})
    ptr = Ref{ego}()
    ax = Cdouble.(vcat(axis...))

    length(ax) != 6 && throw(ErrorMessage(" rotate: flat axis length != 6 "))
    raiseStatus(ccall((:EG_rotate, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}),
                src.obj, Cdouble(angle), ax, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Extrudes the source Object through the distance specified.
If the Object is either a LOOP or WIREBODY the result is a
SHEETBODY. If the source is either a FACE or FACEBODY then
the returned Object is a SOLIDBODY.

Parameters
----------
src:
    the source Body Object (not SHEETBODY or SOLIDBODY)

dist:
    the distance to extrude

dir:
    the vector that is the extrude direction (3 in length)

Returns
-------
result: the resultant BODY object (type is one greater than the
        input source object)
"""
function extrude(src::Cego, dist::Real, dir::vNum)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_extrude, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}),
                  src.obj, Cdouble(dist), Cdouble.(dir), ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Sweeps the source Object through the "spine" specified.
The spine can be either an EDGE, LOOP or WIREBODY.
If the source Object is either a LOOP or WIREBODY the result is a SHEETBODY.
If the source is either a FACE or FACEBODY then the returned Object is a SOLIDBODY.

Note: This does not always work as expected...

Parameters
----------
src:
    the source Body Object (not SHEETBODY or SOLIDBODY)

spine:
    the Object used as guide curve segment(s) to sweep the source through

mode:
    0 - CorrectedFrenet      1 - Fixed

    2 - Frenet               3 - ConstantNormal

    4 - Darboux              5 - GuideAC

    6 - GuidePlan            7 - GuideACWithContact

    8 - GuidePlanWithContact 9 - DiscreteTrihedron

Returns
-------
result: the returned resultant Body Object (type is one higher than src)
"""
function sweep(src::Cego, spine::Cego, mode::int)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_sweep, C_egadslib), Cint, (ego, ego, Cint, Ptr{ego}),
                src.obj, spine.obj, mode, ptr))
    return create_Cego(ptr[] ; ctxt =  src.ctxt)
end


"""
Produces a BODY Object (that has the type SOLIDBODY or SHEETBODY)
that goes through the sections by ruled surfaces between each. All
sections must have the same number of Edges (except for NODEs) and
the Edge order in each is used to specify the connectivity.

Parameters
----------
sections:
    ist of FACEBODY, FACE, WIREBODY or LOOP objects to Blend
    FACEBODY/FACEs must have only a single LOOP;
    the first and last sections can be NODEs;
    if the first and last are NODEs and/or FACEBODY/FACEs
    (and the intermediate sections are CLOSED) the result
    will be a SOLIDBODY otherwise a SHEETBODY will be constructed;
    if the first and last sections contain equivalent LOOPS
    (and all sections are CLOSED) a periodic SOLIDBODY is created

Returns
-------
result: the resultant BODY object
"""
function ruled(secs::Vector{Cego})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_ruled, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{ego}),
                Cint(length(secs)), vcat(getfield.(secs,:obj)...), ptr))
    return create_Cego(ptr[] ; ctxt =  secs[1].ctxt)
end


"""
Simply lofts the input Objects to create a BODY Object
(that has the type SOLIDBODY or SHEETBODY). Cubic BSplines are used.
All sections must have the same number of Edges (except for NODEs)
and the Edge order in each (defined in a CCW manner) is used to
specify the loft connectivity.

Parameters
----------
sections:
    list of FACEBODY, FACE, WIREBODY or LOOP objects to Blend
    FACEBODY/FACEs must have only a single LOOP;
    the first and last sections can be NODEs;
    if the first and last are NODEs and/or FACEBODY/FACEs
    (and the intermediate sections are CLOSED) the result
    will be a SOLIDBODY otherwise a SHEETBODY will be constructed;
    if the first and last sections contain equivalent LOOPS
    (and all sections are CLOSED) a periodic SOLIDBODY is created

Optional
--------
rc1:
    specifies treatment* at the first section (or Nothing for no treatment)

rc2:
    specifies treatment* at the first section (or Nothing for no treatment)

* for NODEs -- elliptical treatment (8 in length): radius of
curvature1, unit direction, rc2, orthogonal direction;
nSection must be at least 3 (or 4 for treatments at both ends)
for other sections -- setting tangency (4 in length):
magnitude, unit direction for FACEs with 2 or 3 EDGEs -- make
a Wing Tip-like cap: zero, growthFactor (len of 2)

Returns
-------
result: the resultant BODY object
"""
function blend(secs::Vector{Cego} ;
               rc1=nothing::Union{Nothing,vNum},rc2=nothing::Union{Nothing,vNum})
    ptr = Ref{ego}()

    prc1 = rc1 == nothing ? C_NULL : Cdouble.(rc1[1:min(length(rc1),8)])
    prc2 = rc2 == nothing ? C_NULL : Cdouble.(rc2[1:min(length(rc2),8)])

    raiseStatus(ccall((:EG_blend, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ego}),
                Cint(length(secs)), getfield.(secs,:obj),prc1, prc2, ptr))
    return create_Cego(ptr[] ; ctxt =  secs[1].ctxt)
end


"""
Takes as input a Body Tessellation Object and returns an Open EBody fully initialized with Effective}
Objects (that may only contain a single non-effective object). EEdges are automatically merged where
possible (removing Nodes that touch 2 Edges, unless degenerate or marked as ".Keep"). The tessellation
is used (and required) in order that single UV space be constructed for EFaces.

Parameters
----------
tess:
    the input Solid or Sheet Body Tessellation Object (can be quite coarse)

angle:
    angle used to determine if Nodes on open
    Edges of Sheet Bodies can be removed. The dot of the tangents at the Node is
    compared to this angle (in degrees). If the dot is greater than the angle,
    the Node does not disappear. The angle is also used to test Edges to see if
    they can be removed. Edges with normals on both trimmed Faces showing deviations
    greater than the input angle will not disappear.

Returns
-------
ebody: the resultant Open Effective Topology Body Object
"""
function initEBody(tess::Cego, angle::Real)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_initEBody, C_egadslib), Cint, (ego, Cdouble, Ptr{ego}),
                tess.obj, Cdouble(angle), ptr))
    return create_Cego(ptr[] ; ctxt =  tess.ctxt, refs = tess)
end


"""
Modifies the EBody by finding "free" Faces (a single Face in an EFace)
with attrName and the same attribute value(s), thus making a collection of EFaces.
The single attrName is placed on the EFace(s) unless in "Full Attribute" mode, which
then performs attribute merging. This function returns the number of EFaces
possibly constructed (neface). The UVbox can be retrieved via calls to either
getTopology or getRange with the returned appropriate efaces.

Parameters
----------
ebody:
    the input Open Effective Topology Body Object

attrName:
    the attribute name used to collect Faces into an EFaces. The attribute value(s)
    are then matched to make the collections.

Returns
-------
efaces: the list of EFaces constructed
"""
function makeAttrEFaces(EBody::Cego, attrName::str)
    neface = Ref{Cint}(0)
    efaces = Ref{Ptr{ego}}()
    raiseStatus(ccall((:EG_makeAttrEFaces, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}, Ptr{Ptr{ego}}),
                EBody.obj, attrName, neface, efaces))
    arr    = ptr_to_array(efaces[], neface[])
    free(efaces)
    return [create_Cego(arr[j] ; ctxt =  EBody.ctxt, delObj = false) for j = 1: length(arr)]
end


"""
Modifies the EBody by removing the nface "free" Faces and replacing them
with a single eface (returned for convenience, and note that ebody has been updated).
There are no attributes set on eface unless the "Full Attribution" mode is in place.
This constructs a single UV for the faces. The UVbox can be retrieved via calls to either
getTopology or getRange with eface.

Parameters
----------
ebody:
    the input Open Effective Topology Body Object

faces:
    the list of Face Objects used to make eface

Returns
-------
eface: the EFace constructed
"""
function makeEFace(EBody::Cego, faces::Union{Vector{Cego},Cego})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeEFace, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
                 EBody.obj, length(faces), getfield.(vcat([faces]...),:obj), ptr))
    return create_Cego(ptr[] ; ctxt =  EBody.ctxt, delObj = false)
end


"""
Finish an EBody
"""
function finishEBody!(EBody::Cego)
    raiseStatus(ccall((:EG_finishEBody, C_egadslib), Cint, (ego,), EBody.obj))
    EBody.refs == nothing && throw(ErrorMessage(" finishEBody, EBody doesn't have any ref "))
    aux = getfield.(vcat([EBody.refs]...),:refs)
    EBody.refs = vcat(aux...)
end

"""
Returns the evaluated location in the BRep for the Effective Topology entity.

Parameters
----------
eobj:
    the input Effective Topology Object (either EEdge or EFace)

eparam:
    t for EEdges / the [u, v] for an EFace

Returns
-------
obj: the returned source Object in the Body

param: the returned t for an Edge / the returned [u, v] for the Face
"""
function effectiveMap(eobject::Cego, eparam::Union{vNum,Real})
    obj     = Ref{ego}()
    vparam  = Cdouble.(vcat([eparam]...))
    infoOBJ = unsafe_load(eobject.obj)
    if infoOBJ.oclass == EEDGE
        length(vparam) != 1 &&
            (throw(ErrorMessage(" effectiveMap for EEDGE shoudl have a single eparam $vparam")))
    elseif infoOBJ.oclass == EFACE
        length(vparam) != 2 &&
            (throw(ErrorMessage(" effectiveMap for EFACE shoudl have 2 eparams $vparam")))
    else
        throw(ErrorMessage(" effectiveMap oclass shoudl be either EFACE or EEDGE  and not $infoOBJ.oclass"))
    end
    raiseStatus(ccall((:EG_effectiveMap, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}, Ptr{Cdouble}),
                eobject.obj, vparam, obj, Cdouble.(vparam)))
    return (obj = create_Cego(obj[], ctxt =  eobject.ctxt, delObj = false), param = length(vparam) ==2 ? vparam : vparam[1])
end


"""
Returns the list of Edge entities in the source Body that make up the EEdge.
A pointer to an integer list of senses for each Edge is returned as well
as the starting t value in the EEdge (remember that the t will go in the
opposite direction in the Edge if the sense is SREVERSE).

Parameters
----------
eedge:
    the Effective Topology Edge Object

Returns
-------
edges: list of Edges - nedge in length

senses: list of senses - nedge in length

tstart: list of t starting values - nedge in length
"""
function effectiveEdgeList(eedge::Cego)
    nedge  = Ref{Cint}(0)
    edges  = Ref{Ptr{ego}}()
    senses = Ref{Ptr{Cint}}()
    tstart = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_effectiveEdgeList, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}),
                eedge.obj, nedge, edges, senses, tstart))
    arr_e = ptr_to_array(edges[], nedge[])
    arr_s = ptr_to_array(senses[], nedge[])
    arr_t = ptr_to_array(tstart[], nedge[])
    free(edges); free(senses); free(tstart)
    return (edges = arr_e, senses = arr_s, tstart = arr_t)
end

# Access all elements from a C pointer and convert it to ta Julia array
# If dQ != 0, matrix is dim (dQ , dP)
function ptr_to_array(Ptr, dP, dQ=1)
    Ptr == C_NULL && return C_NULL
    arr = [unsafe_load(Ptr, dQ * i + j) for j in 1:dQ, i in 0:dP-1]
    ( (dQ == 1 || (dQ > 1 && dP == 1)) && (arr = vcat(arr...)))
    return arr
end

function raiseStatus(icode)
    if length(egadsERROR) == 0
        Base.open(joinpath(EGADS_WRAP_PATH,FILE_EGADS_CTT)) do io
            for line in eachline(io, keep=true) # keep so the new line isn't chomped
                if occursin("const EGADS_",line)
                    s_val1 = findfirst("(",line)
                    s_val2 = findfirst(")",line)
                    s_ctt = findfirst(" ", line)
                    s_eq  = findfirst(" =", line)
                    s_nam = String(line[s_ctt[1]+1:s_eq[1]-1])
                    aux   = egadsERR_NAMES(parse(Int32,line[s_val1[1]+1:s_val2[1]-1] ), s_nam)
                    push!(egadsERROR,aux)
                    continue
                end
            end
        end
    end
    try @assert icode == EGADS_SUCCESS
    catch
        idx = findall(x-> x == icode, getfield.(egadsERROR,:var))
        msg = "EGADS error: "*egadsERROR[idx[1]].name
        throw(ErrorMessage(msg))
    end
end
#ExportAll.@exportAll
end
