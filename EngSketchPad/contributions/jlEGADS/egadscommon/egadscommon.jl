#=
*
*              jlEGADS --- Julia version of EGADS API
*
*
*      Copyright 2011-2025, Massachusetts Institute of Technology
*      Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
*
=#
__precompile__()

# EGADS env set
if !haskey(ENV, "ESP_ROOT")
    throw(@error " ESP_ROOT must be set -- Please fix the environment...\n")
end

# EGADS paths + constants
ESP_ROOT=ENV["ESP_ROOT"] |> normpath
try @assert isdir(ESP_ROOT)
catch
    throw(@error " ESP_ROOT=$ESP_ROOT is not a directory -- Please fix the environment...\n")
end

const ESP_ROOT_LIB = joinpath(ESP_ROOT, "lib")     |> normpath
@info " ESP_ROOT_LIB $ESP_ROOT_LIB"
try @assert isdir(ESP_ROOT_LIB)
catch
    throw(@error " $ESP_ROOT_LIB is not a directory -- Please fix the environment...\n")
end

const EGADS_COMMON_PATH = normpath(@__DIR__)
const FILE_EGADS_CTT    = joinpath(EGADS_COMMON_PATH,"egads_types_ctts.jl")

# Set egads library name. WARNING ! If you change C_egadslib name, all calls inside egads will fail.
const LIBEGADS        = "C_egadslib"

try @assert isfile(FILE_EGADS_CTT )
catch
    throw(@error "egads_types_ctts.jl file not found in egadscommon folder. Run 'make' inside the egadscommon folder ")
end

include(FILE_EGADS_CTT)
# SET GLOBAL PATHS TO DEPS. CHECK ESP_ROOT and include files exist !!

#egads library
global C_egadslib = "none"


struct egadsERR_NAMES
    var  ::Cint
    name ::String
end

global egadsERROR = Vector{egadsERR_NAMES}()

Base.open(FILE_EGADS_CTT) do io
    for line in eachline(io, keep=true) # keep so the new line isn't chomped
        if occursin("const EGADS_",line)
            s_val1 = findfirst("(",line)
            s_val2 = findfirst(")",line)
            s_ctt = findfirst(" ", line)
            s_eq  = findfirst(" =", line)
            s_nam = String(line[s_ctt[1]+1:s_eq[1]-1])
            aux   = egadsERR_NAMES(parse(Cint,line[s_val1[1]+1:s_val2[1]-1] ), s_nam)
            push!(egadsERROR,aux)
        end
     end
end

# MY DATA TYPES
const str  = Union{String,Char}
const int  = Union{Cint,Int64}
const dbl  = Float64
const vInt = Union{Vector{Int64},Vector{Cint}}
const vDbl = Vector{Float64}
const vStr = Union{Vector{Char},Vector{String}}
const mInt = Union{Matrix{Int64},Matrix{Cint}}
const mDbl = Matrix{Float64}
const mStr = Union{Matrix{String}, Matrix{Char}}

# Check if a variable || vector is empty
emptyVar(a) = all(x -> x ∈ [C_NULL,nothing,Nothing], vcat([a]...))

mutable struct CSys
    val   :: vDbl
    ortho :: mDbl
    function CSys(val::Union{vDbl,vInt}, ortho::Union{Nothing, mDbl}=nothing)
        if ortho == nothing
            ortho = [[0.0, 0.0, 0.0] [1.0, 0.0, 0.0] [0.0, 1.0, 0.0] [0.0, 0.0, 1.0]]
        end
        x = new(val, ortho)
    end
end

"""
Wraps an egads CONTXT to trace memory

Members
----------
Context:
    an egads ego

delObj:
    If True then this ego takes ownership of obj and calls EG_close during garbage collection

"""
mutable struct Context
    ego     :: ego
    _delObj :: Bool
    # Julia garbage colletion will call finilize in random order.
    # So we will use our own reference counting to delete egos
    # in reverse order of creation
    _finalize_called :: Bool
    _refCount        :: Int
    function Context(; obj=ego()::ego, delObj::Bool = false)
        if obj == C_NULL
            ptr = Ref{ego}()
            raiseStatus(ccall((:EG_open, C_egadslib), Cint, (Ptr{ego},), ptr))
            obj = ptr[] ; delObj = true
        else
            info = unsafe_load(obj)
            info.oclass != CONTXT && throw(ErrorException("struct Context oclass needs to be CONTXT !!"))
        end
        x = new(obj,delObj,false,0)
        finalizer(x) do t
            #@async println(" CONTEXT FINALIZE $(t.obj) !!!")
            t._finalize_called = true
            if t._refCount == 0 && t._delObj
                #@async println(" CONTEXT CLOSE $(t.obj) !!!")
                _close!(t)
            end
        end
        return x
    end
end


"""
Wraps an egads ego to trace memory.

Members
----------
obj:
    an egads ego that is not context: CONTEXT has its own struct: Context

Context:
    the Context associated with the ego

refs:
    list of ego's that are used internaly in obj. Used to ensure proper deletion order

delObj:
    If True then this ego takes ownership of obj and calls EG_deleteObject during garbage collection

"""
mutable struct Ego
    _ego    :: ego
    ctxt    :: Context
    _refs   :: Vector{Ego}
    _delObj :: Bool
    # Julia garbage colletion will call finilize in random order.
    # So we will use our own reference counting to delete egos
    # in reverse order of creation
    _finalize_called :: Bool
    _refCount        :: Int
    function Ego(obj::ego; ctxt=nothing::Union{Nothing,Context},
                            refs=nothing::Union{Nothing,Vector{Ego},Ego},
                            delObj=false::Bool)

        # Check for an invalid C_NULL pointer
        obj == C_NULL && throw(ErrorException("struct Ego obj may not be C_NULL !!"))

        # Load the context of the ego if nothing provided
        if ctxt === nothing
            context = Ref{ego}()
            ccall((:EG_getContext, C_egadslib), Cint, (ego, Ptr{ego}), obj, context)
            ctxt = Context(context[], delObj = false)
        end

        # Convert refs to Vector{Ego} if needed
        # and increment reference counters
        if refs === nothing
            refs = Vector{Ego}()
        else
            refs = vcat([refs]...)
        end

        # Create the new object and increment the reference counters
        x = new(obj, ctxt, refs, delObj, false, 0)
        for ref in x._refs
            ref._refCount += 1
        end
        x.ctxt._refCount += 1

        finalizer(x) do t
            #@async println("Finalizing  : $(t._ego).\n")
            # Decrement the reference count
            for ref in t._refs
                ref._refCount -= 1
            end
            t.ctxt._refCount -= 1

            # Mark that finlize for t has been called
            t._finalize_called = true

            # Delete this object
            if t._refCount == 0 && t._delObj
                #@async println("deleteObject!  : $(t._ego).\n")
                _deleteObject!(t)
            end

            # Delete any objects whos finalizer was called prematurely
            for ref in t._refs
                if ref._finalize_called && ref._refCount == 0 && ref._delObj
                    #@async println("ref deleteObject!  : $(ref._ego).\n")
                    _deleteObject!(ref)
                end
            end
            if t.ctxt._finalize_called && t.ctxt._refCount == 0 && t.ctxt._delObj
                #@async println("ref _close!  : $(t.ctxt._ego).\n")
                _close!(t.ctxt)
            end
        end
        return x
    end
end


###
# Only for internal usage
#
# Free memory in pointer. Called whenever a EG_() function has freeable arguments
function _free(pointer)
    # Must use the EG_free for Windows
    ccall((:EG_free, C_egadslib), Cvoid, (Ptr{Cvoid},), pointer[])
end

#
# deleteObject! and _close! should never be called externally
function _deleteObject!(x::Ego)
    icode = ccall((:EG_deleteObject, C_egadslib), Cint, (ego,), x._ego)
    x._ego = C_NULL
    (icode < 0) && raiseStatus(icode)
    return icode
end

function _close!(obj::Context)
    # Check that garbage collection is working properly
    oclass = Ref{Cint}(0)
    mtype  = Ref{Cint}(0)
    topObj = Ref{ego}()
    prev   = Ref{ego}()
    next   = Ref{ego}()
    raiseStatus(ccall((:EG_getInfo, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}),
                    obj.ego, oclass, mtype, topObj, prev, next))
    if "$(@__MODULE__)" == "egadslite"
        while next[] != C_NULL
            raiseStatus(ccall((:EG_getInfo, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}),
            next[], oclass, mtype, topObj, prev, next))
            (oclass == TESSELLATION) && @warn "Please report a garbage collection bug!"
        end
    else
        (next[] != C_NULL) && @warn "Please report a garbage collection bug!"
    end
    icode = ccall((:EG_close, C_egadslib), Cint, (ego,), obj.ego)
    obj.ego = C_NULL
    (icode < 0) && raiseStatus(icode)
    return icode
end
###

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
function loadModel(context::Context, bitFlag::int, name::str)
    ptr = Ref{ego}()
    raiseStatus((ccall((:EG_loadModel, C_egadslib), Cint, (ego, Cint, Cstring, Ptr{ego}),
                    context.ego, Cint(bitFlag), Base.cconvert(Cstring,name), ptr)))
    return Ego(ptr[] ; ctxt =   context, delObj = true)
end



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
function saveModel!(model::Ego, name::str, overwrite::Bool = false)
    overwrite && rm(name, force = true)
    raiseStatus(ccall((:EG_saveModel, C_egadslib), Cint, (ego, Cstring), model._ego, name))
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
function getInfo(object::Union{Ego, Context})
    oclass = Ref{Cint}(0)
    mtype  = Ref{Cint}(0)
    topObj = Ref{ego}()
    prev   = Ref{ego}()
    next   = Ref{ego}()
    obj    = typeof(object) == Ego ? object._ego : object.ego
    raiseStatus(ccall((:EG_getInfo, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}),
                    obj, oclass, mtype, topObj, prev, next))
    context = typeof(object) == Ego ? object.ctxt : object
    return (oclass=oclass[], mtype=mtype[],
            topObj= topObj[] == C_NULL ? nothing : Ego(topObj[];ctxt = context, delObj = false),
            prev  = prev[]   == C_NULL ? nothing : Ego(prev[];  ctxt = context, delObj = false),
            next  = next[]   == C_NULL ? nothing : Ego(next[];  ctxt = context, delObj = false))
end


"""
Sets the EGADS verbosity level, the default is 1.

Parameters
----------
outlevel:

    0 <= outlevel <= 3

    0-silent to 3-debug
"""
setOutLevel(context::Context, outLevel::int) = ccall((:EG_setOutLevel, C_egadslib), Cint, (ego, Cint), context.ego, Cint(outLevel))


"""
Returns the number of attributes found with this object

Parameters
----------
object:
    the Object
"""
function attributeNum(object::Ego)
    ptr = Ref{Cint}(0)
    raiseStatus(ccall((:EG_attributeNum, C_egadslib), Cint, (ego, Ptr{Cint}), object._ego, ptr))
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
function attributeRet(object::Ego, name::str)
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    icode = ccall((:EG_attributeRet, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
             object._ego, name, atype, len, ints, reals, str)
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
function attributeGet(object::Ego, index::int)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    name  = Ref{Cstring}()
    raiseStatus(ccall((:EG_attributeGet, C_egadslib), Cint, (ego, Cint, Ptr{Cstring}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
                object._ego, Cint(index), name, atype, len, ints, reals, str))

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
function attributeNumSeq(object::Ego, name::str)
    ptr = Ref{Cint}(0)
    raiseStatus(ccall((:EG_attributeNumSeq, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}), object._ego, name, ptr))
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
function attributeRetSeq(object::Ego, name::str,index::int)
    atype = Ref{Cint}(0)
    len   = Ref{Cint}(0)
    ints  = Ref{Ptr{Cint}}()
    reals = Ref{Ptr{Cdouble}}()
    str   = Ref{Cstring}()
    raiseStatus(ccall((:EG_attributeRetSeq, C_egadslib), Cint, (ego, Cstring, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Cstring}),
                object._ego, name, Cint(index), atype, len, ints, reals, str))
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
function getGeometry(object::Ego)
    oclass  = Ref{Cint}(0)
    mtype   = Ref{Cint}(0)
    refGeom = Ref{ego}()
    rvec    = Ref{Ptr{Cdouble}}()
    ivec    = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getGeometry, C_egadslib), Cint,(ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}),
              object._ego, oclass, mtype, refGeom, ivec, rvec))

    nivec = Ref{Cint}(0)
    nrvec = Ref{Cint}(0)
    ccall((:EG_getGeometryLen, C_egadslib), Cvoid, (ego, Ptr{Cint}, Ptr{Cint}),object._ego, nivec, nrvec)
    li = nivec[]
    lr = nrvec[]

    j_ivec = ((ivec[] != C_NULL && li > 0) ? ptr_to_array(ivec[], li) : nothing)

    j_rvec = ((rvec[] != C_NULL && lr > 0) ? ptr_to_array(rvec[], lr) : nothing)
    _free(ivec) ; _free(rvec)

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
function getRange(object::Ego)
    geom = getInfo(object)
    r = ( (geom.oclass == PCURVE || geom.oclass == CURVE || geom.oclass == EDGE
        || geom.oclass == EEDGE) ? zeros(2) : zeros(4))
    p = Ref{Cint}(0)
    raiseStatus(ccall((:EG_getRange, C_egadslib), Cint,(ego, Ptr{Cdouble}, Ptr{Cint}),
                object._ego, r,p))
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
function arcLength(object::Ego, t1::dbl, t2::dbl)
    alen = Ref{Cdouble}(0.0)
    raiseStatus(ccall((:EG_arcLength, C_egadslib), Cint, (ego, Cdouble, Cdouble, Ptr{Cdouble}),
                object._ego, t1, t2, alen))
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
function evaluate(object::Ego, params::Union{dbl,vDbl, Nothing}=nothing)

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
    param = params == nothing ? C_NULL : vcat(params...)
    raiseStatus(ccall((:EG_evaluate, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object._ego, param, xyz))
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
function invEvaluate(object::Ego, pos::vDbl)
    r = zeros(3)
    p = zeros(2)
    raiseStatus(ccall((:EG_invEvaluate, C_egadslib), Cint,(ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}),
                object._ego,pos, p,r))
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
function curvature(object::Ego,params::Union{dbl,vDbl})
    param  = vcat(params...) #allows feeding pars as a double (not vector)
    geoPtr = unsafe_load(object._ego)
    if geoPtr.oclass == PCURVE
        res = zeros(3)
    elseif (geoPtr.oclass == EDGE  || geoPtr.oclass == EEDGE ||
            geoPtr.oclass == CURVE)
        res = zeros(4)
    elseif (geoPtr.oclass == FACE || geoPtr.oclass == EFACE ||
            geoPtr.oclass == SURFACE)
        res = zeros(8)
    else
        @info "curvature object oclass $geoPtr.oclass not recognised "
        raiseStatus(EGADS_GEOMERR)
    end
    raiseStatus(ccall((:EG_curvature, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object._ego, param, res))
    if length(res) != 8
        return (curvature = res[1], direction = res[2:end])
    else
        return (curvature = res[[1,5]], direction = [res[2:4] res[6:8]]' )
    end
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
function getTopology(topo::Ego)
    geom   = Ref{ego}()
    oclass = Ref{Cint}(0)
    mtype  = Ref{Cint}(0)
    lims   = zeros(4)
    nChild = Ref{Cint}(0)
    sen    = Ref{Ptr{Cint}}(0)
    child  = Ref{Ptr{ego}}()

    raiseStatus(ccall((:EG_getTopology, C_egadslib), Cint, (ego, Ptr{ego},
                 Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}),
                 topo._ego, geom, oclass, mtype,lims, nChild, child, sen))

    nC = (oclass[] == LOOP  && geom[] != C_NULL) ? nChild[] * 2 :
         (oclass[] == MODEL && mtype[] > 0     ) ? mtype[] : nChild[]

    aux = ptr_to_array(child[], nC)

    children = [Ego(unsafe_load(child[], j) ; ctxt =   topo.ctxt, delObj = false) for j = 1:nC]

    senses   = ptr_to_array(sen[], nChild[])

    reals =  oclass[] == NODE ? lims[1:3] : oclass[] == EDGE ? lims[1:2] :
            (oclass[] == FACE || oclass[] == EFACE) ? [lims[[1,3]] lims[[2,4]]] : nothing

    return (geom = geom[] == C_NULL ? nothing : Ego(geom[] ; ctxt =   topo.ctxt, delObj = false, refs = topo),
            oclass = oclass[], mtype=mtype[], reals=reals, senses = senses, children = children)
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
function getBodyTopos(body::Ego, oclass::int; ref=nothing::Union{Ego,Nothing})

    ntopo = Ref{Cint}(0)
    topos = Ref{Ptr{ego}}()
    #ref   = typeof(ref) == Ego ? ref._ego : C_NULL
    sref  = ref == nothing ? C_NULL : ref._ego
    raiseStatus( ccall((:EG_getBodyTopos, C_egadslib), Cint,
                (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{ego}}),
                 body._ego, sref, Cint(oclass), ntopo, topos))

    topo_arr = ptr_to_array(topos[], ntopo[])
    _free(topos)
    return [Ego(topo_arr[j]; ctxt = body.ctxt, delObj = false, refs=body) for j = 1:ntopo[] ]
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
indexBodyTopo(body::Ego,object::Ego) = ccall((:EG_indexBodyTopo, C_egadslib), Cint, (ego, ego), body._ego, object._ego)


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
function objectBodyTopo(body::Ego, oclass::int, index::int)
    ptr = Ref{ego}()
    loadBody = unsafe_load(body._ego)
    raiseStatus(ccall((:EG_objectBodyTopo, C_egadslib), Cint, (ego, Cint, Cint, Ptr{ego}),
                body._ego, Cint(oclass), Cint(index), ptr))
    return Ego(ptr[] ; ctxt =   body.ctxt, delObj = false, refs=body)
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
function getBoundingBox(topo::Ego)
    bbox = zeros(6)
    raiseStatus(ccall((:EG_getBoundingBox, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo._ego, bbox))
    return reshape(bbox, 2,3)
end


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
function inTopology(object::Ego, xyz::vDbl)
    s = ccall((:EG_inTopology, C_egadslib), Cint, (ego, Ptr{Cdouble}), object._ego, xyz)
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
function inFace(face::Ego, uv::vDbl)
    s = ccall((:EG_inFace, C_egadslib), Cint, (ego, Ptr{Cdouble}), face._ego, uv)
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
function getEdgeUV(face::Ego, edge::Ego, sense::int, t::dbl)
    uv = zeros(2)
    raiseStatus(ccall((:EG_getEdgeUV, C_egadslib), Cint, (ego, ego, Cint, Cdouble, Ptr{Cdouble}),
                face._ego, edge._ego, Cint.(sense), t, uv))
    return uv
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
function getBody(object::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_getBody, C_egadslib), Cint, (ego, Ptr{ego}),object._ego,ptr))
    return Ego(ptr[] ; ctxt =   object.ctxt, delObj = false)
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
function tolerance(topo::Ego)
    tol = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_tolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo._ego, tol))
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
function getTolerance(topo::Ego)
    tol = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_getTolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                topo._ego, tol))
    return tol[]
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
function getWindingAngle(edge::Ego, t::dbl)
    ptr = Ref{Cdouble}(0)
    raiseStatus(ccall((:EG_getWindingAngle, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}),
                edge._ego, t, ptr))
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
function makeTessGeom(geom::Ego, limits::vDbl, sizes::vInt)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeTessGeom, C_egadslib), Cint,
                (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{ego}),
                geom._ego,vcat(limits...), Cint.(sizes), ptr))
    return Ego(ptr[] ; ctxt =   geom.ctxt, refs = geom, delObj = true)
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
function getTessGeom(tess::Ego)
    sizes = zeros(Cint,2)
    xyz   = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_getTessGeom, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cdouble}}),
                tess._ego, sizes, xyz))
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
function makeTessBody(body::Ego, parms::vDbl)
    ptr   = Ref{ego}()
    raiseStatus(ccall((:EG_makeTessBody, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
            body._ego, parms,ptr))
    return Ego(ptr[] ; ctxt =   body.ctxt, refs = body, delObj = true)
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
function quadTess(tess::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_quadTess, C_egadslib), Cint, (ego, Ptr{ego}), tess._ego, ptr))
    return Ego(ptr[] ; ctxt =   tess.ctxt, refs = tess, delObj = true)
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
remakeTess!(tess::Ego, facedg::Vector{Ego}, parms::vDbl) = raiseStatus(ccall((:EG_remakeTess, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{Cdouble}),
                                                            tess._ego, Cint(length(facedg)), getfield.(facedg,:_ego), parms))



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
function locateTessBody(tess::Ego, ifaces::vInt, uv::Union{vDbl, mDbl}; mapped = false::Bool)

    npt     = length(ifaces)
    weights = zeros(3 * npt)
    itri    = mapped ? zeros(Cint, npt) : C_NULL
    idx     = collect(Cint, vcat(ifaces...))
    raiseStatus(ccall((:EG_locateTessBody, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}),
                tess._ego, npt, idx, vcat(uv...), itri, weights))
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
function getTessEdge(tess::Ego, eIndex::int)
    len = Ref{Cint}(0)
    xyz = Ref{Ptr{Cdouble}}()
    t   = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_getTessEdge, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}),
                tess._ego, Cint(eIndex), len, xyz, t))
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
function getTessFace(tess::Ego, fIndex::int)
    len    = Ref{Cint}(0)
    ntri   = Ref{Cint}(0)

    xyz    = Ref{Ptr{Cdouble}}()
    uv     = Ref{Ptr{Cdouble}}()
    ptype  = Ref{Ptr{Cint}}()
    pindex = Ref{Ptr{Cint}}()
    tric   = Ref{Ptr{Cint}}()
    tris   = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessFace, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}),
                tess._ego, Cint(fIndex), len, xyz, uv, ptype, pindex, ntri, tris, tric))
    xyz  = ptr_to_array(xyz[] , len[],3)
    uv   = ptr_to_array(uv[]  , len[],2)
    tris = int.(ptr_to_array(tris[],ntri[],3))
    tric = int.(ptr_to_array(tric[],ntri[],3))
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
function getTessLoops(tess::Ego, fIndex::int)
    ptr1 = Ref{Cint}(0)
    ptr2 = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessLoops, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}),
                tess._ego, Cint(fIndex),ptr1,ptr2))
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
makeQuads!(tess::Ego, qparms::vDbl, fIndex::int) = raiseStatus(ccall((:EG_makeQuads, C_egadslib), Cint, (ego, Ptr{Cdouble}, Cint),
                                                        tess._ego,qparms, Cint(fIndex)))

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
function getTessQuads(tess::Ego)
    nface = Ref{Cint}(0)
    fList = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_getTessQuads, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cint}}),
                tess._ego, nface,fList))
    arr   = ptr_to_array(fList[], nface[])
    _free(fList)
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
function getQuads(tess::Ego, fIndex::int)
    len    = Ref{Cint}(0)
    xyz    = Ref{Ptr{Cdouble}}()
    uv     = Ref{Ptr{Cdouble}}()
    ptype  = Ref{Ptr{Cint}}()
    pindex = Ref{Ptr{Cint}}()
    npatch = Ref{Cint}(0)
    raiseStatus(ccall((:EG_getQuads, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}),
                 tess._ego, Cint(fIndex), len, xyz, uv, ptype, pindex, npatch))
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
function getPatch(tess::Ego, fIndex::int, pIndex::int)
    nu     = Ref{Cint}(0)
    nv     = Ref{Cint}(0)
    ipts   = Ref{Ptr{Cint}}(0)
    bounds = Ref{Ptr{Cint}}(0)
    raiseStatus(ccall((:EG_getPatch, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}),
                tess._ego, Cint(fIndex), Cint(pIndex), nu, nv, ipts, bounds))
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
moveEdgeVert!(tess::Ego, eIndex::int, vIndex::int, t::dbl) = raiseStatus(ccall((:EG_moveEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cdouble),
                                                                                                 tess._ego, Cint(eIndex), Cint(vIndex),t))

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
deleteEdgeVert!(tess::Ego, eIndex::int, vIndex::int, dir::int) = raiseStatus(ccall((:EG_deleteEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cint),
                                                                                                                 tess._ego, Cint(eIndex), Cint(vIndex), Cint(dir)))

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
insertEdgeVerts!(tess::Ego, eIndex::int, vIndex::int,
                 ts::Union{dbl, vDbl}) = raiseStatus(ccall((:EG_insertEdgeVerts, C_egadslib), Cint, (ego, Cint, Cint, Cint, Ptr{Cdouble}),
                                                     tess._ego, Cint(eIndex), Cint(vIndex), Cint(length(vcat([ts]...))), vcat([ts]...)))


"""
Opens an existing Tessellation Object for replacing EDGE/FACE discretizations.

Parameters
----------
tess:
    the Tessellation Object to open
"""
openTessBody(tess::Ego) = raiseStatus(ccall((:EG_openTessBody, C_egadslib), Cint, (ego,), tess._ego))


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
function initTessBody(body::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_initTessBody, C_egadslib), Cint, (ego, Ptr{ego}), body._ego, ptr))
    return Ego(ptr[] ; ctxt =   body.ctxt, refs = body, delObj = true)
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
function statusTessBody(tess::Ego)
    ptr = Ref{ego}()
    np  = Ref{Cint}(0)
    st  = Ref{Cint}(0)
    stat = ccall((:EG_statusTessBody, C_egadslib), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}),
                tess._ego, ptr, st, np)
    (stat != EGADS_SUCCESS && stat != EGADS_OUTSIDE) && raiseStatus(stat)
    return ( body = Ego(ptr[] ; ctxt =   tess.ctxt, delObj = false), stat = st[], npts = np[])
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
setTessEdge!(tess::Ego, eIndex::int, xyz::mDbl, t::vDbl) = raiseStatus(ccall((:EG_setTessEdge, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}),
                                      tess._ego, Cint(eIndex), Cint(length(t)), vcat(xyz...), t))


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
function setTessFace!(tess::Ego, fIndex::int, xyz::mDbl, uv::mDbl, tris::mInt)
    len1 = Cint(length(vec(uv))   / 2)
    len2 = Cint(length(vec(tris)) / 3)
    return raiseStatus(ccall((:EG_setTessFace, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}),
                        tess._ego, Cint(fIndex), len1, vcat(xyz...),vcat(uv...), len2, Cint.(vcat(tris...))))
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
finishTess!(tess::Ego, params::vDbl) = raiseStatus(ccall((:EG_finishTess, C_egadslib), Cint,
                                         (ego, Ptr{Cdouble}), tess._ego, params))


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
function localToGlobal(tess::Ego, ind::int, locl::int)
    ptr  = Ref{Cint}(0)
    raiseStatus(ccall((:EG_localToGlobal, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}),
                tess._ego, Cint(ind), Cint(locl), ptr));
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
function getGlobal(tess::Ego, gbl::int)
    ptr1 = Ref{Cint}(0)
    ptr2 = Ref{Cint}(0)
    xyz  = zeros(3)
    raiseStatus(ccall((:EG_getGlobal, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}),
                tess._ego, Cint(gbl), ptr1, ptr2, xyz))
    return (ptype = ptr1[], pindex = ptr2[], xyz = xyz)
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
function effectiveMap(eobject::Ego, eparam::Union{vDbl,dbl})
    obj     = Ref{ego}()
    iparam  = vcat([eparam]...)
    oparam  = zeros(2)
    infoOBJ = unsafe_load(eobject._ego)
    if infoOBJ.oclass == EEDGE
        length(iparam) != 1 &&
            throw(ErrorException(" effectiveMap for EEDGE shoudl have a single eparam $iparam"))
    elseif infoOBJ.oclass == EFACE
        length(iparam) != 2 &&
            throw(ErrorException(" effectiveMap for EFACE shoudl have 2 eparams $iparam"))
    else
        throw(ErrorException(" effectiveMap oclass shoudl be either EFACE or EEDGE  and not $(infoOBJ.oclass)"))
    end
    raiseStatus(ccall((:EG_effectiveMap, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}, Ptr{Cdouble}),
                eobject._ego, iparam, obj, oparam))
    return (obj = Ego(obj[], ctxt = eobject.ctxt, delObj = false), param = infoOBJ.oclass == EEDGE ? oparam[1] : oparam )
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
function effectiveEdgeList(eedge::Ego)
    nedge  = Ref{Cint}(0)
    edges  = Ref{Ptr{ego}}()
    senses = Ref{Ptr{Cint}}()
    tstart = Ref{Ptr{Cdouble}}()
    raiseStatus(ccall((:EG_effectiveEdgeList, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}),
                eedge._ego, nedge, edges, senses, tstart))
    arr_e = ptr_to_array(edges[], nedge[])
    arr_s = ptr_to_array(senses[], nedge[])
    arr_t = ptr_to_array(tstart[], nedge[])
    _free(edges); _free(senses); _free(tstart)
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
    try @assert icode == EGADS_SUCCESS
    catch
        idx = findall(x-> x == icode, getfield.(egadsERROR,:var))
        msg = "EGADS error: "*egadsERROR[idx[1]].name
        throw(ErrorException(msg))
    end
end
