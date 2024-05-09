#=
*
*              jlEGADS --- Julia version of EGADS API
*
*
*      Copyright 2011-2024, Massachusetts Institute of Technology
*      Licensed under The GNU Lesser General Public License, version 2.1
*      See http://www.opensource.org/licenses/lgpl-2.1.php
*      Written by: Julia Docampo Sanchez
*      Email: julia.docampo@bsc.es
*
*
=#

__precompile__()

module egads

function __init__()
    if "egadslite" in string.(Base.loaded_modules_array())
        throw(ErrorException("egadslite has already been loaded! egads and egadslite cannot be loaded at the same time." ))
    end
end

using Libdl

# Dependencies
const EGADS_COMMON = joinpath(@__DIR__, "..", "..", "egadscommon", "egadscommon.jl")
if !isfile(EGADS_COMMON)
    throw(ErrorException("egadscommon.jl file not found !! got $EGADS_COMMON"))
end

include(EGADS_COMMON)

#Set library path

const egadsLIB = Base.Sys.isapple()   ? "libegads.dylib" :
                 Base.Sys.islinux()   ? "libegads.so"    :
                 Base.Sys.iswindows() ? "egads.dll"      :
                 error("Unsupported operating system!")

C_egadslib = Libdl.find_library(egadsLIB, [ESP_ROOT_LIB])
if isempty(C_egadslib)
    throw(ErrorException(" EGADS library $egadsLIB not found !!"))
end



"""
Resets the Context's owning thread ID to the thread calling this function.

Parameters
----------
context:
    the Context Object to update
"""
updateThread!(context::Context) = raiseStatus(ccall((:EG_updateThread, C_egadslib), Cint, (ego,), context.ego))


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
function makeTransform(context::Context, mat::mDbl)
    mat2 = vcat(mat...)
    ptr  = Ref{ego}()
    raiseStatus(ccall((:EG_makeTransform, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
                context.ego, mat2, ptr))
        return Ego(ptr[] ; ctxt = context, delObj = true)
end


"""
Returns the [4][3] transformation information.

Parameters
----------
xform:
    the Transformation Object.

"""
function getTransform(xform::Ego)
    mat = zeros(12)
    raiseStatus(ccall((:EG_getTransformation, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                xform._ego , mat))
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
function flipObject(object::Ego)
 ptr = Ref{ego}()
 raiseStatus(ccall((:EG_flipObject, C_egadslib), Cint, (ego, Ptr{ego}), object._ego , ptr))
 return Ego(ptr[] ; ctxt =  object.ctxt, delObj = true)
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
function copyObject(object::Ego; other=nothing::Union{Ego,Nothing})
    ptr   = Ref{ego}()
    oform = emptyVar(other) ? C_NULL : typeof(other) == Ego ? other._ego  : other[1:3]
    raiseStatus(ccall((:EG_copyObject, C_egadslib), Cint, (ego, Ptr{Cvoid}, Ptr{ego}),
                 object._ego , oform, ptr))
    return Ego(ptr[], ctxt = object.ctxt, delObj = true)
end


"""
Copy an Object to the specified Context
        
This is useful in multithreaded settings when you wish to copy an object that exists in a different
Context/thread. Use copyObject when copying from the object's context/thread to the context specified by other.

Parameters
----------
context:
    the context to copy the Object into

object:
    the Object to copy

Returns
-------
newObject: the resultant new ego
"""
function contextCopy(context::Context; object::Ego)
    ptr   = Ref{ego}()
    raiseStatus(ccall((:EG_contextCopy, C_egadslib), Cint, (ego, ego, Ptr{ego}),
            context._ego, object._ego, ptr))
    return Ego(ptr[], ctxt = context, delObj = true)
end


"""
Create a stream of data serializing the objects in the Model.

Returns
-------

stream:  the byte-stream
"""
function exportModel(model::Ego)
    nbyte  = Ref{Cint}(0)
    stream = Ref{Ptr{UInt8}}()
    raiseStatus(ccall((:EG_exportModel, C_egadslib), Cint,(ego, Ptr{Cint}, Ptr{Ptr{UInt8}}),
                model._ego, nbyte, stream))
    return str = [unsafe_load(stream[], j) for j = 1:nbyte[]]
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
function attributeAdd!(object::Ego, name::str, data)
    ints  = C_NULL ; reals = C_NULL ; chars = C_NULL
    atype = nothing
    if typeof(data) == CSys
        atype = ATTRCSYS
        reals = data.val
        len   = length(data.val)
    else
        dataF = vcat([data]...)
        len = Cint(length(dataF))
        if typeof(dataF) == vDbl
            atype = ATTRREAL
            reals = deepcopy(dataF)
        elseif typeof(dataF) <: vInt
            atype = ATTRINT
            ints  = Cint.(dataF)
        elseif typeof(dataF) <: vStr
            atype = ATTRSTRING
            chars = length(dataF) == 1 ? data : dataF
        else
            throw(@error("attributeAdd Data = $data --> $(typeof(data))  is not all int, double, or str"))
        end
    end
    raiseStatus(ccall((:EG_attributeAdd, C_egadslib), Cint, (ego, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cstring),
                 object._ego , name, atype, len, ints, reals, chars) )
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
function attributeDel!(object::Ego; name=nothing::Union{Nothing,String})
    cname = emptyVar(name) ? C_NULL : string(name)
    raiseStatus(ccall((:EG_attributeDel, C_egadslib), Cint, (ego, Cstring), object._ego , cname))
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
attributeDup!(src::Ego, dst::Ego) = raiseStatus(ccall((:EG_attributeDup, C_egadslib), Cint, (ego, ego), src._ego , dst._ego ))


"""
Sets the attribution mode for the Context.

Parameters
----------
context:
    the Context Object

attrFlag:
    the mode flag: 0 - the default scheme, 1 - full attribution node
"""
setFullAttrs!(context::Context, attrFlag::int) = raiseStatus(ccall((:EG_setFullAttrs, C_egadslib), Cint, (ego, Cint), context.ego, Cint(attrFlag)))


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
function attributeAddSeq!(object::Ego, name::str, data)
    ints  = C_NULL ; reals = C_NULL ; chars = C_NULL
    atype = nothing
    if typeof(data) <: CSys
        atype = ATTRCSYS
        reals = data.val
        len = length(reals)
    else
        dataF = vcat([data]...)
        len = length(dataF)
        if typeof(dataF) == vDbl
            atype = ATTRREAL
            reals = deepcopy(dataF)
        elseif typeof(dataF) <: vInt
            atype = ATTRINT
            ints  = Cint.(dataF)
        elseif typeof(dataF) == vStr
            atype = ATTRSTRING
            chars = length(dataF) == 1 ? data : dataF
        else
            throw(ErrorMessage(" attributeAddSeq! Data = $data  is not all int, double, or str"))
        end
    end
    s = ccall((:EG_attributeAddSeq, C_egadslib), Cint, (ego, Cstring, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Cstring),
               object._ego , name, atype, len, ints, reals, chars)
    s < 0 && raiseStatus(s)
    return s
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
function makeGeometry(context::Context, oclass::int, mtype::int, reals::Union{dbl,mDbl, vDbl};
                      rGeom=nothing::Union{Nothing,Ego},ints=nothing::Union{Nothing,vInt})

    reals = vcat([reals]...)
    if (oclass != CURVE && oclass != PCURVE && oclass != SURFACE)
        throw(ErrorMessage("makeGeometry only accepts CURVE, PCURVE or SURFACE"))
    end
    geo  = rGeom === nothing ? C_NULL : rGeom._ego
    vint = ints  === nothing ? C_NULL : Cint.(ints)

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
    else
        if oclass == PCURVE
            if   mtype == LINE          nreal = 4
            elseif mtype == CIRCLE      nreal = 7
            elseif mtype == ELLIPSE     nreal = 8
            elseif mtype == PARABOLA    nreal = 7
            elseif mtype == HYPERBOLA   nreal = 8
            elseif mtype == TRIMMED     nreal = 2
            elseif mtype == OFFSET      nreal = 1
            else throw(ErrorMessage("makeGeometry Unknown mtype $mtype"))
            end
        elseif oclass == CURVE
            if   mtype == LINE          nreal = 6
            elseif mtype == CIRCLE      nreal = 10
            elseif mtype == ELLIPSE     nreal = 11
            elseif mtype == PARABOLA    nreal = 10
            elseif mtype == HYPERBOLA   nreal = 11
            elseif mtype == TRIMMED     nreal = 2
            elseif mtype == OFFSET      nreal = 4
            else throw(ErrorMessage("makeGeometry Unknown mtype $mtype"))
            end
        elseif oclass == SURFACE
            if   mtype == PLANE         nreal = 9
            elseif mtype == SPHERICAL   nreal = 10
            elseif mtype == CONICAL     nreal = 14
            elseif mtype == CYLINDRICAL nreal = 13
            elseif mtype == TOROIDAL    nreal = 14
            elseif mtype == OFFSET      nreal = 1
            elseif mtype == TRIMMED     nreal = 4
            elseif mtype == EXTRUSION   nreal = 3
            elseif mtype == REVOLUTION  nreal = 6
            else throw(ErrorMessage("makeGeometry Unknown mtype $mtype"))
            end
        end

        if nreal != lr
            throw(ErrorMessage("makeGeometry length(reals) == $lr must be $nreal for oclass $oclass and mtype $mtype"))
        end
    end

    raiseStatus(ccall((:EG_makeGeometry, C_egadslib), Cint, (ego, Cint, Cint, ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}),
                context.ego, Cint(oclass), Cint(mtype), geo, vint, reals, ptr))
    return Ego(ptr[] ; ctxt =  context, refs = rGeom, delObj = true)
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
function skinning(curves::Vector{Ego}, degree::int)
    ptr     = Ref{ego}()
    raiseStatus(ccall((:EG_skinning, C_egadslib), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}),
                length(curves), getfield.(curves,:_ego ), Cint(degree), ptr))
    return Ego(ptr[]; ctxt =  curves[1].ctxt, delObj = true)
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
function approximate(context::Context, sizes::vInt, xyzs::mDbl ;
                     deg=0::Int, tol=1.e-08::dbl)
    ptr      = Ref{ego}()
    raiseStatus(ccall((:EG_approximate, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}),
                context.ego, Cint(deg), tol, Cint.(sizes), vcat(xyzs...), ptr))
    return Ego(ptr[] ; ctxt =  context, delObj = true)
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
function fitTriangles(context::Context, xyzs::Union{vDbl,mDbl}, tris::Union{vInt,mInt} ;
                      tric=nothing::Union{Nothing,vInt}, tol=1.e-08::dbl)

    vint = tric === nothing ? C_NULL : Cint.(tric)

    flat_xyz  = vcat(xyzs...)
    flat_tris = Cint.(vcat(tris...))
    ntris     = Cint(length(flat_tris)/ 3)
    nxyz      = Cint(length(flat_xyz) / 3)

    if vint != C_NULL
        ntric = length(vint)
        Int(ntric / 3 ) != ntris &&
            (throw(ErrorMessage("fitTriangles tric needs to be integer vector of length $ntris != $ntric")))
    end
    ptr       = Ref{ego}()
    raiseStatus(ccall((:EG_fitTriangles, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{ego}),
                context.ego, nxyz, flat_xyz, ntris, flat_tris, vint, tol, ptr))
    return Ego(ptr[] ; ctxt =  context, delObj = true)
end


"""
Insert additional Knots into a BSpline.

Parameters
----------
bspline:
    a B-Spline Curve or Surface

minDegU:
    the minimum U-degree in the resulting B-Spline (Degree for Curve) 

Us:
    the U-Knot values to insert (may be repeated) (Knots for Curve)

minDegV:
    the minimum V-degree in the resulting B-Spline (Must be 0 for Curve) 

Vs:
    the V-Knot values to insert (may be repeated) (Must be None for Curve)

Returns
-------
The BSpline resultant ego
"""
function addKnots(bspline::Ego, minDegU::int, Us::Union{vDbl,Nothing}, minDegV::int, Vs::Union{vDbl,Nothing})
    nU = Us === nothing ? Cint(0) : Cint(length(Us))
    nV = Vs === nothing ? Cint(0) : Cint(length(Vs))
    Us = Us === nothing ? C_NULL : vcat(Us...)
    Vs = Vs === nothing ? C_NULL : vcat(Vs...)
    result = Ref{ego}()
    raiseStatus(ccall((:EG_addKnots, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Cint, Cint, Ptr{Cdouble}, Ptr{ego}),
                bspline._ego, 
                Cint(minDegU), nU, Us, 
                Cint(minDegV), nV, Vs, 
                result))
    return Ego(result[] ; ctxt =  bspline.ctxt, delObj = true)
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
function isoCline(surface::Ego, iUV::int, value::dbl)
    c = Ref{ego}()
    raiseStatus(ccall((:EG_isoCline, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{ego}),
                surface._ego , Cint(iUV), value, c))
    return Ego(c[] ; ctxt =  surface.ctxt, delObj = true)
end


"""
Check if a PCurve is an isocline

Parameters
----------
pcurve:
    the PCurve Object

Returns
-------
isIso:
    True of False if the PCurve is an isocline
    
iUV:
    the type of isocline: UISO (0) constant U or VISO (1) constant V

value:
    the value for the isocline

fwd:
    direction of the PCurve (SFORWARD (1) - Forward, SREVERSE (-1) - Reversed)
"""
function isIsoPCurve(pcurve::Ego)
    iUV   = Ref{Cint}(-1)
    value = Ref{Cdouble}()
    fwd   = Ref{Cint}(0)
    status = ccall((:EG_isIsoPCurve, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}),
             pcurve._ego , iUV, value, fwd)
    if (status < 0) 
        raiseStatus(status)
    end
    return status == EGADS_SUCCESS, iUV[], value[], fwd[]
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
function otherCurve(surface::Ego, curve::Ego ; tol::dbl=0.0)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_otherCurve, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{ego}),
                surface._ego , curve._ego , tol,ptr))
    return Ego(ptr[] ; ctxt =  surface.ctxt, delObj = true)
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
isSame(object1::Ego, object2::Ego) = return EGADS_SUCCESS == ccall((:EG_isSame, C_egadslib), Cint, (ego, ego),object1._ego ,object2._ego )


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
function convertToBSpline(object::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_convertToBSpline, C_egadslib), Cint, (ego, Ptr{ego}), object._ego , ptr))
    return Ego(ptr[] ; ctxt =  object.ctxt, delObj = (object._ego  != ptr[] ))
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
function convertToBSplineRange(geom::Ego, range::Union{vDbl,mDbl})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_convertToBSplineRange, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}),
                geom._ego , vcat(range...), ptr))
    return Ego(ptr[] ; ctxt =  geom.ctxt, delObj = true)
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
function makeTopology!(context::Context, oclass::int;
                       mtype=0::int, geom=nothing::Union{Nothing, Ego},
                       reals=nothing::Union{Nothing,vDbl,mDbl},
                       senses=nothing::Union{Nothing,vInt},
                       children=nothing::Union{Nothing, Ego, Vector{Ego}})
    oclass   = Cint(oclass)
    mtype    = Cint(mtype)
    geo      = geom   === nothing ? C_NULL : geom._ego
    rea      = reals  === nothing ? C_NULL : vcat(reals...)
    sen      = senses === nothing ? C_NULL : Cint.(vcat(senses...))
    childObj = C_NULL ; nchild = 0
    childC   = vcat([children]...)

    if childC != [nothing]
        nchild = length(childC)
        childObj = getfield.(childC,:_ego )
    end

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
    (oclass == LOOP && geo != C_NULL) && (nchild = Cint(nchild/2))
    nsenses = sen == C_NULL ? 0 : length(sen)
    (nchild > 0 && nsenses > 0 && nchild != nsenses) &&
        throw(ErrorMessage(" Number of children is inconsistent with senses $nchild != $nsenses"))


     raiseStatus(ccall((:EG_makeTopology, C_egadslib), Cint,
     (ego, ego, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{ego}, Ptr{Cint}, Ptr{ego}),
      context.ego, geo, oclass, mtype, rea, nchild, childObj, sen, topo_ptr))

    (geo != C_NULL) && push!(refs,geom)

    if oclass == BODY
        refs = nothing
    elseif oclass == MODEL
        refs = nothing
        setfield!.(childC, :_delObj, false)
    end

    topo = Ego(topo_ptr[] ; ctxt =  context, refs = refs, delObj = true)

    # The model takes ownership of bodies, so the bodies must reference the model
    if oclass == MODEL
        for ref in vcat([children]...)
            push!(ref._refs, topo)
        end
    end

    return topo
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
function makeFace(object::Ego, mtype::int ; rdata::Union{mDbl,vDbl,Nothing} = nothing)
    limits = typeof(rdata) == mDbl ? vcat(rdata'...) :
             typeof(rdata) == vDbl ? deepcopy(rdata) : C_NULL

    ptr    = Ref{ego}()
    icode  = ccall((:EG_makeFace, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}),
                    object._ego , Cint(mtype), limits, ptr)

    icode < 0 && raiseStatus(icode)

    return Ego(ptr[]; ctxt =  object.ctxt, refs = object, delObj = true)
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
    !! The edges used in loop are converted to Ego()

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
function makeLoop!(edges::Union{Ego, Vector{Ego}};
                   geom = nothing::Union{Nothing,Ego},
                   toler = 0.0::dbl)

    geo  = geom == nothing ? C_NULL : geom._ego
    nedge    = length(edges)
    loop_ptr = Ref{ego}()

    # Create a copy of the edge pointers
    piv  = findall(x -> !emptyVar(x), edges )
    # Create a copy of the edge pointers
    ePtr = [edges[piv[i]]._ego  for i = 1:length(piv) ]

    ctxt = edges[piv[1]].ctxt

    # This function makes ePtr .== C_NULL !!
    s    = ccall((:EG_makeLoop, C_egadslib), Cint,(Cint, Ptr{ego}, ego, Cdouble, Ptr{ego}),
                nedge, ePtr, geo, toler, loop_ptr)
    pivT        = findall(x -> emptyVar(x), ePtr)
    pivF        = findall(x -> !emptyVar(x), ePtr)
    loop = Ego(loop_ptr[]; ctxt = ctxt, refs = edges[pivT], delObj = true)
    s < EGADS_SUCCESS && raiseStatus(s)
    return loop, edges[pivF]
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
function makeSolidBody(context::Context, stype::int, data::Union{vDbl,vInt})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeSolidBody, C_egadslib), Cint, (ego, Cint,Ptr{Cdouble}, Ptr{ego}),
                context.ego, Cint(stype), Cdouble.(data), ptr))
    return Ego(ptr[]; ctxt =  context, delObj = true)
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
    area = Ref{Cdouble}(0.0)
    lims = emptyVar(limits) ? C_NULL : vcat(limits'...)
    info = getInfo(object)
    if info.oclass == SURFACE && (lims == C_NULL || length(lims) != 4)
        throw(ErrorMessage("getArea for surfaces needs 4 limits (u = limits[1,:] v = limits[2,:] ) "))
    end
    raiseStatus(ccall((:EG_getArea, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}),
                object._ego , lims, area))
    return area[]
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
function getMassProperties(object::Ego)
    props = zeros(14)
    raiseStatus(ccall((:EG_getMassProperties, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                object._ego , props))
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
isEquivalent(topo1::Ego, topo2::Ego) = return EGADS_SUCCESS == ccall((:EG_isEquivalent, C_egadslib), Cint, (ego, ego), topo1._ego , topo2._ego )


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
function sewFaces(objects::Vector{Ego}; toler=0.0::dbl, flag = true::Bool)

    ptr   = Ref{ego}()
    raiseStatus(ccall((:EG_sewFaces, C_egadslib), Cint, (Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}),
                 length(objects),getfield.(objects,:_ego ), toler, flag, ptr))
    return Ego(ptr[] ; ctxt =  objects[1].ctxt, refs = objects, delObj = true)
end


"""
Creates a non-manifold Wire Body

Parameters
----------
objects:
List of Edge Objects to make the Wire Body

toler:
Node tolerance to connect Edges (0.0 indicates the use of the Nodes directly)

Returns
-------
The resultant Wire Body Object
"""
function makeNmWireBody(objects::Vector{Ego}; toler=0.0::dbl)
    nobj  = length(objects)
    objs  = getfield.(objects,:_ego )
    wbody = Ref{ego}()
    raiseStatus(ccall((:EG_makeNmWireBody, C_egadslib), Cint, (Cint, Ptr{ego}, Cdouble, Ptr{ego}), nobj, objs, toler, wbody))
    return Ego(wbody[] ; ctxt =  objects[1].ctxt,refs = objects, delObj = true)
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

objects:
    list of FACE tuple pairs, where the first must be a FACE in the BODY
    and second is either the FACE to use as a replacement or a Nothing which
    indicates that the FACE is to be removed from the BODY

Returns
-------
result: the resultant BODY object, either a SHEETBODY or a SOLIDBODY
        (where the input was a SOLIDBODY and all FACEs are replaced
        in a way that the LOOPs match up)
"""
function replaceFaces(body::Ego, objects)

    obj_flat = vcat(objects...)
    nfaces   = Cint(length(obj_flat) / 2)
    (length(obj_flat) != 2 * nfaces) &&
        (throw(ErrorMessage("replaceFaces object is not a pair ")))
    pfaces = Array{ego}(undef, 2 * nfaces)
    for i = 1: 2 * nfaces
        pfaces[i] = emptyVar(obj_flat[i]) ? C_NULL : obj_flat[i]._ego
    end
    res    = Ref{ego}()
    info   = getInfo(body)
    raiseStatus(ccall((:EG_replaceFaces, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
                       body._ego ,nfaces, pfaces,res))
    return Ego(res[] ; ctxt =  body.ctxt, delObj = true)
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
function matchBodyEdges(body1::Ego, body2::Ego ; toler=0.0::dbl)
    nMatch  = Ref{Cint}(0)
    matches = Ref{Ptr{Cint}}()
    v1 = unsafe_load(body1._ego )
    v2 = unsafe_load(body2._ego )

    (v1.oclass != BODY || v2.oclass != BODY) &&
        (throw(ErrorMessage("matchBodyFaces Both objects must have object class BODY")))
    raiseStatus(ccall((:EG_matchBodyEdges, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}),
                body1._ego , body2._ego , toler, nMatch, matches))
    arr  = ptr_to_array(matches[],nMatch[],2)
    _free(matches)
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
function matchBodyFaces(body1::Ego, body2::Ego; toler=0.0::dbl)
    nMatch  = Ref{Cint}(0)
    matches = Ref{Ptr{Cint}}()
    v1 = unsafe_load(body1._ego )
    v2 = unsafe_load(body2._ego )

    (v1.oclass != BODY || v2.oclass != BODY) &&
        (throw(ErrorMessage("matchBodyFaces Both objects must have object class BODY")))

    raiseStatus(ccall((:EG_matchBodyFaces, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}),
                body1._ego , body2._ego , toler, nMatch, matches))
    arr = ptr_to_array(matches[],nMatch[],2)
    _free(matches)
    return arr
end


"""
Checks for topological equivalence between the the BODY src and the BODY dst.
If necessary, produces a mapping (indices in src which map to dst) and
places these as attributes on the resultant BODY mapped (named .nMap, .eMap and .fMap).
Unlike mapBody, mapBody2 also works on FACEBODYs and WIREBODYs.

Parameters
----------
src:
    the source Body Object (not a WIREBODY)

fAttr:
    the FACE attribute string used to map FACEs

eAttr:
    the EDGE attribute string used to map EDGEs

dst:
    destination Body Object

Returns
-------
Note: It is the responsibility of the caller to have
      uniquely attributed all FACEs in both src and dst
      to aid in the mapping.
"""
function mapBody2(src::Ego, fAttr::str, eAttr::str, dst::Ego)
    raiseStatus(ccall((:EG_mapBody2, C_egadslib), Cint, (ego, Cstring, Cstring, ego),
                src._ego, fAttr, eAttr, dst._ego))
end


"""
Checks for topological equivalence between the the BODY src and the BODY dst.
If necessary, produces a mapping (indices in src which map to dst) and
places these as attributes on dst (named .nMap, .eMap and .fMap).

Parameters
----------
src:
    the source Body Object

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
function mapBody(src::Ego, dst::Ego, fAttr::str)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_mapBody, C_egadslib), Cint, (ego, ego, Cstring, Ptr{ego}),
                src._ego , dst._ego , fAttr, ptr))
    if emptyVar(ptr[])
        return nothing
    else
        return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
    end
end


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
function mapTessBody(tess::Ego, body::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_mapTessBody, C_egadslib), Cint, (ego, ego, Ptr{ego}),
                tess._ego , body._ego , ptr))
    return Ego(ptr[] ; ctxt =  tess.ctxt, refs = body, delObj = true)
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
function tessMassProps(tess::Ego)
    props = zeros(14)
    raiseStatus(ccall((:EG_tessMassProps, C_egadslib), Cint, (ego, Ptr{Cdouble}),
                tess._ego , props))
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
function generalBoolean(src::Ego, tool::Ego, oper::int ; tol=0.0::dbl)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_generalBoolean, C_egadslib), Cint, (ego, ego,Cint, Cdouble, Ptr{ego}),
                src._ego , tool._ego , Cint(oper), tol, ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function fuseSheets(src::Ego, tool::Ego)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_fuseSheets, C_egadslib), Cint, (ego, ego, Ptr{ego}),
                src._ego , tool._ego , ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function intersection(src::Ego, tool::Ego)
    nEdge = Ref{Cint}(0)
    pairs = Ref{Ptr{ego}}()
    model = Ref{ego}()
    raiseStatus(ccall((:EG_intersection, C_egadslib), Cint, (ego, ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{ego}),
                src._ego , tool._ego , nEdge, pairs, model))

    if nEdge[] == 0
        _free(pairs)
        return (pairs = Vector{Ego}(), model = nothing)
    end
    arr = ptr_to_array(pairs[], nEdge[] * 2)
    _free(pairs)
    return (pairs = reshape([Ego(arr[j]; ctxt =  src.ctxt, delObj = false) for j =1:length(arr)],2,Int(nEdge[])),
            model = Ego(model[]; ctxt = src.ctxt, delObj = true) )
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
function imprintBody(src::Ego, pairs::Matrix{Ego})
    ps     = Cint(size(pairs,2))
    #fpairs = vcat(ego.(getfield.(pairs,:_ego ))...)
    fpairs = vcat(getfield.(pairs,:_ego )...)
    ptr    = Ref{ego}()
    raiseStatus(ccall((:EG_imprintBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
               src._ego , ps, fpairs, ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function filletBody(src::Ego, edges::Vector{Ego}, radius::dbl)
    nedge = length(edges)
    res   = Ref{ego}()
    map   = Ref{Ptr{Cint}}()
    raiseStatus(ccall((:EG_filletBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}),
                src._ego , nedge, getfield.(edges,:_ego ),radius, res,map))
    Cres  = Ego(res[] ; ctxt =  src.ctxt, delObj = true)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(map[], length(bodyT), 2)
    _free(map)
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
function chamferBody(src::Ego, edges::Vector{Ego}, faces::Vector{Ego}, dis1::dbl, dis2::dbl)
    n     = length(edges)
    n    != length(faces) && (throw(ErrorMessage("edges and faces must be same length")))
    fmap  = Ref{Ptr{Cint}}()
    res   = Ref{ego}()
    raiseStatus(ccall((:EG_chamferBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}, Cdouble, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}),
                src._ego , n, getproperty.(edges,:_ego ),getproperty.(faces,:_ego ), dis1, dis2, res, fmap))
    Cres  = Ego(res[]; ctxt =  src.ctxt, delObj = true)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(fmap[], length(bodyT),2)
    _free(fmap)
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
function hollowBody(src::Ego, faces::Union{Ego,Vector{Ego}, Nothing}, offset::dbl, join::int)
    Cface = emptyVar(faces) ? C_NULL : getfield.(vcat([faces]...),:_ego )
    nface = Cface == C_NULL ? 0 : length(Cface)
    fmap  = Ref{Ptr{Cint}}()
    res   = Ref{ego}()
    raiseStatus(ccall((:EG_hollowBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}, Ptr{Ptr{Cint}}),
                src._ego , nface, Cface,offset, Cint(join), res, fmap))
    Cres  = Ego(res[]; ctxt =  src.ctxt, delObj = true)
    bodyT = getBodyTopos(Cres, FACE)
    arr   = ptr_to_array(fmap[], length(bodyT), 2)
    _free(fmap)
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
function rotate(src::Ego, angle::Union{int,dbl}, axis::Union{vDbl,mDbl})
    ptr = Ref{ego}()
    ax  = vcat(axis...)
    length(ax) != 6 && throw(ErrorMessage(" rotate: flat axis length != 6 "))
    raiseStatus(ccall((:EG_rotate, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}),
                src._ego , Cdouble(angle),ax, ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function extrude(src::Ego, dist::dbl, dir::vDbl)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_extrude, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}),
                  src._ego ,dist, dir, ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function sweep(src::Ego, spine::Ego, mode::int)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_sweep, C_egadslib), Cint, (ego, ego, Cint, Ptr{ego}),
                src._ego , spine._ego , Cint(mode), ptr))
    return Ego(ptr[] ; ctxt =  src.ctxt, delObj = true)
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
function ruled(secs::Vector{Ego})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_ruled, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{ego}),
                Cint(length(secs)), vcat(getfield.(secs,:_ego )...), ptr))
    return Ego(ptr[] ; ctxt =  secs[1].ctxt, delObj = true)
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
function blend(secs::Vector{Ego} ;
               rc1=nothing::Union{Nothing,vDbl},rc2=nothing::Union{Nothing,vDbl})
    ptr = Ref{ego}()

    prc1 = rc1 == nothing ? C_NULL : rc1[1:min(length(rc1),8)]
    prc2 = rc2 == nothing ? C_NULL : rc2[1:min(length(rc2),8)]

    raiseStatus(ccall((:EG_blend, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ego}),
                Cint(length(secs)), getfield.(secs,:_ego ),prc1, prc2, ptr))
    return Ego(ptr[] ; ctxt =  secs[1].ctxt, delObj = true)
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
function initEBody(tess::Ego, angle::dbl)
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_initEBody, C_egadslib), Cint, (ego, Cdouble, Ptr{ego}),
                tess._ego , angle, ptr))
    return Ego(ptr[] ; ctxt =  tess.ctxt, refs = tess, delObj = true)
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
function makeAttrEFaces(EBody::Ego, attrName::str)
    neface = Ref{Cint}(0)
    efaces = Ref{Ptr{ego}}()
    raiseStatus(ccall((:EG_makeAttrEFaces, C_egadslib), Cint, (ego, Cstring, Ptr{Cint}, Ptr{Ptr{ego}}),
                EBody._ego , attrName, neface, efaces))
    arr    = ptr_to_array(efaces[], neface[])
    _free(efaces)
    return [Ego(arr[j] ; ctxt =  EBody.ctxt, delObj = false) for j = 1: length(arr)]
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
function makeEFace(EBody::Ego, faces::Union{Vector{Ego},Ego})
    ptr = Ref{ego}()
    raiseStatus(ccall((:EG_makeEFace, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}),
                 EBody._ego , length(faces), getfield.(vcat([faces]...),:_ego ), ptr))
    return Ego(ptr[] ; ctxt =  EBody.ctxt, delObj = false)
end


"""
Finish an EBody
"""
function finishEBody!(EBody::Ego)
    raiseStatus(ccall((:EG_finishEBody, C_egadslib), Cint, (ego,), EBody._ego ))
    EBody._refs === nothing && throw(ErrorMessage(" finishEBody, EBody doesn't have any ref "))

    # the EBody now references the body referenced by the tessellation
    refs = Vector{Ego}()
    for ref_tess in EBody._refs
        ref_tess._refCount -= 1
        for ref_body in ref_tess._refs
            ref_body._refCount += 1
            push!(refs, ref_body)
        end
    end
    EBody._refs = refs
end

end
