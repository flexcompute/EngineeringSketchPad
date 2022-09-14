struct egObject
    magicnumber::Cint
    oclass::Cshort
    mtype::Cshort
    attrs::Ptr{Cvoid}
    blind::Ptr{Cvoid}
    topObj::Ptr{egObject}
    tref::Ptr{egObject}
    prev::Ptr{egObject}
    next::Ptr{egObject}
end

const ego = Ptr{egObject}

struct KEY
    keys::NTuple{3, Cint}
end

struct DATA
    close::Cint
    xyz::NTuple{3, Cdouble}
end

struct ENTRY
    key::KEY
    data::DATA
end

struct element
    item::ENTRY
    next::Ptr{element}
end

const ELEMENT = element

struct triVert
    type::Cint
    edge::Cint
    index::Cint
    loop::Cint
    xyz::NTuple{3, Cdouble}
    uv::NTuple{2, Cdouble}
end

struct triTri
    indices::NTuple{3, Cint}
    neighbors::NTuple{3, Cint}
    mid::NTuple{3, Cdouble}
    area::Cdouble
    mark::Cshort
    close::Cshort
    hit::Cshort
    count::Cshort
end

struct triSeg
    indices::NTuple{2, Cint}
    neighbor::Cint
    edge::Cint
    index::Cint
end

struct triStruct
    face::Ptr{Cint}
    fIndex::Cint
    orUV::Cint
    planar::Cint
    phase::Cint
    VoverU::Cdouble
    maxlen::Cdouble
    chord::Cdouble
    dotnrm::Cdouble
    accum::Cdouble
    edist2::Cdouble
    eps2::Cdouble
    devia2::Cdouble
    minlen::Cdouble
    range::NTuple{4, Cdouble}
    qparm::NTuple{3, Cdouble}
    uvs::Ptr{Cdouble}
    orCnt::Cint
    maxPts::Cint
    mverts::Cint
    nverts::Cint
    verts::Ptr{triVert}
    mtris::Cint
    ntris::Cint
    tris::Ptr{triTri}
    msegs::Cint
    nsegs::Cint
    segs::Ptr{triSeg}
    nfrvrts::Cint
    mframe::Cint
    nframe::Cint
    frame::Ptr{Cint}
    mloop::Cint
    nloop::Cint
    loop::Ptr{Cint}
    lens::NTuple{5, Cint}
    numElem::Cint
    tfi::Cint
    hashTab::Ptr{Ptr{ELEMENT}}
end

struct Front
    sleft::Cint
    i0::Cint
    i1::Cint
    sright::Cint
    snew::Cshort
    mark::Cshort
end

struct fillArea
    mfront::Cint
    nfront::Cint
    mpts::Cint
    npts::Cint
    pts::Ptr{Cint}
    nsegs::Cint
    segs::Ptr{Cint}
    front::Ptr{Front}
end

struct connect
    node1::Cint
    node2::Cint
    tri::Ptr{Cint}
    thread::Cint
end

struct EMPtess
    mutex::Ptr{Cvoid}
    master::Clong
    _end::Cint
    index::Cint
    ignore::Cint
    silent::Cint
    mark::Ptr{Cint}
    tess::Ptr{Cint}
    btess::Ptr{Cint}
    body::Ptr{Cint}
    faces::Ptr{Ptr{Cint}}
    edges::Ptr{Ptr{Cint}}
    params::Ptr{Cdouble}
    tparam::Ptr{Cdouble}
    qparam::NTuple{3, Cdouble}
    ptr::Ptr{Cvoid}
end

struct var"##Ctag#360"
    data::NTuple{8, UInt8}
end

function Base.getproperty(x::Ptr{var"##Ctag#360"}, f::Symbol)
    f === :integer && return Ptr{Cint}(x + 0)
    f === :integers && return Ptr{Ptr{Cint}}(x + 0)
    f === :real && return Ptr{Cdouble}(x + 0)
    f === :reals && return Ptr{Ptr{Cdouble}}(x + 0)
    f === :string && return Ptr{Ptr{Cchar}}(x + 0)
    return getfield(x, f)
end

function Base.getproperty(x::var"##Ctag#360", f::Symbol)
    r = Ref{var"##Ctag#360"}(x)
    ptr = Base.unsafe_convert(Ptr{var"##Ctag#360"}, r)
    fptr = getproperty(ptr, f)
    GC.@preserve r unsafe_load(fptr)
end

function Base.setproperty!(x::Ptr{var"##Ctag#360"}, f::Symbol, v)
    unsafe_store!(getproperty(x, f), v)
end

struct egAttr
    name::Ptr{Cchar}
    type::Cint
    length::Cint
    vals::var"##Ctag#360"
end

struct egAttrSeq
    root::Ptr{Cchar}
    nSeq::Cint
    attrSeq::Ptr{Cint}
end

struct egAttrs
    nattrs::Cint
    attrs::Ptr{egAttr}
    nseqs::Cint
    seqs::Ptr{egAttrSeq}
end

struct egCntxt
    outLevel::Cint
    fixedKnots::Cint
    fullAttrs::Cint
    tess::NTuple{2, Cdouble}
    signature::Ptr{Ptr{Cchar}}
    usrPtr::Ptr{Cvoid}
    threadID::Clong
    mutex::Ptr{Cvoid}
    pool::Ptr{egObject}
    last::Ptr{egObject}
end

struct egFconn
    index::Cint
    nface::Cint
    faces::Ptr{Cint}
    tric::Ptr{Cint}
end

struct egBary
    tri::Cint
    w::NTuple{2, Cdouble}
end

struct egTess1D
    obj::Ptr{egObject}
    nodes::NTuple{2, Cint}
    faces::NTuple{2, egFconn}
    xyz::Ptr{Cdouble}
    t::Ptr{Cdouble}
    _global::Ptr{Cint}
    npts::Cint
end

struct egPatch
    ipts::Ptr{Cint}
    bounds::Ptr{Cint}
    nu::Cint
    nv::Cint
end

struct egTess2D
    mKnots::Ptr{egObject}
    xyz::Ptr{Cdouble}
    uv::Ptr{Cdouble}
    _global::Ptr{Cint}
    ptype::Ptr{Cint}
    pindex::Ptr{Cint}
    bary::Ptr{egBary}
    frame::Ptr{Cint}
    frlps::Ptr{Cint}
    tris::Ptr{Cint}
    tric::Ptr{Cint}
    patch::Ptr{egPatch}
    npts::Cint
    nframe::Cint
    nfrlps::Cint
    ntris::Cint
    npatch::Cint
    tfi::Cint
end

struct egTessel
    src::Ptr{egObject}
    xyzs::Ptr{Cdouble}
    tess1d::Ptr{egTess1D}
    tess2d::Ptr{egTess2D}
    globals::Ptr{Cint}
    params::NTuple{6, Cdouble}
    tparam::NTuple{2, Cdouble}
    nGlobal::Cint
    nEdge::Cint
    nFace::Cint
    nu::Cint
    nv::Cint
    done::Cint
end

struct egEMap
    nobjs::Cint
    objs::Ptr{Ptr{egObject}}
end

struct egEdVert
    edge::Ptr{egObject}
    curve::Cint
    npts::Cint
    ts::Ptr{Cdouble}
    dstart::NTuple{3, Cdouble}
    dend::NTuple{3, Cdouble}
end

struct egEEseg
    iedge::Cint
    sense::Cint
    nstart::Ptr{egObject}
    tstart::Cdouble
    tend::Cdouble
end

struct egEEdge
    sedges::Ptr{egEdVert}
    nsegs::Cint
    segs::Ptr{egEEseg}
    trange::NTuple{2, Cdouble}
    nodes::NTuple{2, Ptr{egObject}}
end

struct egEdgeUV
    edge::Ptr{egObject}
    sense::Cint
    npts::Cint
    iuv::Ptr{Cint}
end

struct egELoop
    eedges::egEMap
    senses::Ptr{Cint}
    area::Cdouble
    nedge::Cint
    edgeUVs::Ptr{egEdgeUV}
end

struct egEPatch
    face::Ptr{egObject}
    tol::Cdouble
    start::Cint
    nuvs::Cint
    ndeflect::Cint
    ntris::Cint
    uvtris::Ptr{Cint}
    uvtric::Ptr{Cint}
    uvs::Ptr{Cdouble}
    deflect::Ptr{Cdouble}
end

struct egEFace
    sedges::Ptr{egEdVert}
    npatch::Cint
    patches::Ptr{egEPatch}
    eloops::egEMap
    senses::Ptr{Cint}
    trmap::Ptr{Cint}
    uvmap::Ptr{Cvoid}
    range::NTuple{4, Cdouble}
    last::Cint
end

struct egEShell
    efaces::egEMap
end

struct egEBody
    ref::Ptr{egObject}
    eedges::egEMap
    eloops::egEMap
    efaces::egEMap
    eshells::egEMap
    senses::Ptr{Cint}
    angle::Cdouble
    done::Cint
    nedge::Cint
    edges::Ptr{egEdVert}
end

# Skipping MacroDefinition: __ProtoExt__ __HOST_AND_DEVICE__ extern

const EGADS_EXTRAPOL = Int32(-37)
const EGADS_EFFCTOBJ = Int32(-36)
const EGADS_UVMAP = Int32(-35)
const EGADS_SEQUERR = Int32(-34)
const EGADS_CNTXTHRD = Int32(-33)
const EGADS_READERR = Int32(-32)
const EGADS_TESSTATE = Int32(-31)
const EGADS_EXISTS = Int32(-30)
const EGADS_ATTRERR = Int32(-29)
const EGADS_TOPOCNT = Int32(-28)
const EGADS_OCSEGFLT = Int32(-27)
const EGADS_BADSCALE = Int32(-26)
const EGADS_NOTORTHO = Int32(-25)
const EGADS_DEGEN = Int32(-24)
const EGADS_CONSTERR = Int32(-23)
const EGADS_TOPOERR = Int32(-22)
const EGADS_GEOMERR = Int32(-21)
const EGADS_NOTBODY = Int32(-20)
const EGADS_WRITERR = Int32(-19)
const EGADS_NOTMODEL = Int32(-18)
const EGADS_NOLOAD = Int32(-17)
const EGADS_RANGERR = Int32(-16)
const EGADS_NOTGEOM = Int32(-15)
const EGADS_NOTTESS = Int32(-14)
const EGADS_EMPTY = Int32(-13)
const EGADS_NOTTOPO = Int32(-12)
const EGADS_REFERCE = Int32(-11)
const EGADS_NOTXFORM = Int32(-10)
const EGADS_NOTCNTX = Int32(-9)
const EGADS_MIXCNTX = Int32(-8)
const EGADS_NODATA = Int32(-7)
const EGADS_NONAME = Int32(-6)
const EGADS_INDEXERR = Int32(-5)
const EGADS_MALLOC = Int32(-4)
const EGADS_NOTOBJ = Int32(-3)
const EGADS_NULLOBJ = Int32(-2)
const EGADS_NOTFOUND = Int32(-1)
const EGADS_SUCCESS = Int32(0)
const EGADS_OUTSIDE = Int32(1)
const CHUNK = Int32(256)
const MAXELEN = Int32(2048)
const DEGENUV = 1.0e-13

const EGADSMAJOR = Int32(1)
const EGADSMINOR = Int32(21)
const EGADSPROP ="  EGADSprop:Revision(1.21)" 
const MAGIC = Int32(98789)
const MTESSPARAM = Int32(2)
const CONTXT = Int32(0)
const TRANSFORM = Int32(1)
const TESSELLATION = Int32(2)
const NIL = Int32(3)
const EMPTY = Int32(4)
const REFERENCE = Int32(5)
const PCURVE = Int32(10)
const CURVE = Int32(11)
const SURFACE = Int32(12)
const NODE = Int32(20)
const EDGE = Int32(21)
const LOOP = Int32(22)
const FACE = Int32(23)
const SHELL = Int32(24)
const BODY = Int32(25)
const MODEL = Int32(26)
const EEDGE = Int32(31)
const ELOOPX = Int32(32)
const EFACE = Int32(33)
const ESHELL = Int32(34)
const EBODY = Int32(35)
const LINE = Int32(1)
const CIRCLE = Int32(2)
const ELLIPSE = Int32(3)
const PARABOLA = Int32(4)
const HYPERBOLA = Int32(5)
const TRIMMED = Int32(6)
const BEZIER = Int32(7)
const BSPLINE = Int32(8)
const OFFSET = Int32(9)
const PLANE = Int32(1)
const SPHERICAL = Int32(2)
const CYLINDRICAL = Int32(3)
const REVOLUTION = Int32(4)
const TOROIDAL = Int32(5)
const CONICAL = Int32(10)
const EXTRUSION = Int32(11)
const SREVERSE = Int32(-1)
const SINNER = Int32(-1)
const NOMTYPE = Int32(0)
const SFORWARD = Int32(1)
const SOUTER = Int32(1)
const ONENODE = Int32(1)
const TWONODE = Int32(2)
const OPEN = Int32(3)
const CLOSED = Int32(4)
const DEGENERATE = Int32(5)
const WIREBODY = Int32(6)
const FACEBODY = Int32(7)
const SHEETBODY = Int32(8)
const SOLIDBODY = Int32(9)
const ATTRINT = Int32(1)
const ATTRREAL = Int32(2)
const ATTRSTRING = Int32(3)
const ATTRCSYS = Int32(12)
const ATTRPTR = Int32(13)
const SUBTRACTION = Int32(1)
const INTERSECTION = Int32(2)
const FUSION = Int32(3)
const SPLITTER = Int32(4)
const BOX = Int32(1)
const SPHERE = Int32(2)
const CONE = Int32(3)
const CYLINDER = Int32(4)
const TORUS = Int32(5)
const UISO = Int32(0)
const VISO = Int32(1)
const NODEOFF = Int32(1)
const EDGEOFF = Int32(2)
const FACEDUP = Int32(3)
const FACECUT = Int32(4)
const FACEOFF = Int32(5)
