###########################################################################
#                                                                         #
# pyEGADS --- Python version of EGADS API                                 #
#                                                                         #
#                                                                         #
#      Copyright 2011-2022, Massachusetts Institute of Technology         #
#      Licensed under The GNU Lesser General Public License, version 2.1  #
#      See http://www.opensource.org/licenses/lgpl-2.1.php                #
#                                                                         #
#   This API was inspired by the egads4py package:                        #
#      https://github.com/smdogroup/egads4py                              #
#                                                                         #
###########################################################################

import ctypes
from ctypes import POINTER, c_int, c_size_t, c_double, c_void_p, c_char, c_char_p
import os
import sys
import weakref
import itertools
from sys import version_info

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

# load the shared library
if sys.platform.startswith('darwin'):
    _egads = ctypes.CDLL(_ESP_ROOT + "/lib/libegads.dylib")
elif sys.platform.startswith('linux'):
    _egads = ctypes.CDLL(_ESP_ROOT + "/lib/libegads.so")
elif sys.platform.startswith('win32'):
    if version_info.major == 3 and version_info.minor < 8:
        _egads = ctypes.CDLL(_ESP_ROOT + "\\lib\\egads.dll")
    else:
        _egads = ctypes.CDLL(_ESP_ROOT + "\\lib\\egads.dll", winmode=0)
else:
    raise IOError("Unknown platform: " + sys.platform)

__all__ = ["revision", "Context", "ego", "c_ego",
           "csystem", "c_to_py", "EGADSError"]

# =============================================================================
# Taken from weakref to be used with older versions of Python
if sys.version_info[:2] < (3, 4):
    class finalize:
        """Class for finalization of weakrefable objects

        finalize(obj, func, *args, **kwargs) returns a callable finalizer
        object which will be called when obj is garbage collected. The
        first time the finalizer is called it evaluates func(*arg, **kwargs)
        and returns the result. After this the  is dead, and
        calling it just returns None.

        When the program exits any remaining finalizers for which the
        atexit attribute is true will be run in reverse order of creation.
        By default atexit is true.
        """

        # Finalizer objects don't have any state of their own.  They are
        # just used as keys to lookup _Info objects in the registry.  This
        # ensures that they cannot be part of a ref-cycle.

        __slots__ = ()
        _registry = {}
        _shutdown = False
        _index_iter = itertools.count()
        _dirty = False
        _registered_with_atexit = False

        class _Info:
            __slots__ = ("weakref", "func", "args", "kwargs", "atexit", "index")

        def __init__(self, obj, func, *args, **kwargs):
            if not self._registered_with_atexit:
                # We may register the exit function more than once because
                # of a thread race, but that is harmless
                import atexit
                atexit.register(self._exitfunc)
                finalize._registered_with_atexit = True
            info = self._Info()
            info.weakref = weakref.ref(obj, self)
            info.func = func
            info.args = args
            info.kwargs = kwargs or None
            info.atexit = True
            info.index = next(self._index_iter)
            self._registry[self] = info
            finalize._dirty = True

        def __call__(self, _=None):
            """If alive then mark as dead and return func(*args, **kwargs);
            otherwise return None"""
            info = self._registry.pop(self, None)
            if info and not self._shutdown:
                return info.func(*info.args, **(info.kwargs or {}))

        def detach(self):
            """If alive then mark as dead and return (obj, func, args, kwargs);
            otherwise return None"""
            info = self._registry.get(self)
            obj = info and info.weakref()
            if obj is not None and self._registry.pop(self, None):
                return (obj, info.func, info.args, info.kwargs or {})

        def peek(self):
            """If alive then return (obj, func, args, kwargs);
            otherwise return None"""
            info = self._registry.get(self)
            obj = info and info.weakref()
            if obj is not None:
                return (obj, info.func, info.args, info.kwargs or {})

        @property
        def alive(self):
            """Whether finalizer is alive"""
            return self in self._registry

        @property
        def atexit(self):
            """Whether finalizer should be called at exit"""
            info = self._registry.get(self)
            return bool(info) and info.atexit

        @atexit.setter
        def atexit(self, value):
            info = self._registry.get(self)
            if info:
                info.atexit = bool(value)

        def __repr__(self):
            info = self._registry.get(self)
            obj = info and info.weakref()
            if obj is None:
                return '<%s object at %#x; dead>' % (type(self).__name__, id(self))
            else:
                return '<%s object at %#x; for %r at %#x>' % \
                    (type(self).__name__, id(self), type(obj).__name__, id(obj))

        @classmethod
        def _select_for_exit(cls):
            # Return live finalizers marked for exit, oldest first
            L = [(f,i) for (f,i) in cls._registry.items() if i.atexit]
            L.sort(key=lambda item:item[1].index)
            return [f for (f,i) in L]

        @classmethod
        def _exitfunc(cls):
            # At shutdown invoke finalizers for which atexit is true.
            # This is called once all other non-daemonic threads have been
            # joined.
            reenable_gc = False
            try:
                if cls._registry:
                    import gc
                    if gc.isenabled():
                        reenable_gc = True
                        gc.disable()
                    pending = None
                    while True:
                        if pending is None or finalize._dirty:
                            pending = cls._select_for_exit()
                            finalize._dirty = False
                        if not pending:
                            break
                        f = pending.pop()
                        try:
                            # gc is disabled, so (assuming no daemonic
                            # threads) the following is the only line in
                            # this function which might trigger creation
                            # of a new finalizer
                            f()
                        except Exception:
                            sys.excepthook(*sys.exc_info())
                        assert f not in cls._registry
            finally:
                # prevent any more finalizers from executing during shutdown
                finalize._shutdown = True
                if reenable_gc:
                    gc.enable()
else:
    finalize = weakref.finalize
    
# =============================================================================
# Decode function to play nice with Pyhon 2.7
def _decode(data):
    if version_info.major > 2 and isinstance(data, bytes):
        return data.decode()
    return data

# =============================================================================
# define union for c_egAttr structure
class union_egAttr(ctypes.Union):
    pass

union_egAttr.__slots__ = [
    'integer',
    'integers',
    'real',
    'reals',
    'string',
]
union_egAttr._fields_ = [
    ('integer', c_int),
    ('integers', POINTER(c_int)),
    ('real', c_double),
    ('reals', POINTER(c_double)),
    ('string', c_char_p),
]

# =============================================================================
# define c_egAttr structure
class c_egAttr(ctypes.Structure):
    pass

c_egAttr.__slots__ = [
    'name',
    'type',
    'length',
    'vals',
]
c_egAttr._fields_ = [
    ('name', c_char_p),
    ('type', c_int),
    ('length', c_int),
    ('vals', union_egAttr),
]

# =============================================================================
# define c_egAttrSeq structure
class c_egAttrSeq(ctypes.Structure):
    pass

c_egAttrSeq.__slots__ = [
    'root',
    'nSeq',
    'attrSeq',
]
c_egAttrSeq._fields_ = [
    ('root', c_char_p),
    ('nSeq', c_int),
    ('attrSeq', POINTER(c_int)),
]

# =============================================================================
# define c_egAttrSeq structure
class c_egAttrs(ctypes.Structure):
    pass

c_egAttrs.__slots__ = [
    'nattrs',
    'attrs',
    'nseqs',
    'seqs',
]
c_egAttrs._fields_ = [
    ('nattrs', c_int),
    ('attrs', POINTER(c_egAttr)),
    ('nseqs', c_int),
    ('seqs', POINTER(c_egAttrSeq)),
]

# =============================================================================
# define c_egObject structure
class c_egObject(ctypes.Structure):
    pass

c_egObject.__slots__ = [
    'magicnumber',
    'oclass',
    'mtype',
    'attrs',
    'blind',
    'topObj',
    'tref',
    'prev',
    'next',
]

c_egObject._fields_ = [
    ("magicnumber", ctypes.c_int),           # must be set to validate the object
    ("oclass", ctypes.c_short),              # object Class
    ("mtype", ctypes.c_short),               # member Type
    ("attrs", ctypes.POINTER(None)),         # object Attributes or Reference
    ("blind", ctypes.POINTER(None)),         # blind pointer to object data
    ("topObj", ctypes.POINTER(c_egObject)),  # top of the hierarchy or context (if top)
    ("tref", ctypes.POINTER(c_egObject)),    # threaded list of references
    ("prev", ctypes.POINTER(c_egObject)),    # back pointer
    ("next", ctypes.POINTER(c_egObject)),    # forward pointer
]

# define the c_ego (ego) pointer
c_ego = ctypes.POINTER(c_egObject)

# =============================================================================
#
# egads.h functions
#
# =============================================================================

# egads functions declarations
_egads.EG_free.argtypes = [c_void_p]
_egads.EG_free.restype = None

_egads.EG_revision.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_char_p)]
_egads.EG_revision.restype = None

_egads.EG_open.argtypes = [POINTER(c_ego)]
_egads.EG_open.restype = c_int

_egads.EG_loadModel.argtypes = [c_ego, c_int, c_char_p, POINTER(c_ego)]
_egads.EG_loadModel.restype = c_int

_egads.EG_saveModel.argtypes = [c_ego, c_char_p]
_egads.EG_saveModel.restype = c_int

_egads.EG_exportModel.argtypes = [c_ego, POINTER(c_size_t), POINTER(c_char_p)]
_egads.EG_exportModel.restype = c_int

_egads.EG_deleteObject.argtypes = [c_ego]
_egads.EG_deleteObject.restype = c_int

_egads.EG_makeTransform.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
_egads.EG_makeTransform.restype = c_int

_egads.EG_getTransformation.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_getTransformation.restype = c_int

_egads.EG_getContext.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_getContext.restype = c_int

_egads.EG_setOutLevel.argtypes = [c_ego, c_int]
_egads.EG_setOutLevel.restype = c_int

_egads.EG_updateThread.argtypes = [c_ego]
_egads.EG_updateThread.restype = c_int

_egads.EG_getInfo.argtypes = [c_ego, POINTER(c_int), POINTER(c_int), POINTER(c_ego), POINTER(c_ego), POINTER(c_ego)]
_egads.EG_getInfo.restype = c_int

_egads.EG_copyObject.argtypes = [c_ego, POINTER(None), POINTER(c_ego)]
_egads.EG_copyObject.restype = c_int

_egads.EG_flipObject.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_flipObject.restype = c_int

_egads.EG_close.argtypes = [c_ego]
_egads.EG_close.restype = c_int

_egads.EG_setUserPointer.argtypes = [c_ego, POINTER(None)]
_egads.EG_setUserPointer.restype = c_int

_egads.EG_getUserPointer.argtypes = [c_ego, POINTER(POINTER(None))]
_egads.EG_getUserPointer.restype = c_int

_egads.EG_setFullAttrs.argtypes = [c_ego, c_int]
_egads.EG_setFullAttrs.restype = c_int

_egads.EG_attributeAdd.argtypes = [c_ego, c_char_p, c_int, c_int, POINTER(c_int), POINTER(c_double), c_char_p]
_egads.EG_attributeAdd.restype = c_int

_egads.EG_attributeDel.argtypes = [c_ego, c_char_p]
_egads.EG_attributeDel.restype = c_int

_egads.EG_attributeNum.argtypes = [c_ego, POINTER(c_int)]
_egads.EG_attributeNum.restype = c_int

_egads.EG_attributeGet.argtypes = [c_ego, c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egads.EG_attributeGet.restype = c_int

_egads.EG_attributeRet.argtypes = [c_ego, c_char_p, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egads.EG_attributeRet.restype = c_int

_egads.EG_attributeDup.argtypes = [c_ego, c_ego]
_egads.EG_attributeDup.restype = c_int

_egads.EG_attributeAddSeq.argtypes = [c_ego, c_char_p, c_int, c_int, POINTER(c_int), POINTER(c_double), c_char_p]
_egads.EG_attributeAddSeq.restype = c_int

_egads.EG_attributeNumSeq.argtypes = [c_ego, c_char_p, POINTER(c_int)]
_egads.EG_attributeNumSeq.restype = c_int

_egads.EG_attributeRetSeq.argtypes = [c_ego, c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egads.EG_attributeRetSeq.restype = c_int

_egads.EG_getGeometry.argtypes = [c_ego, POINTER(c_int), POINTER(c_int), POINTER(c_ego), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
_egads.EG_getGeometry.restype = c_int

_egads.EG_getGeometryLen.argtypes = [c_ego, POINTER(c_int), POINTER(c_int)]
_egads.EG_getGeometryLen.restype = None

_egads.EG_makeGeometry.argtypes = [c_ego, c_int, c_int, c_ego, POINTER(c_int), POINTER(c_double), POINTER(c_ego)]
_egads.EG_makeGeometry.restype = c_int

_egads.EG_getRange.argtypes = [c_ego, POINTER(c_double), POINTER(c_int)]
_egads.EG_getRange.restype = c_int

_egads.EG_evaluate.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
_egads.EG_evaluate.restype = c_int

_egads.EG_invEvaluate.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_invEvaluate.restype = c_int

_egads.EG_invEvaluateGuess.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_invEvaluateGuess.restype = c_int

_egads.EG_arcLength.argtypes = [c_ego, c_double, c_double, POINTER(c_double)]
_egads.EG_arcLength.restype = c_int

_egads.EG_curvature.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
_egads.EG_curvature.restype = c_int

_egads.EG_approximate.argtypes = [c_ego, c_int, c_double, POINTER(c_int), POINTER(c_double), POINTER(c_ego)]
_egads.EG_approximate.restype = c_int

_egads.EG_fitTriangles.argtypes = [c_ego, c_int, POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), c_double, POINTER(c_ego)]
_egads.EG_fitTriangles.restype = c_int

_egads.EG_otherCurve.argtypes = [c_ego, c_ego, c_double, POINTER(c_ego)]
_egads.EG_otherCurve.restype = c_int

_egads.EG_isSame.argtypes = [c_ego, c_ego]
_egads.EG_isSame.restype = c_int

_egads.EG_isoCline.argtypes = [c_ego, c_int, c_double, POINTER(c_ego)]
_egads.EG_isoCline.restype = c_int

_egads.EG_convertToBSpline.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_convertToBSpline.restype = c_int

_egads.EG_convertToBSplineRange.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
_egads.EG_convertToBSplineRange.restype = c_int

_egads.EG_skinning.argtypes = [c_int, POINTER(c_ego), c_int, POINTER(c_ego)]
_egads.EG_skinning.restype = c_int

_egads.EG_tolerance.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_tolerance.restype = c_int

_egads.EG_getTolerance.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_getTolerance.restype = c_int

_egads.EG_setTolerance.argtypes = [c_ego, c_double]
_egads.EG_setTolerance.restype = c_int

_egads.EG_getTopology.argtypes = [c_ego, POINTER(c_ego), POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(POINTER(c_int))]
_egads.EG_getTopology.restype = c_int

_egads.EG_makeTopology.argtypes = [c_ego, c_ego, c_int, c_int, POINTER(c_double), c_int, POINTER(c_ego), POINTER(c_int), POINTER(c_ego)]
_egads.EG_makeTopology.restype = c_int

_egads.EG_makeLoop.argtypes = [c_int, POINTER(c_ego), c_ego, c_double, POINTER(c_ego)]
_egads.EG_makeLoop.restype = c_int

_egads.EG_getArea.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
_egads.EG_getArea.restype = c_int

_egads.EG_makeFace.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_ego)]
_egads.EG_makeFace.restype = c_int

_egads.EG_getBodyTopos.argtypes = [c_ego, c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_ego))]
_egads.EG_getBodyTopos.restype = c_int

_egads.EG_indexBodyTopo.argtypes = [c_ego, c_ego]
_egads.EG_indexBodyTopo.restype = c_int

_egads.EG_objectBodyTopo.argtypes = [c_ego, c_int, c_int, POINTER(c_ego)]
_egads.EG_objectBodyTopo.restype = c_int

_egads.EG_inTopology.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_inTopology.restype = c_int

_egads.EG_inFace.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_inFace.restype = c_int

_egads.EG_getEdgeUV.argtypes = [c_ego, c_ego, c_int, c_double, POINTER(c_double)]
_egads.EG_getEdgeUV.restype = c_int

_egads.EG_getEdgeUVs.argtypes = [c_ego, c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double)]
_egads.EG_getEdgeUVs.restype = c_int

_egads.EG_getPCurve.argtypes = [c_ego, c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
_egads.EG_getPCurve.restype = c_int

_egads.EG_getWindingAngle.argtypes = [c_ego, c_double, POINTER(c_double)]
_egads.EG_getWindingAngle.restype = c_int

_egads.EG_getBody.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_getBody.restype = c_int

_egads.EG_makeSolidBody.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_ego)]
_egads.EG_makeSolidBody.restype = c_int

_egads.EG_getBoundingBox.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_getBoundingBox.restype = c_int

_egads.EG_getMassProperties.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_getMassProperties.restype = c_int

_egads.EG_isEquivalent.argtypes = [c_ego, c_ego]
_egads.EG_isEquivalent.restype = c_int

_egads.EG_sewFaces.argtypes = [c_int, POINTER(c_ego), c_double, c_int, POINTER(c_ego)]
_egads.EG_sewFaces.restype = c_int

_egads.EG_replaceFaces.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_replaceFaces.restype = c_int

_egads.EG_mapBody.argtypes = [c_ego, c_ego, c_char_p, POINTER(c_ego)]
_egads.EG_mapBody.restype = c_int

_egads.EG_matchBodyEdges.argtypes = [c_ego, c_ego, c_double, POINTER(c_int), POINTER(POINTER(c_int))]
_egads.EG_matchBodyEdges.restype = c_int

_egads.EG_matchBodyFaces.argtypes = [c_ego, c_ego, c_double, POINTER(c_int), POINTER(POINTER(c_int))]
_egads.EG_matchBodyFaces.restype = c_int

_egads.EG_setTessParam.argtypes = [c_ego, c_int, c_double, POINTER(c_double)]
_egads.EG_setTessParam.restype = c_int

_egads.EG_makeTessGeom.argtypes = [c_ego, POINTER(c_double), POINTER(c_int), POINTER(c_ego)]
_egads.EG_makeTessGeom.restype = c_int

_egads.EG_getTessGeom.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_double))]
_egads.EG_getTessGeom.restype = c_int

_egads.EG_makeTessBody.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
_egads.EG_makeTessBody.restype = c_int

_egads.EG_remakeTess.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_double)]
_egads.EG_remakeTess.restype = c_int

_egads.EG_finishTess.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_finishTess.restype = c_int

_egads.EG_mapTessBody.argtypes = [c_ego, c_ego, POINTER(c_ego)]
_egads.EG_mapTessBody.restype = c_int

_egads.EG_locateTessBody.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(c_double)]
_egads.EG_locateTessBody.restype = c_int

_egads.EG_getTessEdge.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double))]
_egads.EG_getTessEdge.restype = c_int

_egads.EG_getTessFace.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_int))]
_egads.EG_getTessFace.restype = c_int

_egads.EG_getTessLoops.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_int))]
_egads.EG_getTessLoops.restype = c_int

_egads.EG_getTessQuads.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_int))]
_egads.EG_getTessQuads.restype = c_int

_egads.EG_makeQuads.argtypes = [c_ego, POINTER(c_double), c_int]
_egads.EG_makeQuads.restype = c_int

_egads.EG_getQuads.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
_egads.EG_getQuads.restype = c_int

_egads.EG_getPatch.argtypes = [c_ego, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_int))]
_egads.EG_getPatch.restype = c_int

_egads.EG_quadTess.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_quadTess.restype = c_int

_egads.EG_insertEdgeVerts.argtypes = [c_ego, c_int, c_int, c_int, POINTER(c_double)]
_egads.EG_insertEdgeVerts.restype = c_int

_egads.EG_deleteEdgeVert.argtypes = [c_ego, c_int, c_int, c_int]
_egads.EG_deleteEdgeVert.restype = c_int

_egads.EG_moveEdgeVert.argtypes = [c_ego, c_int, c_int, c_double]
_egads.EG_moveEdgeVert.restype = c_int

_egads.EG_openTessBody.argtypes = [c_ego]
_egads.EG_openTessBody.restype = c_int

_egads.EG_initTessBody.argtypes = [c_ego, POINTER(c_ego)]
_egads.EG_initTessBody.restype = c_int

_egads.EG_statusTessBody.argtypes = [c_ego, POINTER(c_ego), POINTER(c_int), POINTER(c_int)]
_egads.EG_statusTessBody.restype = c_int

_egads.EG_setTessEdge.argtypes = [c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double)]
_egads.EG_setTessEdge.restype = c_int

_egads.EG_setTessFace.argtypes = [c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int)]
_egads.EG_setTessFace.restype = c_int

_egads.EG_localToGlobal.argtypes = [c_ego, c_int, c_int, POINTER(c_int)]
_egads.EG_localToGlobal.restype = c_int

_egads.EG_getGlobal.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double)]
_egads.EG_getGlobal.restype = c_int

_egads.EG_saveTess.argtypes = [c_ego, c_char_p]
_egads.EG_saveTess.restype = c_int

_egads.EG_loadTess.argtypes = [c_ego, c_char_p, POINTER(c_ego)]
_egads.EG_loadTess.restype = c_int

_egads.EG_tessMassProps.argtypes = [c_ego, POINTER(c_double)]
_egads.EG_tessMassProps.restype = c_int

_egads.EG_fuseSheets.argtypes = [c_ego, c_ego, POINTER(c_ego)]
_egads.EG_fuseSheets.restype = c_int

_egads.EG_generalBoolean.argtypes = [c_ego, c_ego, c_int, c_double, POINTER(c_ego)]
_egads.EG_generalBoolean.restype = c_int

_egads.EG_solidBoolean.argtypes = [c_ego, c_ego, c_int, POINTER(c_ego)]
_egads.EG_solidBoolean.restype = c_int

_egads.EG_intersection.argtypes = [c_ego, c_ego, POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(c_ego)]
_egads.EG_intersection.restype = c_int

_egads.EG_imprintBody.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_imprintBody.restype = c_int

_egads.EG_filletBody.argtypes = [c_ego, c_int, POINTER(c_ego), c_double, POINTER(c_ego), POINTER(POINTER(c_int))]
_egads.EG_filletBody.restype = c_int

_egads.EG_chamferBody.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego), c_double, c_double, POINTER(c_ego), POINTER(POINTER(c_int))]
_egads.EG_chamferBody.restype = c_int

_egads.EG_hollowBody.argtypes = [c_ego, c_int, POINTER(c_ego), c_double, c_int, POINTER(c_ego), POINTER(POINTER(c_int))]
_egads.EG_hollowBody.restype = c_int

_egads.EG_extrude.argtypes = [c_ego, c_double, POINTER(c_double), POINTER(c_ego)]
_egads.EG_extrude.restype = c_int

_egads.EG_rotate.argtypes = [c_ego, c_double, POINTER(c_double), POINTER(c_ego)]
_egads.EG_rotate.restype = c_int

_egads.EG_sweep.argtypes = [c_ego, c_ego, c_int, POINTER(c_ego)]
_egads.EG_sweep.restype = c_int

_egads.EG_loft.argtypes = [c_int, POINTER(c_ego), c_int, POINTER(c_ego)]
_egads.EG_loft.restype = c_int

_egads.EG_blend.argtypes = [c_int, POINTER(c_ego), POINTER(c_double), POINTER(c_double), POINTER(c_ego)]
_egads.EG_blend.restype = c_int

_egads.EG_ruled.argtypes = [c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_ruled.restype = c_int

_egads.EG_initEBody.argtypes = [c_ego, c_double, POINTER(c_ego)]
_egads.EG_initEBody.restype = c_int

_egads.EG_finishEBody.argtypes = [c_ego]
_egads.EG_finishEBody.restype = c_int

_egads.EG_makeEFace.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_makeEFace.restype = c_int

_egads.EG_makeAttrEFaces.argtypes = [c_ego, c_char_p, POINTER(c_int), POINTER(POINTER(c_ego))]
_egads.EG_makeAttrEFaces.restype = c_int

_egads.EG_effectiveMap.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego), POINTER(c_double)]
_egads.EG_effectiveMap.restype = c_int

_egads.EG_effectiveEdgeList.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
_egads.EG_effectiveEdgeList.restype = c_int

_egads.EG_inTriExact.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_inTriExact.restype = c_int

# =============================================================================
#
# egads_dot.h functions
#
# =============================================================================

_egads.EG_setGeometry_dot.argtypes = [c_ego, c_int, c_int, POINTER(c_int), POINTER(c_double), POINTER(c_double)]
_egads.EG_setGeometry_dot.restype = c_int

_egads.EG_hasGeometry_dot.argtypes = [c_ego]
_egads.EG_hasGeometry_dot.restype = c_int

_egads.EG_copyGeometry_dot.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), c_ego]
_egads.EG_copyGeometry_dot.restype = c_int

_egads.EG_evaluate_dot.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_evaluate_dot.restype = c_int

_egads.EG_approximate_dot.argtypes = [c_ego, c_int, c_double, POINTER(c_int), POINTER(c_double), POINTER(c_double)]
_egads.EG_approximate_dot.restype = c_int

_egads.EG_skinning_dot.argtypes = [c_ego, c_int, POINTER(c_ego)]
_egads.EG_skinning_dot.restype = c_int

_egads.EG_makeSolidBody_dot.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_double)]
_egads.EG_makeSolidBody_dot.restype = c_int

_egads.EG_setRange_dot.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_double)]
_egads.EG_setRange_dot.restype = c_int

_egads.EG_getRange_dot.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_int)]
_egads.EG_getRange_dot.restype = c_int

_egads.EG_tessMassProps_dot.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_tessMassProps_dot.restype = c_int

_egads.EG_extrude_dot.argtypes = [c_ego, c_ego, c_double, c_double, POINTER(c_double), POINTER(c_double)]
_egads.EG_extrude_dot.restype = c_int

_egads.EG_ruled_dot.argtypes = [c_ego, c_int, POINTER(c_ego)]
_egads.EG_ruled_dot.restype = c_int

_egads.EG_blend_dot.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egads.EG_blend_dot.restype = c_int

# =============================================================================

# Extract EGADS error codes
_egads_error_codes = {}
globalDict = globals()
with open(os.path.join(_ESP_ROOT,"include","egadsErrors.h")) as fp:
    lines = fp.readlines()
    for line in lines:
        define = line.split()
        if len(define) < 3: continue
        if define[0] != "#define": continue
        _egads_error_codes[int(define[2])] = define[1]
        globalDict[define[1]] = int(define[2])
        __all__.append(define[1])
del fp
del lines

# Extract PRM error codes
with open(os.path.join(_ESP_ROOT,"include","prm.h")) as fp:
    lines = fp.readlines()
    for line in lines:
        define = line.split()
        if len(define) < 3: continue
        if define[0] != "#define": continue
        if int(define[2]) < 0:
            _egads_error_codes[int(define[2])] = define[1]
        globalDict[define[1]] = int(define[2])
        __all__.append(define[1])
del fp
del lines

# Extract EGADS constants and them to the globals() dictionary
with open(os.path.join(_ESP_ROOT,"include","egadsTypes.h")) as fp:
    lines = fp.readlines()
    for line in lines:
        define = line.split()
        if len(define) < 3: continue
        if define[0] != "#define": continue
        if define[1] == "EGADSPROP":
            globalDict[define[1]] = ' '.join(define[2:])
            continue
        globalDict[define[1]] = int(define[2])
        __all__.append(define[1])
del fp
del lines
del globalDict

# =============================================================================

class EGADSError(Exception):

    # Constructor or Initializer
    def __init__(self, string, status = None):
        self.string = string
        self.status = status

    def __str__(self):
        return self.string

# =============================================================================

def _raiseStatus(stat):
    if stat == EGADS_SUCCESS: return
    errmsg = 'EGADS error code: %s'%(_egads_error_codes[stat])
    raise EGADSError(errmsg, stat)

# =============================================================================

def revision():
    """
    Revision - return EGADS revision

    Returns
    -------
        imajor:  major version number
        iminor:  minor version number
        OCCrev:  OpenCASCADE version string
    """
    imajor = c_int()
    iminor = c_int()
    string = ctypes.c_char_p()

    stat = _egads.EG_revision(ctypes.byref(imajor), ctypes.byref(iminor), ctypes.byref(string))
    if stat: _raiseStatus(stat)

    return (imajor.value, iminor.value, _decode(string.value))

# =============================================================================

def free(ptr):
    """
    Free memory allocated with EG_alloc or EG_realloc returned via ctypes

    Parameters
    ----------
    ptr:
        the ctypes memory poimter to deallocate
    """
    _egads.EG_free(ptr)

# =============================================================================

def c_to_py(c_obj, deleteObject=False, context=None):
    """
    Creates a Python class for an exiting c_ego object.
    
    Parameters
    ----------
    c_obj:
        the c_ego instance
        
    deleteObject:
        if True the class instance will call EG_deleteObject or EG_close 
        for the provided c_ego during garbage collection

    context:
        a Context instance associated with c_obj
    """
    if not isinstance(c_obj, c_ego):
        raise EGADSError('c_obj must be a c_ego instance')

    if context is not None:
        if not isinstance(context, Context):
            raise EGADSError('context must be a Context instance')

        c_context = c_ego()
        stat = _egads.EG_getContext(c_obj, ctypes.byref(c_context))
        if stat: _raiseStatus(stat)
        if ctypes.addressof(context._context.contents) != ctypes.addressof(c_context.contents):
            raise EGADSError("context inconsistent with c_obj!")

    if c_obj.contents.oclass == CONTXT:
        return Context(c_obj, deleteObject=deleteObject)
    else:
        return ego(c_obj, context, deleteObject=deleteObject)

# =============================================================================

class csystem(list):
    """
    Represents a coordinate system attribute
    """
    pass

# =============================================================================

def _close_context(context):
    stat = _egads.EG_close(context)
    if stat: _raiseStatus(stat)

class Context:
    """
    Wrapper to represent and EGADS CONTEXT ego
    """
#=============================================================================-
    def __init__(self, c_context=None, deleteObject=False):
        """
        Constructor only intended for internal use.

        Parameters
        ----------
        c_context:
            a c_ego egads context

        deleteObject:
            If True c_context will be closed automatically
        """
        self._finalize = None
        self._context = None
        if c_context is not None:
            if not isinstance(c_context, c_ego):
                errmsg = 'context must be a c_ego'
                raise EGADSError(errmsg)

            if c_context.contents.oclass != CONTXT:
                errmsg = 'context must be an EGADS CONTXT'
                raise EGADSError(errmsg)

            self._context = c_context
            if deleteObject:
                self._finalize = finalize(self, _close_context, self._context)
        else:
            self._context = c_ego()
            stat = _egads.EG_open(ctypes.byref(self._context))
            if stat: _raiseStatus(stat)
            self._finalize = finalize(self, _close_context, self._context)

#=============================================================================-
    def py_to_c(self, takeOwnership=False):
        """
        Returns the c_ego pointer
        
        Parameters
        ----------
        takeOwnership:
            if True, pyEGADS will no longer close the returned c_ego Context during garbage collection
        """
        if takeOwnership:
            self._finalize.detach()
            self._finalize = None
        return self._context

#=============================================================================-
    def __del__(self):
        # close the context if owned by this instance
        if self._finalize is not None:
            self._finalize()

#=============================================================================-
    def __eq__(self, obj):
        """
        checks pointer equality of underlying c_ego memory address
        """
        return isinstance(obj, Context) and ctypes.addressof(obj._context.contents) == ctypes.addressof(self._context.contents)

#=============================================================================-
    def __ne__(self, obj):
        return not self == obj

#=============================================================================-
    def setOutLevel(self, outlevel=1):
        """
        Sets the EGADS verbosity level, the default is 1.

        Parameters
        ----------
        outlevel: 0 <= outlevel <= 3
            0-silent to 3-debug

        Returns
        -------
        The old outLevel
        """
        stat = _egads.EG_setOutLevel(self._context, c_int(outlevel))
        if stat < 0: _raiseStatus(stat)
        return stat

#=============================================================================-
    def setFullAttrs(self, attrFlag):
        """
        Sets the attribution mode for the Context.

        Parameters
        ----------
        attrFlag:
            the mode flag: 0 - the default scheme, 1 - full attribution node
        """
        stat = _egads.EG_setFullAttrs(self._context, c_int(attrFlag))
        if stat: _raiseStatus(stat)

#=============================================================================-
    def updateThread(self):
        """
        Resets the Context's owning thread ID to the thread calling this function.
        """
        stat = _egads.EG_updateThread(self._context)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def getInfo(self):
        """
        Return information about the context

        Returns
        -------
        oclass:
            object class type

        mtype:
            object sub-type

        topRef:
            is the top level BODY/MODEL that owns the object or context (if top)

        prev:
            is the previous object in the threaded list (None at CONTEXT)

        next:
            is the next object in the list (None is the end of the list)
        """
        oclass = c_int()
        mtype  = c_int()
        topObj = c_ego()
        prev   = c_ego()
        next   = c_ego()
        stat = _egads.EG_getInfo(self._context, ctypes.byref(oclass), ctypes.byref(mtype),
                          ctypes.byref(topObj), ctypes.byref(prev), ctypes.byref(next))
        if stat: _raiseStatus(stat)
        return oclass.value, mtype.value, ego(topObj, self), ego(prev, self), ego(next, self)

#=============================================================================-
    def makeTransform(self, xform):
        """
        Makes a Transformation Object from a [3][4] translation/rotation matrix.
        The rotation portion [3][3] must be "scaled" orthonormal (orthogonal with a
        single scale factor).

        Returns
        -------
        resultant new transformation ego
        """
        T = (ctypes.c_double * 12)()
        for i in range(3):
            T[4*i:4*i+4] = xform[i]

        obj = c_ego()
        stat = _egads.EG_makeTransform(self._context, T, ctypes.byref(obj))
        if stat: _raiseStatus(stat)

        return ego(obj, self, deleteObject=True)

#=============================================================================-
    def makeGeometry(self, oclass, mtype, reals,
                     ints=None, geom=None):
        """
        Creates a geometric object:

        Parameters
        ----------
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
            content and length depends on the oclass/mtype.

        ints:
            is a pointer to the block of integer information. Required for
            either BEZIER or BSPLINE.

        geom:
            The reference geometry object (if none use None)

        Returns
        -------
        Resultant new geometry ego
        """
        if not (oclass == CURVE or oclass == PCURVE or oclass == SURFACE):
            errmsg = 'makeGeometry only accepts CURVE, PCURVE or SURFACE'
            raise ValueError(errmsg)

        ref = geom._obj if geom is not None else None

        ivec = None
        if ints is not None:
            ivec = (c_int * len(ints))()
            for i in range(len(ints)): ivec[i] = ints[i]

        rvec = None
        if reals is not None:
            if mtype == BSPLINE:

                if oclass == PCURVE or oclass == CURVE:
                    flag   = ints[0]
                    nCP    = ints[2]
                    nKnots = ints[3]

                    d = 2 if oclass == PCURVE else 3
                    nrvec = nKnots+d*nCP
                    if flag == 2: nrvec += nCP

                    rvec = (c_double * nrvec)()

                    rvec[:nKnots] = reals[0][:]
                    for i in range(nCP):
                        rvec[nKnots+d*i:nKnots+d*(i+1)] = reals[1][i][:]

                    if flag == 2:
                        rvec[nKnots+d*nCP:] = reals[2][:]
                else:
                    flag    = ints[0]
                    nCPu    = ints[2]
                    nUKnots = ints[3]
                    nCPv    = ints[5]
                    nVKnots = ints[6]
                    nKnots = nUKnots + nVKnots

                    nrvec = nKnots+3*nCPu*nCPv
                    if flag == 2: nrvec += nCPu*nCPv

                    rvec = (c_double * nrvec)()

                    rvec[      0:nUKnots] = reals[0][:]
                    rvec[nUKnots: nKnots] = reals[1][:]
                    for i in range(nCPu*nCPv):
                        rvec[nKnots+3*i:nKnots+3*(i+1)] = reals[2][i][:]

                    if flag == 2:
                        rvec[nKnots+3*nCPu*nCPv:] = reals[3][:]

            elif mtype == BEZIER:
                if oclass == PCURVE:
                    rvec = (c_double * (2*len(reals)))()
                    for i in range(len(reals)):
                        rvec[2*i:2*i+2] = reals[i]
                else:
                    rvec = (c_double * (3*len(reals)))()
                    for i in range(len(reals)):
                        rvec[3*i:3*i+3] = reals[i]
            else:
                rvec = (c_double * len(reals))()
                for i in range(len(reals)): rvec[i] = reals[i]

        obj = c_ego()
        stat = _egads.EG_makeGeometry(self._context, oclass, mtype, ref,
                                      ivec, rvec, ctypes.byref(obj))
        if stat: _raiseStatus(stat)

        return ego(obj, self, deleteObject=True, refs=geom)

#=============================================================================-
    def makeTopology(self, oclass, mtype=0, geom=None, reals=None,
                     children=None, senses=None):
        """
        Creates and returns a topological object:

        Parameters
        ----------
        oclass:
            either NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

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
            may be None except for:
            NODE which contains the [x,y,z] location
            EDGE is [t-min t-max] (the parametric bounds)
            FACE is [u-min, u-max, v-min, v-max] (the parametric bounds)

        children:
            chldrn a list of children objects (nchild in length)
            if LOOP and has reference SURFACE, then 2*nchild in length
            (PCURVES follow)

        senses:
            a list of integer senses for the children (required for FACEs
            & LOOPs only). For LOOPs, the senses are SFORWARD or SREVERSE
            for each EDGE. For FACEs, the senses are SOUTER or SINNER for
            (may be None for 1 child).

        Returns
        -------
        resultant new topological ego
        """

        errmsg = None
        if oclass == NODE:
            if len(reals) != 3:
                errmsg = 'NODE must be have data of length 3. len(data) = ' + str(len(data))
        if oclass == EDGE:
            if mtype != TWONODE and mtype != ONENODE and mtype != CLOSED:
                errmsg = 'EDGE must be TWONODE, ONENODE or CLOSED'
            if len(reals) != 2:
                errmsg = 'EDGE must be have data of length 2. len(data) = ' + str(len(data))
        elif oclass == LOOP:
            if mtype != OPEN and mtype != CLOSED:
                errmsg = 'LOOP must be OPEN or CLOSED'
        elif oclass == FACE:
            if mtype != SFORWARD and mtype != SREVERSE:
                errmsg = 'FACE must be SFORWARD or SREVERSE'
        elif oclass == SHELL:
            if mtype != OPEN and mtype != CLOSED:
                errmsg = 'SHELL must be OPEN or CLOSED'
        elif oclass == BODY:
            if (mtype != WIREBODY and mtype != FACEBODY and
                mtype != SHEETBODY and mtype != SOLIDBODY):
                errmsg = 'BODY must be WIREBODY, FACEBODY'
                errmsg += ', SHEETBODY or SOLIDBODY'

        if errmsg is not None:
            raise ValueError(errmsg)

        data = None if reals is None else (c_double * len(reals))()
        if oclass == NODE:
            data[0] = reals[0]
            data[1] = reals[1]
            data[2] = reals[2]
        elif oclass == EDGE:
            data[0] = reals[0]
            data[1] = reals[1]


        # Store references to children so they will be delted after the new ego
        refs = []

        nchildren = 0
        pchildren = None
        if children is not None:
            refs = list(children)
            nchildren = len(children)
            pchildren = (c_ego * nchildren)()
            for i in range(nchildren):
                pchildren[i] = children[i]._obj
        
        # Account for effective topology 
        if oclass == MODEL:
            mtype = nchildren
            nchildren = 0
            for child in children:
                if child._obj.contents.oclass == BODY:
                    nchildren += 1

        geo = None if geom is None else geom._obj
        if geom is not None: refs.append(geom)

        nsenses = 0
        psenses = None
        if oclass == LOOP or (oclass == FACE and nchildren > 1):
            nsenses = len(senses)
            psenses = (c_int * nsenses)()
            for i in range(nsenses):
                psenses[i] = senses[i]

        # children contains both edges and p-curves
        if oclass == LOOP and geom is not None:
            nchildren = int(nchildren/2)

        if nchildren > 0 and nsenses > 0 and nchildren != nsenses:
            errmsg = 'Number of children is inconsistent with senses'
            raise ValueError(errmsg)

        obj = c_ego()
        stat = _egads.EG_makeTopology(self._context, geo, oclass, mtype,
                                      data, c_int(nchildren), pchildren,
                                      psenses, ctypes.byref(obj))
        if stat: _raiseStatus(stat)

        # BODY make copies all children, so no need to reference the children
        if oclass == BODY: refs=None

        # MODEL take ownership of BODYs and the BODYs can nolonger be deleted outside the model
        if oclass == MODEL:
            for ref in children:
                if ref._finalize is not None:
                    ref._finalize.detach()
                    ref._finalize = None

        return ego(obj, self, deleteObject=True, refs=refs)

#=============================================================================-
    def makeSolidBody(self, stype, data):
        """
        Creates a simple SOLIDBODY. Can be either a box, cylinder,
        sphere, cone, or torus.

        Parameters
        ----------
        stype: BOX, SPHERE, CONE, CYLINDER or TORUS

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
        rvec = (ctypes.c_double * len(data))()
        for i in range(len(data)): rvec[i] = data[i]

        body = c_ego()
        stat = _egads.EG_makeSolidBody(self._context, stype, rvec, ctypes.byref(body))
        if stat: _raiseStatus(stat)

        return ego(body, self, deleteObject=True)

#=============================================================================-
    def approximate(self, sizes, xyzs, mDeg=0, tol=1.e-8):
        """
        Computes and returns the resultant geometry object created by
        approximating the data by a BSpline (OCC or EGADS method).

        Parameters
        ----------
        sizes:
            a vector of 2 integers that specifies the size and dimensionality of
            the data. If the second is zero, then a CURVE is fit and the first
            integer is the length of the number of [x,y,z] triads. If the second
            integer is nonzero then the input data reflects a 2D map.

        xyzs:
            the data to fit (list of [x,y,z] triads)

        mDeg:
            the maximum degree used by OCC [3-8], or cubic by EGADS [0-2]
            0 - fixes the bounds and uses natural end conditions
            1 - fixes the bounds and maintains the slope input at the bounds
            2 - fixes the bounds & quadratically maintains the slope at 2 nd order

        tol:
            is the tolerance to use for the BSpline approximation procedure,
            zero for a SURFACE fit (OCC).

        Returns
        -------
        the approximated (or fit) BSpline resultant ego
        """
        psizes = (c_int * 2)()
        psizes[0] = sizes[0]
        psizes[1] = sizes[1]

        length = sizes[0] if sizes[1] <= 0 else sizes[0] * sizes[1]

        # Allocate the points
        pxyz = (c_double * (3*length))()
        for i in range(length):
            pxyz[3*i  ] = xyzs[i][0]
            pxyz[3*i+1] = xyzs[i][1]
            pxyz[3*i+2] = xyzs[i][2]

        geom = c_ego()
        stat = _egads.EG_approximate(self._context, c_int(mDeg), c_double(tol), psizes,
                                     pxyz, ctypes.byref(geom))
        if stat: _raiseStatus(stat)

        return ego(geom, self, deleteObject=True)

#=============================================================================-
    def fitTriangles(self, xyzs, tris, tric=None, tol=1e-8):
        """
        Computes and returns the resultant geometry object created by approximating
        the triangulation by a BSpline surface.

        Parameters
        ----------
        xyzs:
            the coordinates to fit (list of [x,y,z] triads)

        tris:
            triangle indices (1 bias) (list of 3-tuples)

        tric:
            neighbor triangle indices (1 bias) -- 0 or (-) at bounds
            None -- will compute (len(tris) in length, if not None)

        tol:
            the tolerance to use for the BSpline approximation procedure

        Returns
        -------
        The approximated (or fit) BSpline resultant ego
        """
        npnt = len(xyzs)
        pxyz = (c_double * (3*npnt))()
        ntri = len(tris)
        ptris = (c_int * (3*ntri))()

        for i in range(npnt):
            pxyz[3*i:3*i+3] = xyzs[i][:]

        for i in range(ntri):
            ptris[3*i:3*i+3] = tris[i][:]

        ptric = None
        if tric is not None:
            ptric = (c_int * (3*ntri))()
            for i in range(ntri):
                ptric[3*i:3*i+3] = tric[i][:]

        geom = c_ego()
        stat = _egads.EG_fitTriangles(self._context, c_int(npnt), pxyz,
                                      c_int(ntri), ptris, ptric,
                                      c_double(tol), ctypes.byref(geom))
        if stat: _raiseStatus(stat)

        return ego(geom, self, deleteObject=True)

#=============================================================================-
    def loadModel(self, name, bitFlag=0):
        """
        Loads a MODEL object from a file

        Parameters
        ----------
        name:
            path of file to load (with extension - case insensitive):
            igs/iges    IGES file
            stp/step    STEP file
            brep        native OpenCASCADE file
            egads       native file format with persistent Attributes, splits ignored)

        bitFlag:
            Options (additive):
                1 Don't split closed and periodic entities
                2 Split to maintain at least C1 in BSPLINEs
                4 Don't try maintaining Units on STEP read (always millimeters) 8 Try to merge Edges and Faces (with same geometry)
                16 Load unattached Edges as WireBodies (stp/step & igs/iges)

        Returns
        -------
        the MODEL ego
        """
        filename = name.encode() if isinstance(name,str) else name

        model = c_ego()
        stat = _egads.EG_loadModel(self._context, c_int(bitFlag), filename, ctypes.byref(model))
        if stat: _raiseStatus(stat)

        return ego(model, self, deleteObject=True)

#=============================================================================-

def _del_ego(obj):
    stat = _egads.EG_deleteObject(obj)
    if stat:
        _raiseStatus(stat)

class ego:
    """
    A python wrapper class for an EGADS c_ego object.
    """
#=============================================================================-
    def __init__(self, obj, context=None, deleteObject=False, refs=None):
        """
        Constructor for internal use. An ego class should not be constructed directly.

        Parameters
        ----------
        obj:
            an egads c_ego (any object except CONTEXT)

        context:
            the Context associated with the ego

        deleteObject:
            If True then this ego takes ownershop of obj and calls EG_deleteObject during garbage collection

        refs:
            list of ego's that are used internaly in obj. Used to ensure proper deletion order
        """
        self._obj = obj
        self._refs = refs # any children referenced in this instance

        self._finalize = None
        if deleteObject:
           self._finalize = finalize(self, _del_ego, self._obj)

        if context is None:
            c_context = c_ego()
            stat = _egads.EG_getContext(obj, ctypes.byref(c_context))
            if stat: _raiseStatus(stat)

            self.context = Context(c_context)
        else:
            self.context = context

        if not isinstance(obj, c_ego):
            raise EGADSError("obj must be a c_ego instance")

#=============================================================================-
    def py_to_c(self, takeOwnership=False):
        """
        Returns the c_ego pointer
        
        Parameters
        ----------
        takeOwnership:
            if True, pyEGADS will no longer delete the returned c_ego during garbage collection
        """
        if takeOwnership:
            self._finalize.detach()
            self._finalize = None
        return self._obj

#=============================================================================-
    def __eq__(self, obj):
        """
        checks pointer equality of underlying c_ego memory address
        """
        return isinstance(obj, ego) and ctypes.addressof(obj._obj.contents) == ctypes.addressof(self._obj.contents)

#=============================================================================-
    def __ne__(self, obj):
        return not self == obj

#=============================================================================-
    def __bool__(self):
        return self._obj.__bool__()

#=============================================================================-
    def __nonzero__(self):
        return self._obj.__nonzero__()

#=============================================================================-
    def __del__(self):
        #if hasattr(self, "name"):
        #    print("delete", self.name, self._finalize, self._refs, self.context)
        #else:
        #    print("delete no name", self._finalize, self._refs, self.context)
        if self._finalize is not None:
            self._finalize()

        # refernces can now be deleted
        self._obj  = None
        self._refs = None

#=============================================================================-
    def isSame(self, obj):
        """
        Compares two objects for geometric equivalence.

        Parameters
        ----------
        obj:
            an object to compare with
        """
        stat = _egads.EG_isSame(self._obj, obj._obj)
        return stat == EGADS_SUCCESS

#=============================================================================-
    def isEquivalent(self, obj):
        """
        Compares two topological objects for equivalence.

        Parameters
        ----------
        obj:
            an object to compare with
        """
        stat = _egads.EG_isEquivalent(self._obj, obj._obj)
        return stat == EGADS_SUCCESS

#=============================================================================-
    def getInfo(self):
        """
        Return information about the object

        Returns
        -------
        oclass:
            object class type

        mtype:
            object sub-type

        topRef:
            is the top level BODY/MODEL that owns the object or context (if top)

        prev:
            is the previous object in the threaded list (None at CONTEXT)

        next:
            is the next object in the list (None is the end of the list)
        """
        oclass = c_int()
        mtype  = c_int()
        topObj = c_ego()
        prev   = c_ego()
        next   = c_ego()
        stat = _egads.EG_getInfo(self._obj, ctypes.byref(oclass), ctypes.byref(mtype),
                          ctypes.byref(topObj), ctypes.byref(prev), ctypes.byref(next))
        if stat: _raiseStatus(stat)
        return oclass.value, mtype.value, ego(topObj, self.context), ego(prev, self.context), ego(next, self.context)

#=============================================================================-
    def saveModel(self, name, overwrite=False):
        """
        Writes the geometry to disk. Only writes analytic geometric data for anything but EGADS output.
        Will not overwrite an existing file of the same name unless overwrite = True.

        Parameters
        ----------
        self:
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
        filename = name.encode() if isinstance(name,str) else str
        if overwrite and os.path.exists(filename):
            os.remove(filename)
        stat = _egads.EG_saveModel(self._obj, filename)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def exportModel(self):
        """
        Create a stream of data serializing the objects in the Model.

        Returns
        -------
        the byte-stream
        """
        nbyte = c_size_t()
        pstream = c_char_p()
        stat = _egads.EG_exportModel(self._obj, ctypes.byref(nbyte), ctypes.byref(pstream))
        if stat: _raiseStatus(stat)

        stream = pstream.value[:nbyte.value]
        _egads.EG_free(pstream)

        return stream

#=============================================================================-
    def getTransform(self):
        """
        Returns the [3][4] transformation information.
        """
        T = (c_double * 12)()
        stat = _egads.EG_getTransformation(self._obj, T)
        if stat: _raiseStatus(stat)
        return tuple((T[4*i], T[4*i+1], T[4*i+2], T[4*i+3]) for i in range(3))

#=============================================================================-
    def copyObject(self, other = None):
        """
        Creates a new EGADS Object by copying and optionally transforming the input object.
        If the Object is a Tessellation Object, then other can be a vector of
        displacements that is 3 times the number of vertices of doubles in length.

        Parameters
        ----------
        other:
            the transformation or context object (an ego) --
            None for a strict copy
            a displacement vector for TESSELLATION Objects only
            (number of global indices by 3 doubles in length)

        Returns
        -------
        the resultant new ego
        """

        oform = None
        if other is not None:
            if isinstance(other, ego):
                oform = other._obj
            else:
                oform = (c_double*3)()
                oform[:] = other[:3]

        obj = c_ego()
        stat = _egads.EG_copyObject(self._obj, oform, ctypes.byref(obj))
        if stat: _raiseStatus(stat)
        return ego(obj, self.context, deleteObject=True)

    def __copy__(self):
        return self.copyObject()

    def __deepcopy__(self, memo):
        obj = self.__copy__()
        memo[id(self)] = obj
        return obj

#=============================================================================-
    def flipObject(self):
        """
        Creates a new EGADS Object by copying and reversing the input object.
        Can be Geometry (flip the parameterization) or Topology (reverse the sense).
        Not for Node, Body or Model. Surfaces reverse only the u parameter.

        Returns
        -------
        the resultant new flipped ego
        """
        obj = c_ego()
        stat = _egads.EG_flipObject(self._obj, ctypes.byref(obj))
        if stat: _raiseStatus(stat)

        return ego(obj, self.context, deleteObject=True)

#=============================================================================-
    def getGeometry(self):
        """
        Returns information about the geometric object.
        Notes: ints is returned for either mtype = BEZIER or BSPLINE.

        Returns
        -------
        oclass:
            the returned Object Class: PCURVE, CURVE or SURFACE

        mtype:
            the returned Member Type (depends on oclass)

        reals:
            the returned pointer to real data used to describe the geometry

        ints:
            the returned pointer to integer information (None if none)

        rGeom:
            the returned reference Geometry Object (None if none)
        """
        oclass = c_int()
        mtype = c_int()
        ivec = POINTER(c_int)()
        rvec = POINTER(c_double)()
        c_geom = c_ego()
        stat = _egads.EG_getGeometry(self._obj, ctypes.byref(oclass), ctypes.byref(mtype), ctypes.byref(c_geom),
                                     ctypes.byref(ivec), ctypes.byref(rvec))
        if stat: _raiseStatus(stat)

        nivec = c_int()
        nrvec = c_int()
        _egads.EG_getGeometryLen(self._obj, ctypes.byref(nivec), ctypes.byref(nrvec))

        oclass = oclass.value
        mtype = mtype.value

        geom = None
        if c_geom:
            geom = ego(c_geom)

        ints = None
        if ivec:
            ints = [ivec[i] for i in range(nivec.value)]
            _egads.EG_free(ivec)

        if mtype == BSPLINE:
            if oclass == PCURVE or oclass == CURVE:
                flag   = ints[0]
                nCP    = ints[2]
                nKnots = ints[3]
                reals = [list(rvec[:nKnots])]
                if oclass == PCURVE:
                    reals.append([(rvec[nKnots+2*i], rvec[nKnots+2*i+1]) for i in range(nCP)])
                if oclass == CURVE:
                    reals.append([(rvec[nKnots+3*i], rvec[nKnots+3*i+1], rvec[nKnots+3*i+2]) for i in range(nCP)])
                if flag == 2:
                    reals.append(list(rvec[nrvec-nCP:nrvec]))
            else:
                flag    = ints[0]
                nCPu    = ints[2]
                nUKnots = ints[3]
                nCPv    = ints[5]
                nVKnots = ints[6]
                nKnots = nUKnots + nVKnots
                reals = [list(rvec[       :nUKnots]),
                         list(rvec[nUKnots: nKnots])]
                reals.append([(rvec[nKnots+3*i], rvec[nKnots+3*i+1], rvec[nKnots+3*i+2]) for i in range(nCPu*nCPv)])
                if flag == 2:
                    reals.append(list(rvec[nrvec-nCPu*nCPv:nrvec]))

        elif mtype == BEZIER:
            if oclass == PCURVE:
                nrvec = int(nrvec.value/2)
                reals = [(rvec[2*i], rvec[2*i+1]) for i in range(nrvec)]
            else:
                nrvec = int(nrvec.value/3)
                reals = [(rvec[3*i], rvec[3*i+1], rvec[3*i+2]) for i in range(nrvec)]
        else:
            reals = [rvec[i] for i in range(nrvec.value)]

        _egads.EG_free(rvec)

        return oclass, mtype, reals, ints, geom

#=============================================================================-
    def getRange(self):
        """
        Returns the valid range of the object: may be one of PCURVE,
        CURVE, SURFACE, EDGE or FACE

        Returns
        -------
        range:
            for PCURVE, CURVE or EDGE returns 2 values, t-start and t-end
            for SURFACE or FACE returns 4 values, u-min, u-max, v-min and v-max

        periodic:
            0 for non-periodic, 1 for periodic in t or u 2 for periodic in
            v (or-able)
        """
        r = (c_double * 4)()
        periodic = c_int()
        oclass = c_int()
        stat = _egads.EG_getRange(self._obj, r, ctypes.byref(periodic))
        if stat: _raiseStatus(stat)

        oclass = self._obj.contents.oclass
        if (oclass == PCURVE or
            oclass == CURVE or
            oclass == EDGE or
            oclass == EEDGE):
            return (r[0], r[1]), periodic.value
        elif (oclass == SURFACE or
              oclass == FACE or
              oclass == EFACE):
            return (r[0], r[1], r[2], r[3]), periodic.value

        return None

#=============================================================================-
    def getBoundingBox(self):
        """
        Computes the Cartesian bounding box around the object:

        Returns
        -------
        bounding box [(x_min,y_min,z_min),
                      (x_max,y_max,z_max)]
        """
        data = (c_double * 6)()
        stat = _egads.EG_getBoundingBox(self._obj, data)
        if stat: _raiseStatus(stat)
        return [(data[0], data[1], data[2]), (data[3], data[4], data[5])]

#=============================================================================-
    def evaluate(self, params):
        """
        Returns the result of evaluating the object at a parameter
        point. May be used for PCURVE, CURVE, SURFACE, EDGE or FACE.

        In the return: X = [x,y,z]

        Parameters
        ----------
        params:
            The parametric location
            For NODE: None
            For PCURVE, CURVE, EDGE, EEDGE: The t-location
            For SURFACE, FACE, EFACE: The (u, v) coordiantes

        Returns
        -------
            For NODE: X
            For PCURVE, CURVE, EDGE: X, X_t and X_tt
            For SURFACE, FACE: X, [X_u, X_v], [X_uu, X_uv, X_vv]
        """
        param = (c_double * 2)()
        oclass = self._obj.contents.oclass
        try:
            if (oclass == EDGE or
                oclass == EEDGE or
                oclass == CURVE or
                oclass == PCURVE):
                if isinstance(params, (int, float)):
                    param[0] = params
                else:
                    param[0] = params[0]
            elif (oclass == FACE or
                  oclass == EFACE or
                  oclass == SURFACE):
                param[0] = params[0]
                param[1] = params[1]
        except:
            errmsg = 'Failed to convert parameter value: params = ' + str(params)
            raise ValueError(errmsg)

        r = (c_double * 18)()
        stat = _egads.EG_evaluate(self._obj, param, r)
        if stat: _raiseStatus(stat)

        if oclass == NODE:
            return (r[0], r[1], r[2])
        elif oclass == PCURVE:
            return (r[0], r[1]), (r[2], r[3]), (r[4], r[5])
        elif (oclass == EDGE or
              oclass == EEDGE or
              oclass == CURVE):
            return (r[0], r[1], r[2]), (r[3], r[4], r[5]), (r[6], r[7], r[8])
        elif (oclass == FACE or
              oclass == EFACE or
              oclass == SURFACE):
            return (r[0], r[1], r[2]), \
                   ((r[3], r[4], r[5]), (r[6], r[7], r[8])), \
                   ((r[9], r[10], r[11]),
                    (r[12], r[13], r[14]),
                    (r[15], r[16], r[17]))
        return None

#=============================================================================-
    def invEvaluate(self, pos):
        """
        Returns the result of inverse evaluation on the object.
        For topology the result is limited to inside the EDGE/FACE valid bounds.

        Parameters
        ----------
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
        xp = (c_double * 3)()
        params = (c_double * 2)()
        result = (c_double * 3)()

        oclass = self._obj.contents.oclass
        if oclass == PCURVE:
            xp[0] = pos[0]
            xp[1] = pos[1]
        else:
            xp[0] = pos[0]
            xp[1] = pos[1]
            xp[2] = pos[2]

        stat = _egads.EG_invEvaluate(self._obj, xp, params, result)
        if stat: _raiseStatus(stat)

        if oclass == PCURVE:
            return params[0], [result[0], result[1]]
        elif (oclass == EDGE or
              oclass == CURVE):
            return params[0], [result[0], result[1], result[2]]
        elif (oclass == FACE or
              oclass == SURFACE):
            return [params[0], params[1]], [result[0], result[1], result[2]]
        return None, None

#=============================================================================-
    def getTolerance(self):
        """
        Returns the internal tolerance defined for the object
        """
        tol = c_double()
        stat = _egads.EG_getTolerance(self._obj, ctypes.byref(tol))
        if stat: _raiseStatus(stat)
        return tol.value

#=============================================================================-
    def tolerance(self):
        """
        Returns the maximum tolerance defined for the object's hierarchy
        """
        tol = c_double()
        stat = _egads.EG_tolerance(self._obj, ctypes.byref(tol))
        if stat: _raiseStatus(stat)
        return tol.value

#=============================================================================-
#    depricated
#    def setTolerance(self, tol):
#        """
#        Sets the internal tolerance defined for the object. Useful for SBOs.
#        """
#        stat = _egads.EG_setTolerance(self._obj, c_double(tol))
#        if stat: _raiseStatus(stat)

#=============================================================================-
    def getBody(self):
        """
        Get the body that this object is contained within

        Returns
        -------
        The BODY ego
        """
        body = c_ego()
        stat = _egads.EG_getBody(self._obj, ctypes.byref(body))
        if stat: _raiseStatus(stat)
        return ego(body, self.context) if body else None

#=============================================================================-
    def arcLength(self, t1, t2):
        """
        Get the arc-length of the  dobject

        Parameters
        ----------
        t1:
            The starting t value

        t2:
            The end t value

        Returns
        -------
        the arc-length
        """
        length = c_double()
        stat = _egads.EG_arcLength(self._obj, c_double(t1), c_double(t2), ctypes.byref(length))
        if stat: _raiseStatus(stat)
        return length.value

#=============================================================================-
    def getArea(self, limits=None):
        """
        Computes the surface area from a LOOP, a surface or a FACE.
        When a LOOP is used a planar surface is fit
        and the resultant area can be negative if the orientation of
        the fit is opposite of the LOOP.

        Parameters
        ----------
        limits:
            may be None except must contain the limits for a surface (4 words)

        Returns
        -------
        the surface area
        """
        lim = None
        area = c_double()
        if self._obj.contents.oclass == SURFACE:
            lim = (c_double * 4)()
            lim[0] = limits[0]
            lim[1] = limits[1]
            lim[2] = limits[2]
            lim[3] = limits[3]
        stat = _egads.EG_getArea(self._obj, lim, ctypes.byref(area))
        if stat: _raiseStatus(stat)
        return area.value

#=============================================================================-
    def curvature(self, params):
        """
        Returns the curvature and principle directions/tangents

        Parameters
        ----------
        params:
            The parametric location
            For PCURVE, CURVE, EDGE, EEDGE: The t-location
            For SURFACE, FACE, EFACE: The (u, v) coordiantes

        Returns
        -------
        For PCURVE: curvature, [dir.x, dir.y]
        For CURVE, EDGE: curvature, [dir.x ,dir.y, dir.z]
        For SURFACE, FACE: curvature1, [dir1.x, dir1.y, dir1.z], curvature2, [dir2.x, dir2.y, dir2.z]
        """
        param = (c_double * 2)()
        r = (c_double * 8)()
        oclass = self._obj.contents.oclass
        try:
            if (oclass == EDGE or
                oclass == EEDGE or
                oclass == CURVE or
                oclass == PCURVE):
                if isinstance(params, (int, float)):
                    param[0] = params
                else:
                    param[0] = params[0]
            elif (oclass == FACE or
                  oclass == EFACE or
                  oclass == SURFACE):
                param[0] = params[0]
                param[1] = params[1]
        except:
            errmsg = 'Failed to convert parameter value: p = ' + str(p)
            raise ValueError(errmsg)

        stat = _egads.EG_curvature(self._obj, param, r)
        if stat: _raiseStatus(stat)

        if oclass == PCURVE:
            return r[0], (r[1], r[2])
        elif (oclass == EDGE or
              oclass == EEDGE or
              oclass == CURVE):
            return r[0], (r[1], r[2], r[3])
        elif (oclass == FACE or
              oclass == EFACE or
              oclass == SURFACE):
            return r[0], (r[1], r[2], r[3]), \
                   r[4], (r[5], r[6], r[7])
        return None

#=============================================================================-
    def getMassProperties(self):
        """
        Computes and returns the physical and inertial properties of a topological object
        that can be EDGE, LOOP, FACE, SHELL or BODY

        Returns
        -------
        volume, surface area (length for EDGE, LOOP or WIREBODY) center of gravity (3)
        inertia matrix at CoG (9)
        """
        lim = None
        data = (c_double*14)()

        stat = _egads.EG_getMassProperties(self._obj, data)

        volume = data[0]
        areaOrLen = data[1]
        CG = tuple(data[2:5])
        I = tuple(data[5:14])

        return volume, areaOrLen, CG, I

#=============================================================================-
    def attributeAdd(self, aname, data):
        """
        Adds an attribute to the object. If an attribute exists with
        the name it is overwritten with the new information.

        Parameters
        ----------
        name:
            The name of the attribute. Must not contain a space or other
            special characters

        data:
            attribute data to add to the ego
            must be all intergers, all floats, a string, or a csystem
        """
        length = c_int(0)
        ints = None
        reals = None
        chars = None
        name = aname.encode() if isinstance(aname,str) else aname

        if all(isinstance(x, float) for x in data) or isinstance(data, csystem):
            atype = ATTRCSYS if isinstance(data, csystem) else ATTRREAL
            length = len(data)
            reals = (c_double * length)()
            for i in range(length):
                reals[i] = float(data[i])
        elif all(isinstance(x, int) for x in data):
            atype = ATTRINT
            length = len(data)
            ints = (c_int * length)()
            for i in range(length):
                ints[i] = data[i]
        elif isinstance(data, str):
            atype = ATTRSTRING
            chars = data.encode() if isinstance(data,str) else data
        else:
            errmsg = "Data list is not all int, double, or str\n"
            errmsg = "data = " + str(data)
            raise EGADSError(errmsg)

        stat = _egads.EG_attributeAdd(self._obj, name, atype, length,
                                      ints, reals, chars)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def attributeDel(self, aname=None):
        """
        Deletes an attribute from the object.

        Parameters
        ----------
        name:
            Name of the attribute to delete
            If the aname is None then
            all attributes are removed from this object
        """
        name = None
        if aname is not None:
            name = aname.encode() if isinstance(aname,str) else aname

        stat = _egads.EG_attributeDel(self._obj, name)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def attributeNum(self):
        """
        Returns the number of attributes found with this object
        """
        num = c_int()
        stat = _egads.EG_attributeNum(self._obj, ctypes.byref(num))
        if stat: _raiseStatus(stat)
        return num.value

#=============================================================================-
    def attributeGet(self, index):
        """
        Retrieves a specific attribute from the object
        """
        name = c_char_p()
        atype = c_int()
        length = c_int()
        ints = POINTER(c_int)()
        reals = POINTER(c_double)()
        chars = ctypes.c_char_p()

        stat = _egads.EG_attributeGet(self._obj, c_int(index), ctypes.byref(name),
                                      ctypes.byref(atype), ctypes.byref(length),
                                      ctypes.byref(ints), ctypes.byref(reals), ctypes.byref(chars))
        if stat: _raiseStatus(stat)

        if atype.value == ATTRINT:
            res = [ints[i] for i in range(length.value)]
        elif atype.value == ATTRREAL:
            res = [reals[i] for i in range(length.value)]
        elif atype.value == ATTRSTRING:
            res = _decode(chars.value)
        elif atype.value == ATTRCSYS:
            res = csystem([reals[i] for i in range(length.value)])
            res.ortho = [[reals[length.value + 3*j+i] for i in range(3)] for j in range(4)]

        return _decode(name.value), res

#=============================================================================-
    def attributeRet(self, aname):
        """
        Retrieves a specific attribute from the object
        """
        atype = c_int()
        length = c_int()
        ints = POINTER(c_int)()
        reals = POINTER(c_double)()
        chars = c_char_p()
        name = aname.encode() if isinstance(aname,str) else aname

        stat = _egads.EG_attributeRet(self._obj, name, ctypes.byref(atype), ctypes.byref(length),
                               ctypes.byref(ints), ctypes.byref(reals), ctypes.byref(chars))
        if stat == EGADS_NOTFOUND:
            return None
        elif stat:
            _raiseStatus(stat)

        if atype.value == ATTRINT:
            return [ints[i] for i in range(length.value)]
        elif atype.value == ATTRREAL:
            return [reals[i] for i in range(length.value)]
        elif atype.value == ATTRCSYS:
            res = csystem([reals[i] for i in range(length.value)])
            res.ortho = [[reals[length.value + 3*j+i] for i in range(3)] for j in range(4)]
            return res
        elif atype.value == ATTRSTRING:
            return _decode(chars.value)

        return None

#=============================================================================-
    def attributeDup(self, dup):
        """
        Removes all attributes from the destination object, then copies the
        attributes from the source
        """
        stat = _egads.EG_attributeDup(self._obj, dup._obj)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def attributeAddSeq(self, aname, data):
        """
        Add additional attributes with the same root name to a sequence.

        Parameters
        ----------
        name:
            The name of the attribute. Must not contain a space or other
            special characters

        data:
            attribute data to add to the ego
            must be all intergers, all floats, a string, or a csystem

        Returns
        -------
        the new number of entries in the sequence
        """
        length = c_int(0)
        ints = None
        reals = None
        chars = None
        name = aname.encode() if isinstance(aname,str) else aname

        if all(isinstance(x, float) for x in data) or isinstance(data, csystem):
            atype = ATTRCSYS if isinstance(data, csystem) else ATTRREAL
            length = len(data)
            reals = (c_double * length)()
            for i in range(length):
                reals[i] = float(data[i])
        elif all(isinstance(x, int) for x in data):
            atype = ATTRINT
            length = len(data)
            ints = (c_int * length)()
            for i in range(length):
                ints[i] = data[i]
        elif isinstance(data, str):
            atype = ATTRSTRING
            chars = data.encode() if isinstance(data,str) else data
        else:
            errmsg = "Data list is not all int, double, or str\n"
            errmsg = "data = " + str(data)
            raise EGADSError(errmsg)

        stat = _egads.EG_attributeAddSeq(self._obj, name, atype, length,
                                         ints, reals, chars)
        if stat < 0: _raiseStatus(stat)

        return stat

#=============================================================================-
    def attributeNumSeq(self, aname):
        """
        Returns the number of named sequenced attributes found on this object.

        Parameters
        ----------
        name:
            The name of the attribute. Must not contain a space or other
            special characters
        """
        name = aname.encode() if isinstance(aname,str) else aname

        num = c_int()
        stat = _egads.EG_attributeNumSeq(self._obj, name, ctypes.byref(num))
        if stat: _raiseStatus(stat)
        return num.value

#=============================================================================-
    def attributeRetSeq(self, aname, index):
        """
        Retrieves a specific attribute from the object
        """
        atype = c_int()
        length = c_int()
        ints = POINTER(c_int)()
        reals = POINTER(c_double)()
        chars = c_char_p()
        name = aname.encode() if isinstance(aname,str) else aname

        stat = _egads.EG_attributeRetSeq(self._obj, name, c_int(index), ctypes.byref(atype), ctypes.byref(length),
                               ctypes.byref(ints), ctypes.byref(reals), ctypes.byref(chars))
        if stat == EGADS_NOTFOUND:
            return None
        elif stat:
            _raiseStatus(stat)

        if atype.value == ATTRINT:
            return [ints[i] for i in range(length.value)]
        elif atype.value == ATTRREAL:
            return [reals[i] for i in range(length.value)]
        elif atype.value == ATTRCSYS:
            res = csystem([reals[i] for i in range(length.value)])
            res.ortho = [[reals[length.value + 3*j+i] for i in range(3)] for j in range(4)]
            return res
        elif atype.value == ATTRSTRING:
            return _decode(chars.value)

        return None

#=============================================================================-
    def otherCurve(self, iCrv, tol=0.0):
        """
        Computes and returns the other curve that matches the input curve.
        If the input curve is a PCURVE, the output is a 3D CURVE (and vice versa).

        Parameters
        ----------
        iCrv:
            the input PCURVE or CURVE/EDGE object

        tol:
            is the tolerance to use when fitting the output curve

        Returns
        -------
        the returned approximated resultant curve object
        """
        oCrv = c_ego()
        stat = _egads.EG_otherCurve(self._obj, iCrv._obj, c_double(tol), ctypes.byref(oCrv))
        if stat: _raiseStatus(stat)

        return ego(oCrv, self.context, deleteObject=True)

#=============================================================================-
    def isoCline(self, iUV, value):
        """
        Computes from the input Surface and returns the isocline curve.

        Parameters
        ----------
        iUV:
            the type of isocline: UISO (0) constant U or VISO (1) constant V

        value:
            the value used for the isocline

        Returns
        -------
        the returned resultant curve object
        """
        oCrv = c_ego()
        stat = _egads.EG_isoCline(self._obj, c_int(iUV), c_double(value), ctypes.byref(oCrv))
        if stat: _raiseStatus(stat)

        return ego(oCrv, self.context, deleteObject=True)

#=============================================================================-
    def convertToBSpline(self):
        """
        Computes and returns the BSpline representation of the geometric object
        which can be a PCURVE, CURVE, EDGE, SURFACE or FACE

        Returns
        -------
        the returned approximated resultant BSPLINE object
        can be the input ego when a BSPLINE and the same limits
        """
        bspline = c_ego()
        stat = _egads.EG_convertToBSpline(self._obj, ctypes.byref(bspline))
        if stat: _raiseStatus(stat)

        return ego(bspline, self.context, deleteObject=True)

#=============================================================================-
    def convertToBSplineRange(self, range):
        """
        Computes and returns the BSpline representation of the geometric object
        which can be a PCURVE, CURVE, or SURFACE

        Required when converting Geometry Objects with infinite range.

        Parameters
        ----------
        range:
            t range (2) or [u, v] box (4) to limit the conversion

        Returns
        -------
        the returned approximated resultant BSPLINE object
        can be the input ego when a BSPLINE and the same limits
        """
        bspline = c_ego()
        alen = len(range)
        prange = (c_double*alen)()
        prange[:] = range[:]

        stat = _egads.EG_convertToBSplineRange(self._obj, prange, ctypes.byref(bspline))
        if stat: _raiseStatus(stat)

        return ego(bspline, self.context, deleteObject=True)

#=============================================================================-
    def getTopology(self):
        """
        Returns information about the topological object

        Returns
        -------
        oclass:
            is NODE, EGDE, LOOP, FACE, SHELL, BODY or MODEL

        mtype:
            for EDGE is TWONODE, ONENODE or DEGENERATE
            for LOOP is OPEN or CLOSED
            for FACE is either SFORWARD or SREVERSE
            for SHELL is OPEN or CLOSED
            BODY is either WIREBODY, FACEBODY, SHEETBODY or SOLIDBODY

        geom:
            The reference geometry object (if none this is returned as None)

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
        geo = c_ego()
        oclass = c_int()
        mtype = c_int()
        limits = (c_double * 4)()
        nchildren = c_int()
        pchildren = POINTER(c_ego)()
        psenses = POINTER(c_int)()
        stat = _egads.EG_getTopology(self._obj, ctypes.byref(geo), ctypes.byref(oclass), ctypes.byref(mtype),
                                     limits, ctypes.byref(nchildren), ctypes.byref(pchildren), ctypes.byref(psenses))
        if stat: _raiseStatus(stat)

        oclass = oclass.value
        mtype = mtype.value
        nchildren = nchildren.value

        nchild = nchildren
        if oclass == LOOP  and geo      : nchild *= 2
        if oclass == MODEL and mtype > 0: nchild = mtype
        children = [ego(c_ego(pchildren[i].contents), self.context) for i in range(nchild)]

        senses = [psenses[i] for i in range(nchildren)] if psenses else None

        reals = []
        if oclass == NODE:
            reals = [limits[0], limits[1], limits[2]]
        elif oclass == EDGE:
            reals = [limits[0], limits[1]]
        elif oclass == FACE or oclass == EFACE:
            reals = [limits[0], limits[1], limits[2], limits[3]]

        geom = ego(c_ego(geo.contents), self.context) if geo else None

        return oclass, mtype, geom, reals, children, senses

#=============================================================================-
    def getBodyTopos(self, oclass, ref=None):
        """
        Returns topologically connected objects:

        Parameters
        ----------
        oclass:
            is NODE, EGDE, LOOP, FACE or SHELL -- must not be the same
            class as ref

        ref:
            reference topological object or NULL. Sets the context for the
            returned objects (i.e. all objects of a class [oclass] in the
            tree looking towards that class from ref). None starts from
            the BODY (for example all NODEs in the BODY)

        Returns
        -------
        the returned number of requested topological objects is a
        returned pointer to the block of objects
        """

        src = c_ego()
        ntopos = c_int(0)
        topos = POINTER(c_ego)()
        if ref is not None:
            src = ref._obj

        stat = _egads.EG_getBodyTopos(self._obj, src, oclass, ctypes.byref(ntopos), ctypes.byref(topos))
        if stat: _raiseStatus(stat)

        tlist = [ego(c_ego(topos[i].contents), self.context) for i in range(ntopos.value)]
        _egads.EG_free(topos)

        return tlist

#=============================================================================-
    def matchBodyEdges(self, body, toler=0.0):
        """
        Examines the EDGEs in one BODY against all of the EDGEs in another.
        If the number of NODEs, the NODE locations, the EDGE bounding boxes
        and the EDGE arc lengths match it is assumed that the EDGEs match.
        A list of pairs of indices are returned.

        Parameters
        ----------
        body:
            body container object

        toler:
            the tolerance used (can be zero to use entity tolerances)

        Returns
        -------
        a list of the tuples of matching indices
        """
        nmatches = c_int(0)
        match = POINTER(c_int)()

        if not (body._obj.contents.oclass == BODY and self._obj.contents.oclass == BODY):
            errmsg = 'Both objects must have object class BODY'
            raise ValueError(errmsg)

        stat = _egads.EG_matchBodyEdges(self._obj, body._obj, c_double(toler),
                                ctypes.byref(nmatches), ctypes.byref(match))
        if stat: _raiseStatus(stat)

        matches = [(match[2*i], match[2*i+1]) for i in range(nmatches.value)]

        _egads.EG_free(match)

        return matches

#=============================================================================-
    def matchBodyFaces(self, body, toler=0.0):
        """
        Examines the FACEs in one BODY against all of the FACEs in
        another. If the number of LOOPs, number of NODEs, the NODE
        locations, the number of EDGEs and the EDGE bounding boxes as
        well as the EDGE arc lengths match it is assumed that the
        FACEs match. A list of pairs of indices are returned.

        Parameters
        ----------
        body:
            body container object

        toler:
            the tolerance used (can be zero to use entity tolerances)

        Returns
        -------
        a list of the tuples of matching indices

        Note: This is useful for the situation where there are
              glancing FACEs and a UNION operation fails (or would
              fail). Simply find the matching FACEs and do not include them
              in a call to EG_sewFaces.
        """
        nmatches = c_int(0)
        match = POINTER(c_int)()

        if not (body._obj.contents.oclass == BODY and self._obj.contents.oclass == BODY):
            errmsg = 'Both objects must have object class BODY'
            raise ValueError(errmsg)

        stat = _egads.EG_matchBodyFaces(self._obj, body._obj, c_double(toler),
                                ctypes.byref(nmatches), ctypes.byref(match))
        if stat: _raiseStatus(stat)

        matches = [(match[2*i], match[2*i+1]) for i in range(nmatches.value)]

        _egads.EG_free(match)

        return matches

#=============================================================================-
    def indexBodyTopo(self, obj):
        """
        Return the (1-based) index of the topological object in the BODY/EBODY

        Parameters
        ----------
        obj:
            the Topology Object in the Body/EBody

        Retruns
        -------
        the index
        """
        stat = _egads.EG_indexBodyTopo(self._obj, obj._obj)
        if stat < EGADS_SUCCESS: _raiseStatus(stat)

        return stat

#=============================================================================-
    def objectBodyTopo(self, oclass, index):
        """
        Returns the topological object (based on index) in the BODY/EBODY

        Parameters
        ----------
        oclass:
            is NODE, EDGE, LOOP, FACE or SHELL

        index:
            is the index (bias 1) of the entity requested

        Returns
        -------
        the topological object in the Body
        """
        obj = c_ego()
        stat = _egads.EG_objectBodyTopo(self._obj, c_int(oclass), c_int(index), ctypes.byref(obj))
        if stat < EGADS_SUCCESS: _raiseStatus(stat)

        return ego(obj, self.context)

#=============================================================================-
    def inTopology(self, xyz):
        """
        Computes whether the point is on or contained within the object.
        Works with EDGEs and FACEs by projection. SHELLs must be CLOSED.

        The object, can be EDGE, FACE, SHELL or SOLIDBODY.

        Parameters
        ----------
        xyz:
            the coordinate location to check

        Returns
        -------
        True if xyz is in the topology, False otherwise
        """
        pxyz = (c_double*3)()
        pxyz[:] = xyz[:3]
        stat = _egads.EG_inTopology(self._obj, pxyz)
        if stat < EGADS_SUCCESS: _raiseStatus(stat)

        return stat == EGADS_SUCCESS

#=============================================================================-
    def inFace(self, uv):
        """
        Computes the result of the [u,v] location in the valid part of the FACE.

        Parameters
        ----------
        uv:
            the parametric location to check

        Returns
        -------
        True if uv is in the FACE, False otherwise
        """
        puv = (c_double*2)()
        puv[:] = uv[:2]
        stat = _egads.EG_inFace(self._obj, puv)
        if stat < EGADS_SUCCESS: _raiseStatus(stat)

        return stat == EGADS_SUCCESS

#     This has been superseeded by generalBoolean
#     def solidBoolean(self, tool, oper):
#         """
#         Performs the Solid Boolean Operations (SBOs) on the source
#         BODY Object (that has the type SOLIDBODY). The tool object
#         types depend on the operation. This supports Intersection,
#         Subtraction and Union. The object must be a SOLIDBODY or
#         MODEL.
#         Note: This may be called with src being a MODEL. In this case
#         tool may be a SOLIDBODY for Intersection or a FACE/FACEBODY
#         for Fusion. The input MODEL may contain anything, but must not
#         have duplicate topology.
#
#         Parameters
#         ----------
#         tool: the tool object:
#         either a SOLIDBODY for all operators or a FACE/FACEBODY for
#         Subtraction.
#         oper: the operation to perform
#
#         Returns
#         -------
#         the resultant MODEL object (this is because there may be
#         multiple bodies from either the subtraction or intersection
#         operation).
#         """
#         new_obj = c_ego()
#         stat = _egads.EG_solidBoolean(self._obj, tool._obj, oper, ctypes.byref(new_obj))
#         if stat: _raiseStatus(stat)
#         return ego(new_obj, self.context)

#=============================================================================-
    def generalBoolean(self, tool, oper, tol=0.0):
         """
         Performs the Boolean Operations on the source BODY Object(s).
         The tool Object(s) are applied depending on the operation.
         This supports Intersection, Splitter, Subtraction and Union.

         Parameters
         ----------
         tool:
             the tool BODY object(s) in the form of a MODEL or BODY

         oper:
             1-SUBTRACTION, 2-INTERSECTION, 3-FUSION, 4-SPLITTER

         tol:
             the tolerance applied to perform the operation.
             If set to zero, the minimum tolerance is applied

         Returns
         -------
         the resultant MODEL object (this is because there may be
         multiple bodies from either the subtraction or intersection operation).

         Note: The BODY Object(s) for src and tool can be of any type but
               the results depend on the operation.
               Only works for OpenCASCADE 7.3.1 and higher.
         """
         model = c_ego()
         stat = _egads.EG_generalBoolean(self._obj, tool._obj, c_int(oper), c_double(tol), ctypes.byref(model))
         if stat: _raiseStatus(stat)

         return ego(model, self.context, deleteObject=True)

#=============================================================================-
    def fuseSheets(self, tool):
         """
         Fuses (unions) two SHEETBODYs resulting in a single SHEETBODY.

         Parameters
         ----------
         tool:
            the tool SHEETBODY object

         Returns
         -------
         the resultant SHEETBODY object
         """
         sheet = c_ego()
         stat = _egads.EG_fuseSheets(self._obj, tool._obj, ctypes.byref(sheet))
         if stat: _raiseStatus(stat)

         return ego(sheet, self.context, deleteObject=True)

#=============================================================================-
    def intersection(self, tool):
        """
        Intersects the source BODY Object (that has the type
        SOLIDBODY, SHEETBODY or FACEBODY) with a surface or
        surfaces. The tool object contains the intersecting geometry
        in the form of a FACEBODY, SHEETBODY, SOLIDBODY or a single
        FACE.

        Parameters
        ----------
        tool:
        the FACE/FACEBODY/SHEETBODY/SOLIDBODY tool object

        Returns
        -------
        pairs:
            List of the FACE/EDGE object pairs

        model:
            The resultant MODEL object which contains the set of WIREBODY
            BODY objects (this is because there may be multiple LOOPS as a
            result of the operation).  The EDGE objects contained within
            the LOOPS have the attributes of the FACE in src responsible
            for that EDGE.

        Note: The EDGE objects contained within the LOOPS have the attributes
              of the FACE in src responsible for that EDGE.
        """
        nobj = c_int()
        faceEdgePairs = POINTER(c_ego)()
        obj = c_ego()
        stat = _egads.EG_intersection(self._obj, tool._obj, ctypes.byref(nobj), ctypes.byref(faceEdgePairs),
                                      ctypes.byref(obj))
        if nobj.value == 0:
            return [], None
        elif stat: _raiseStatus(stat)

        pairs = [(ego(c_ego(faceEdgePairs[2*i  ].contents), self.context), \
                  ego(c_ego(faceEdgePairs[2*i+1].contents), self.context)) for i in range(nobj.value)]

        _egads.EG_free(faceEdgePairs)

        return pairs, ego(obj, self.context, deleteObject=True)

#=============================================================================-
    def imprintBody(self, pairs):
        """
        Imprints EDGE/LOOPs on the source BODY Object (that has the type
        SOLIDBODY, SHEETBODY or FACEBODY). The EDGE/LOOPs are
        paired with the FACEs in the source that will be scribed with the
        EDGE/LOOP.

        Parameter
        ---------
        pairs:
            list of FACE/EDGE and/or FACE/LOOP object pairs to scribe
            2*nObj in len -- can be the output from intersect
            result

        Returns
        -------
        the resultant BODY object (with the same type as the input source
        object, though the splitting of FACEBODY objects results in a
        SHEETBODY)
        """
        nobj = len(pairs)
        objs = (c_ego * (nobj*2))()
        for i in range(nobj):
            objs[2*i  ] = pairs[i][0]._obj
            objs[2*i+1] = pairs[i][1]._obj

        body = c_ego()
        stat = _egads.EG_imprintBody(self._obj, c_int(nobj), objs, ctypes.byref(body))
        if stat: _raiseStatus(stat)

        return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def makeFace(self, mtype, rdata=None):
        """
        Creates a simple FACE from a LOOP or a SURFACE. Also can be
        used to hollow a single LOOPed existing FACE. This function
        creates any required NODEs, EDGEs and LOOPs.

        Parameters
        ----------
        self:
            Either a LOOP (for a planar cap), a SURFACE with [u,v] bounds,
            or a FACE to be hollowed out

        mtype:
            Is either SFORWARD or SREVERSE
            For LOOPs you may want to look at the orientation using
            getArea, ignored when the input object is a FACE

        rdata:
            May be None for LOOPs, but must be the limits for a SURFACE (4
            values), the hollow/offset distance and fillet radius (zero is
            for no fillets) for a FACE input object (2 values)

        Returns
        -------
        the resultant returned topological FACE object (a return of
        EGADS_OUTSIDE is the indication that offset distance was too
        large to produce any cutouts, and this result is the input
        object)
        """
        data = None
        oclass = self._obj.contents.oclass
        if oclass == SURFACE:
            # Limits of the surface
            data = (c_double * 4)()
            data[0] = rdata[0]
            data[1] = rdata[1]
            data[2] = rdata[2]
            data[3] = rdata[3]
        elif oclass == FACE:
            data = (c_double * 2)()
            data[0] = rdata[0]
            data[1] = rdata[1]

        face = c_ego()
        stat = _egads.EG_makeFace(self._obj, mtype, data, ctypes.byref(face))
        if stat and stat != EGADS_OUTSIDE: _raiseStatus(stat)

        return ego(face, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def replaceFaces(self, faces):
         """
         Creates a new SHEETBODY or SOLIDBODY from an input SHEETBODY or SOLIDBODY
         and a list of FACEs to modify. The FACEs are input in pairs where
         the first must be an Object in the BODY and the second either
         a new FACE or None. The None replacement flags removal of the FACE in the BODY.

         Parameters
         ----------
         faces:
             list of FACE tuple pairs, where the first must be a FACE in the BODY
             and second is either the FACE to use as a replacement or a None which
             indicates that the FACE is to be removed from the BODY

         Returns
         -------
         the resultant BODY object, either a SHEETBODY or a SOLIDBODY
         (where the input was a SOLIDBODY and all FACEs are replaced
         in a way that the LOOPs match up)
         """
         nfaces = len(faces)
         pfaces = (c_ego *(2*nfaces))()
         for i in range(nfaces):
             pfaces[2*i  ] = faces[i][0]._obj
             pfaces[2*i+1] = None if faces[i][1] is None else faces[i][1]._obj

         body = c_ego()
         stat = _egads.EG_replaceFaces(self._obj, c_int(nfaces), pfaces, ctypes.byref(body))
         if stat: _raiseStatus(stat)

         return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def filletBody(self, edges, radius):
        """
        Fillets the EDGEs on the source BODY Object (that has the type
        SOLIDBODY or SHEETBODY).

        Parameters
        ----------
        edges:  list of EDGE objects to fillet
        radius: the radius of the fillets created

        Returns
        -------
        body:
            the resultant BODY object (with the same type as the input
            source object)

        maps:
            list of Face mappings (in the result) which includes
            operations and an index to src where
            the Face originated - 2*nFaces in result in length
        """
        nedges = c_int(len(edges))
        edgs = (c_ego * nedges.value)()
        for i in range(nedges.value):
            edgs[i] = edges[i]._obj

        facemap = POINTER(c_int)()

        body = c_ego()
        stat = _egads.EG_filletBody(self._obj, nedges, edgs, c_double(radius), ctypes.byref(body),
                                    ctypes.byref(facemap))
        if stat: _raiseStatus(stat)

        nface = c_int()
        stat = _egads.EG_getBodyTopos(body, None, FACE, ctypes.byref(nface), None)

        maps = facemap[:2*nface.value]

        _egads.EG_free(facemap)

        return ego(body, self.context, deleteObject=True), maps

#=============================================================================-
    def chamferBody(self, edges, faces, dis1, dis2):
        """
        Chamfers the EDGEs on the source BODY Object (that has the type SOLIDBODY or SHEETBODY).

        Parameters
        ----------
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
        body:
            the resultant BODY object (with the same type as the input source object)

        maps:
            list of Face mappings (in the result) which includes operations and an
            index to src where the Face originated - 2*nFaces in result in length
        """
        if len(edges) != len(faces):
            raise EGADSError("edges and faces must be same length")

        nedges = len(edges)
        edgs = (c_ego * nedges)()
        facs = (c_ego * nedges)()
        for i in range(nedges):
            edgs[i] = edges[i]._obj
            facs[i] = faces[i]._obj

        facemap = POINTER(c_int)()

        body = c_ego()
        stat = _egads.EG_chamferBody(self._obj, c_int(nedges), edgs, facs, c_double(dis1), c_double(dis2), ctypes.byref(body),
                                     ctypes.byref(facemap))
        if stat: _raiseStatus(stat)

        nface = c_int()
        stat = _egads.EG_getBodyTopos(body, None, FACE, ctypes.byref(nface), None)

        maps = [None]*2*nface.value
        maps[:] = facemap[:2*nface.value]

        _egads.EG_free(facemap)

        return ego(body, self.context, deleteObject=True), maps

#=============================================================================-
    def hollowBody(self, faces, off, join):
        """
        A hollowed solid is built from an initial SOLIDBODY Object and a set of
        FACEs that initially bound the solid. These FACEs are removed and the
        remaining FACEs become the walls of the hollowed solid with the
        specified thickness. If there are no FACEs specified then the Body
        is offset by the specified distance (which can be negative).

        Parameters
        ----------
        faces:
            list of FACE objects to remove - (None performs an Offset)

        off:
            the wall thickness (offset) of the hollowed result

        join:
            0 - fillet-like corners, 1 - expanded corners

        Returns
        -------
        result:
            the resultant BODY object

        maps:
            list of Face mappings (in the result) which includes operations and an
            index to src where the Face originated - 2*nFaces in result in length
        """
        nfaces = 0
        pfaces = None
        if faces is not None:
            nfaces = len(faces)
            pfaces = (c_ego * nfaces)()
            for i in range(nfaces):
                pfaces[i] = faces[i]._obj

        facemap = POINTER(c_int)()

        body = c_ego()
        stat = _egads.EG_hollowBody(self._obj, c_int(nfaces), pfaces, c_double(off), c_int(join), ctypes.byref(body),
                                    ctypes.byref(facemap))
        if stat: _raiseStatus(stat)

        nface = c_int()
        stat = _egads.EG_getBodyTopos(body, None, FACE, ctypes.byref(nface), None)

        maps = [None]*2*nface.value
        maps[:] = facemap[:2*nface.value]

        _egads.EG_free(facemap)

        return ego(body, self.context, deleteObject=True), maps

#=============================================================================-
    def extrude(self, dist, dir):
        """
        Extrudes the source Object through the distance specified.
        If the Object is either a LOOP or WIREBODY the result is a
        SHEETBODY. If the source is either a FACE or FACEBODY then
        the returned Object is a SOLIDBODY.

        Parameters
        ----------
        dist:
            the distance to extrude

        dir:
            the vector that is the extrude direction (3 in length)

        Returns
        -------
        the resultant BODY object (type is one greater than the
        input source object)
        """
        direction = (c_double*3)()
        direction[:] = dir[:3]
        body = c_ego()
        stat = _egads.EG_extrude(self._obj, c_double(dist), direction, ctypes.byref(body))
        if stat: _raiseStatus(stat)

        return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def rotate(self, angle, axis):
        """
        Rotates the source Object about the axis through the angle
        specified. If the Object is either a LOOP or WIREBODY the
        result is a SHEETBODY. If the source is either a FACE or
        FACEBODY then the returned Object is a SOLIDBODY.

        Parameters
        ----------
        angle:
            the angle to rotate the object through [0-360 Degrees]

        axis:
            a point (on the axis) and a direction (two 3-tubles)

        Returns
        -------
        the resultant BODY object (type is one greater than the input
        source object)
        """
        _axis = (c_double*6)()
        _axis[0:3] = axis[0][:]
        _axis[3:6] = axis[1][:]

        body = c_ego()
        stat = _egads.EG_rotate(self._obj, c_double(angle), _axis, ctypes.byref(body))
        if stat: _raiseStatus(stat)

        return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def sweep(self, spline, mode):
        """
        Sweeps the source Object through the "spine" specified.
        The spine can be either an EDGE, LOOP or WIREBODY.
        If the source Object is either a LOOP or WIREBODY the result is a SHEETBODY.
        If the source is either a FACE or FACEBODY then the returned Object is a SOLIDBODY.

        Note: This does not always work as expected...

        Parameters
        ----------
        spline:
            the Object used as guide curve segment(s) to sweep the source through

        mode:
            0 - CorrectedFrenet      1 - Fixed
            2 - Frenet               3 - ConstantNormal
            4 - Darboux              5 - GuideAC
            6 - GuidePlan            7 - GuideACWithContact
            8 - GuidePlanWithContact 9 - DiscreteTrihedron
        """
        body = c_ego()
        stat = _egads.EG_sweep(self._obj, spline._obj, c_int(mode), ctypes.byref(body))
        if stat: _raiseStatus(stat)

        return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def initTessBody(self):
        """
        Creates an empty (open) discretization object for a Topological BODY Object.

        Returns
        -------
        resultant empty TESSELLATION object where each EDGE in the BODY must be
        filled via a call to EG_setTessEdge and each FACE must be filled with
        invocations of EG_setTessFace. The TESSSELLATION object is considered
        open until all EDGEs have been set (for a WIREBODY) or all FACEs have been
        set (for other Body types).
        """
        tess = c_ego()
        stat = _egads.EG_initTessBody(self._obj, ctypes.byref(tess))
        if stat: _raiseStatus(stat)

        return ego(tess, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def openTessBody(self):
        """
        Opens an existing Tessellation Object for replacing EDGE/FACE discretizations.
        """
        stat = _egads.EG_openTessBody(self._obj)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def setTessEdge(self, eIndex, xyzs, ts):
        """
        Sets the data associated with the discretization of an EDGE
        for an open Body-based Tessellation Object.

        Notes: (1) all vertices must be specified in increasing t.
               (2) the coordinates for the first and last vertex MUST match the
                   appropriate NODE's coordinates.
               (3) problems are reported to Standard Out regardless of the OutLevel.

        Parameters
        ----------
        eIndex:
            the EDGE index (1 bias). The EDGE Objects and number of EDGEs
            can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
            If this EDGE already has assigned data, it is overwritten.

        xyzs:
            the pointer to the set of coordinate data.

        ts:
            the pointer to the parameter values associated with each vertex.
        """
        npnt  = len(xyzs)
        pxyzs = (c_double*(3*npnt))()
        pts   = (c_double*npnt)()

        for i in range(npnt):
            pxyzs[3*i:3*i+3] = xyzs[i][:]

        pts[:] = ts[:]

        stat = _egads.EG_setTessEdge(self._obj, c_int(eIndex), c_int(npnt), pxyzs, pts)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def setTessFace(self, fIndex, xyz, uv, tris):
        """
        Sets the data associated with the discretization of a FACE for an open Body-based Tessellation Object.

        Parameters
        ----------
        fIndex:
            the FACE index (1 bias). The FACE Objects and number of FACEs can be
            retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
            If this FACE already has assigned data, it is overwritten.

        xyz:
            the pointer to the set of coordinate data for all vertices -- 3*len in length.

        uv:
            the pointer to the vertex parameter values -- 2*len in length.

        tris:
            the pointer to triangle vertex indices (1 bias) -- 3*ntri in length.
        """
        npnt  = len(xyz)
        pxyz  = (c_double*(3*npnt))()
        puv   = (c_double*(2*npnt))()
        ntri  = len(tris)
        ptris = (c_int*(3*ntri))()

        for i in range(npnt):
            pxyz[3*i:3*i+3] = xyz[i][:]
            puv[2*i:2*i+2]  = uv[i][:]

        for i in range(ntri):
            ptris[3*i:3*i+3] = tris[i][:]

        stat = _egads.EG_setTessFace(self._obj, c_int(fIndex), c_int(npnt), pxyz, puv, c_int(ntri), ptris)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def statusTessBody(self):
        """
        Returns the status of a TESSELLATION Object.

        Note: Placing the attribute ".mixed" on tess before invoking this
              function allows for tri/quad (2 tris) tessellations
              The type must be ATTRINT and the length is the number of FACEs,
              where the values are the number of quads

        Returns
        -------
        body:
            the returned associated BODY Object.

        stat:
            the state of the tessellation: -1 - closed but warned, 0 - open, 1 - OK, 2 - displaced.

        icode:
            EGADS_SUCCESS -- complete, EGADS_OUTSIDE -- still open.

        npts:
            the number of global points in the tessellation (0 -- open)
        """
        body = c_ego()
        stat = c_int()
        npts = c_int()

        icode = _egads.EG_statusTessBody(self._obj, ctypes.byref(body), ctypes.byref(stat), ctypes.byref(npts))

        return ego(body, self.context), stat.value, icode, npts.value

#=============================================================================-
    def finishTess(self, parms):
        """
        Completes the discretization for specified objects
        for the input TESSELLATION object.

        Note: an open TESSELLATION object is created by EG_initTessBody and
              can be partially filled via EG_setTessEdge and/or EG_setTessFace
              before this function is invoked.

        Parameters
        ----------
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
        params = (c_double*3)()
        params[:] = parms[:3]

        stat = _egads.EG_finishTess(self._obj, params)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def moveEdgeVert(self, eIndex, vIndex, t):
        """
        Moves the position of an EDGE vertex in a Body-based Tessellation Object.
        Will invalidate the Quad patches on any FACEs touching the EDGE.

        Parameters
        ----------
        eIndex:
            the EDGE index (1 bias).

        vIndex:
            the Vertex index in the EDGE (2 - nVert-1)

        t:
            the new parameter value on the EDGE for the point
        """
        stat = _egads.EG_moveEdgeVert(self._obj, c_int(eIndex), c_int(vIndex), c_double(t))
        if stat: _raiseStatus(stat)

#=============================================================================-
    def deleteEdgeVert(self, eIndex, vIndex, dir):
        """
        Deletes an EDGE vertex from a Body-based Tessellation Object.
        Will invalidate the Quad patches on any FACEs touching the EDGE.

        Parameters
        ----------
        eIndex:
            the EDGE index (1 bias).

        vIndex:
            the Vertex index in the EDGE (2 - nVert-1)

        dir:
            the direction to collapse any triangles (either -1 or 1)
        """
        stat = _egads.EG_deleteEdgeVert(self._obj, c_int(eIndex), c_int(vIndex), c_int(dir))
        if stat: _raiseStatus(stat)

#=============================================================================-
    def insertEdgeVerts(self, eIndex, vIndex, ts):
        """
        Inserts vertices into the EDGE discretization of a Body Tessellation Object.
        This will invalidate the Quad patches on any FACEs touching the EDGE.

        Parameters
        ----------
        eIndex:
            the EDGE index (1 bias).

        vIndex:
            the Vertex index in the EDGE (2 - nVert-1)

        ts:
            the t values for the new points. Must be monotonically increasing and
            be greater than the t of vIndex and less than the t of vIndex+1.
        """
        npts = len(ts)
        pts = (c_double*npts)()
        pts[:] = ts[:]

        stat = _egads.EG_insertEdgeVerts(self._obj, c_int(eIndex), c_int(vIndex), c_int(npts), pts)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def makeTessBody(self, parms):
        """
        Creates a discretization object from a Topological BODY Object.

        Parameters
        ----------
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
        the resulting egads tessellation object
        """
        params = (c_double*3)()
        params[:] = parms[:3]

        tess = c_ego()
        stat = _egads.EG_makeTessBody(self._obj, params, ctypes.byref(tess))
        if stat: _raiseStatus(stat)

        return ego(tess, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def remakeTess(self, facedg, parms):
        """
        Redoes the discretization for specified objects from within a BODY TESSELLATION.

        Parameters
        ----------
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
        nobj = len(facedg)
        pfacedg = (c_ego*nobj)()
        for i in range(nobj):
            pfacedg[i] = facedg[i]._obj

        params = (c_double*3)()
        params[:] = parms[:3]

        stat = _egads.EG_remakeTess(self._obj, c_int(nobj), pfacedg, params)
        if stat: _raiseStatus(stat)

#=============================================================================-
    def getEdgeUV(self, edge, t, sense=0):
        """
        Computes on the EDGE/PCURVE to get the appropriate [u,v] on the FACE.

        Parameters
        ----------
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
        the resulting [u,v] evaluated at t
        """
        puv = (c_double*2)()
        stat = _egads.EG_getEdgeUV(self._obj, edge._obj, c_int(sense), c_double(t), puv)
        if stat: _raiseStatus(stat)

        return tuple(puv[:])

#=============================================================================-
    def getTessEdge(self, eIndex):
        """
        Retrieves the data associated with the discretization of an EDGE
        from a Body-based Tessellation Object.

        Parameters
        ----------
        eIndex:
            the EDGE index (1 bias). The EDGE Objects and number of EDGEs
            can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
            A minus refers to the use of a mapped (+) Edge index from applying
            the functions mapBody and mapTessBody.

        Returns
        -------
        xyzs:
            the returned pointer to the set of coordinate data for each vertex

        ts:
            the returned pointer to the parameter values associated with each vertex
        """
        npnt = c_int()
        pxyz = POINTER(c_double)()
        pt = POINTER(c_double)()

        stat = _egads.EG_getTessEdge(self._obj, c_int(eIndex), ctypes.byref(npnt), ctypes.byref(pxyz), ctypes.byref(pt))
        if stat: _raiseStatus(stat)

        xyzs = [(pxyz[3*i], pxyz[3*i+1], pxyz[3*i+2]) for i in range(npnt.value)]
        ts   = list(pt[:npnt.value])

        return xyzs, ts

#=============================================================================-
    def getTessLoops(self, fIndex):
        """
        Retrieves the data for the LOOPs associated with the discretization of a FACE
        from a Body-based Tessellation Object.

        Parameters
        ----------
        fIndex: the FACE index (1 bias). The FACE Objects and number of FACEs
                can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.

        Returns
        -------
        the returned pointer to a vector of the last index (bias 1)
        for each LOOP (nloop in length).

        Notes: (1) all boundary vertices are listed first for any FACE tessellation,
               (2) outer LOOP is ordered in the counter-clockwise direction, and
               (3) inner LOOP(s) are ordered in the clockwise direction.
        """
        nloop = c_int()
        plIndex = POINTER(c_int)()

        stat = _egads.EG_getTessLoops(self._obj, c_int(fIndex),
                                      ctypes.byref(nloop), ctypes.byref(plIndex))
        if stat: _raiseStatus(stat)

        return list(plIndex[:nloop.value])

#=============================================================================-
    def getTessFace(self, fIndex):
        """
        Retrieves the data associated with the discretization of a FACE from a
        Body-based Tessellation Object.

        Parameters
        ----------
        fIndex:
            the FACE index (1 bias). The FACE Objects and number of FACEs
            can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.
            A minus refers to the use of a mapped (+) FACE index (if it exists).

        Returns
        -------
        xyz:
            set of coordinate data for each vertex

        uv:
            parameter values associated with each vertex

        ptype:
            vertex type (-1 - internal, 0 - NODE, >0 EDGE)

        pindex:
            vertex index (-1 internal)

        tris:
            returned pointer to triangle indices, 3 per triangle (1 bias)
            orientation consistent with FACE's mtype

        tric:
            returned pointer to neighbor information, 3 per triangle looking at opposing side:
            triangle (1-ntri), negative is Edge index for an external side
        """
        npnt = c_int()
        pxyz = POINTER(c_double)()
        puv = POINTER(c_double)()
        ptype = POINTER(c_int)()
        pindex = POINTER(c_int)()
        ntri = c_int()
        ptris = POINTER(c_int)()
        ptric = POINTER(c_int)()

        stat = _egads.EG_getTessFace(self._obj, c_int(fIndex), ctypes.byref(npnt), ctypes.byref(pxyz), ctypes.byref(puv),
                                     ctypes.byref(ptype), ctypes.byref(pindex),
                                     ctypes.byref(ntri), ctypes.byref(ptris), ctypes.byref(ptric))
        if stat: _raiseStatus(stat)

        xyz = [(pxyz[3*i], pxyz[3*i+1], pxyz[3*i+2]) for i in range(npnt.value)]
        uv  = [(puv[2*i], puv[2*i+1]) for i in range(npnt.value)]
        type = list(ptype[:npnt.value])
        index = list(pindex[:npnt.value])

        tris = [(ptris[3*i], ptris[3*i+1], ptris[3*i+2]) for i in range(ntri.value)]
        tric = [(ptric[3*i], ptric[3*i+1], ptric[3*i+2]) for i in range(ntri.value)]

        return xyz, uv, type, index, tris, tric

#=============================================================================-
    def localToGlobal(self, ind, locl):
        """
        Perform Local to Global index lookup. Tessellation Object must be closed.

        Parameters
        ----------
        ind:
            the topological index (1 bias) - 0 Node, (-) Edge, (+) Face

        locl:
            the local (or Node) index

        Returns
        -------
        the returned global index
        """
        gbl = c_int()

        stat = _egads.EG_localToGlobal(self._obj, c_int(ind), c_int(locl),
                                       ctypes.byref(gbl))
        if stat: _raiseStatus(stat)

        return gbl.value

#=============================================================================-
    def getGlobal(self, iglobal):
        """
        Returns the point type and index (like from EG_getTessFace)
          with coordinates.

        Parameters
        ----------
        iglobal:
            the global index (1 bias).

        Returns
        -------
        type:
            the point type (-) Face local index, (0) Node, (+) Edge local index

        index:
            the point topological index (1 bias)

        xyz:
            the coordinates at this global index (can be NULL for no return)
        """
        ptype = c_int()
        pindex = c_int()
        pxyz = (c_double*3)()

        stat = _egads.EG_getGlobal(self._obj, c_int(iglobal), ctypes.byref(ptype), ctypes.byref(pindex), pxyz)
        if stat: _raiseStatus(stat)

        return ptype.value, pindex.value, (pxyz[0], pxyz[1], pxyz[2])

#=============================================================================-
    def makeTessGeom(self, limits, sizes):
        """
        Creates a discretization object from a geometry-based Object.

        Parameters
        ----------
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
        """
        plimits = (c_double*4)()
        plimits[:] = limits[:4]
        psizes = (c_int*2)()
        psizes[:] = sizes[:2]

        tess = c_ego()
        stat = _egads.EG_makeTessGeom(self._obj, plimits, psizes, ctypes.byref(tess))
        if stat: _raiseStatus(stat)

        return ego(tess, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def getTessGeom(self):
        """
        Retrieves the data associated with the discretization of a geometry-based Object.

        Returns
        -------
        sizes:
            a returned set of 2 integers that specifies the size and dimensionality
            of the data. If the second is zero, then it is from a CURVE and
            the first integer is the length of the number of [x,y,z] triads.
            If the second integer is nonzero then the input data reflects a 2D map of coordinates.

        xyz:
            the returned pointer to the suite of coordinate data.
        """
        pxyz = POINTER(c_double)()
        psizes = (c_int*2)()

        tess = c_ego()
        stat = _egads.EG_getTessGeom(self._obj, psizes, ctypes.byref(pxyz))
        if stat: _raiseStatus(stat)

        sizes = tuple(psizes[:])
        xyz   = [(pxyz[3*i], pxyz[3*i+1], pxyz[3*i+2]) for i in range(sizes[0]*sizes[1])]

        return sizes, xyz

#=============================================================================-
    def mapBody(self, dst, fAttr):
        """
        Checks for topological equivalence between the the BODY self and the BODY dst.
        If necessary, produces a mapping (indices in src which map to dst) and
        places these as attributes on the resultant BODY mapped (named .nMap, .eMap and .fMap).
        Also may modify BSplines associated with FACEs.

        Parameters
        ----------
        dst:
            destination body object

        fAttr:
            the FACE attribute string used to map FACEs

        Returns
        -------
        the mapped resultant BODY object copied from dst
        If None and icode == EGADS_SUCCESS, dst is equivalent
        and can be used directly in EG_mapTessBody

        Note: It is the responsibility of the caller to have
              uniquely attributed all FACEs in both src and dst
              to aid in the mapping for all but FACEBODYs.
        """
        pfAttr = fAttr.encode() if isinstance(fAttr,str) else fAttr

        body = c_ego()
        stat = _egads.EG_mapBody(self._obj, dst._obj, pfAttr, ctypes.byref(body))
        if stat: _raiseStatus(stat)
        if not body: return None

        return ego(body, self.context, deleteObject=True)

#=============================================================================-
    def getWindingAngle(self, t):
        """
        Computes the Winding Angle along an Edge

        The Winding Angle is measured from one Face ``winding'' 
        around to the other based on the normals.
        An Edge with a single Face always returns 180.0.

        Parameters
        ----------
        t:
            the t-value along the Edge used to compute the Winding Angle

        Returns
        -------
        Winding Angle in degrees (0.0 -- 360.0)
        """
        angle = c_double()
        stat = _egads.EG_getWindingAngle(self._obj, c_double(t), ctypes.byref(angle))
        if stat: _raiseStatus(stat)
 
        return angle.value

#=============================================================================-
    def mapTessBody(self, body):
        """
        Maps the input discretization object to another BODY Object.
        The topologies of the BODY that created the input tessellation
        must match the topology of the body argument (the use of EG_mapBody
        can be used to assist).

        Parameters
        ----------
        body:
            the BODY object (with a matching Topology) used to map the
            tessellation.

        Returns
        -------
        the resultant TESSELLATION object. The triangulation is
        simply copied but the uv and xyz positions reflect
        the input body (above).

        Note: Invoking EG_moveEdgeVert, EG_deleteEdgeVert and/or EG_insertEdgeVerts
              in the source tessellation before calling this routine invalidates
              the ability of EG_mapTessBody to perform its function.
        """
        tess = c_ego()
        stat = _egads.EG_mapTessBody(self._obj, body._obj, ctypes.byref(tess))
        if stat: _raiseStatus(stat)

        return ego(tess, self.context, deleteObject=True, refs=body)

#=============================================================================-
    def locateTessBody(self, ifaces, uv, mapped=False):
        """
        Provides the triangle and the vertex weights for each of the input requests
        or the evaluated positions in a mapped tessellation.

        Parameters
        ----------
        ifaces:
            the face indices for each request - minus index refers to the use of a
            mapped Face index from EG_mapBody and EG_mapTessBody

        uv:
            the UV positions in the face for each request (2*len(ifaces) in length)

        mapped:
            perform mapped evaluations

        Returns
        -------
        weight:
            the vertex weights in the triangle that refer to the requested position
            (any negative weight indicates that the point was extrapolated)
            -or-
            the evaluated position based on the input uvs (when mapped is True)

        itris:
            the resultant 1-bias triangle index (if not mapped)
            if input as NULL then this function will perform mapped evaluations
        """
        npts = len(ifaces)
        pifaces = (c_int*npts)()
        pifaces[:] = ifaces[:]
        puv = (c_double*(2*npts))()
        for i in range(npts):
            puv[2*i:2*i+2] = uv[i][:]

        pitris = None if mapped else (c_int*npts)()

        pweight = (c_double*(3*npts))()

        stat = _egads.EG_locateTessBody(self._obj, c_int(npts), pifaces, puv, pitris, pweight)
        if stat: _raiseStatus(stat)

        weight = [(pweight[3*i], pweight[3*i+1], pweight[3*i+2]) for i in range(npts)]

        if mapped:
            return weight
        else:
            return weight, list(pitris[:])

#=============================================================================-
    def makeQuads(self, qparms, fIndex):
        """
        Creates Quadrilateral Patches for the indicated FACE and updates
        the Body-based Tessellation Object.

        Parameters
        ----------
        qparms:
            a set of 3 parameters that drive the Quadrilateral patching for the
            FACE. Any may be set to zero to indicate the use of the default value:
            parms[0] EDGE matching tolerance expressed as the deviation
                     from an aligned dot product [default: 0.05]
            parms[1] Maximum quad side ratio point count to allow
                    [default: 3.0]
            parms[2] Number of smoothing loops [default: 0.0]

        fIndex:
            the FACE index (1 bias)
        """
        pqparms = (c_double*3)()
        pqparms[:] = qparms[:3]

        stat = _egads.EG_makeQuads(self._obj, pqparms, c_int(fIndex))
        if stat: _raiseStatus(stat)

#=============================================================================-
    def getTessQuads(self):
        """
        Returns a list of Face indices found in the Body-based Tessellation Object
        that has been successfully Quadded via makeQuads.

        Returns
        -------
        fList:
            the returned pointer the Face indices (1 bias) - nface in length
            The Face Objects themselves can be retrieved via getBodyTopos
        """
        nface = c_int()
        pfList = POINTER(c_int)()

        stat = _egads.EG_getTessQuads(self._obj, ctypes.byref(nface), ctypes.byref(pfList))
        if stat: _raiseStatus(stat)

        fList = list(pfList[:nface.value])
        _egads.EG_free(pfList)

        return fList

#=============================================================================-
    def getQuads(self, fIndex):
        """
        Retrieves the data associated with the Quad-patching of a FACE
        from a Body-based Tessellation Object.

        Parameters
        ----------
        fIndex:
            the FACE index (1 bias). The FACE Objects and number of FACEs
            can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.

        Returns
        -------
        xyz:
            set of coordinate data for each vertex

        uv:
            parameter values associated with each vertex

        ptype:
            vertex type (-1 - internal, 0 - NODE, >0 EDGE)

        pindx:
            vertex index (-1 internal)

        npatch:
            number of patches
        """
        npnt = c_int()
        pxyz = POINTER(c_double)()
        puv = POINTER(c_double)()
        ptype = POINTER(c_int)()
        pindx = POINTER(c_int)()
        npatch = c_int()

        stat = _egads.EG_getQuads(self._obj, c_int(fIndex), ctypes.byref(npnt), ctypes.byref(pxyz), ctypes.byref(puv),
                                  ctypes.byref(ptype), ctypes.byref(pindx), ctypes.byref(npatch))
        if stat: _raiseStatus(stat)

        npnt = npnt.value

        xyz  = [(pxyz[3*i], pxyz[3*i+1], pxyz[3*i+2]) for i in range(npnt)]
        uv   = [(puv[2*i], puv[2*i+1]) for i in range(npnt)]
        type = list(ptype[:npnt])
        indx = list(pindx[:npnt])

        return xyz, uv, type, indx, npatch.value

#=============================================================================-
    def getPatch(self, fIndex, pIndex):
        """
        Retrieves the data associated with the Patch of a FACE
        from a Body-based Tessellation Object.

        Parameters
        ----------
        fIndex:
            the FACE index (1 bias). The FACE Objects and number of FACEs
            can be retrieved via EG_getBodyTopos and/or EG_indexBodyTopo.

        pIndex:
            the patch index (1-npatch from EG_getQuads)

        Returns
        -------
        n1:
            patch size in the first direction (indexed by i)

        n2:
            patch size in the second direction (indexed by j)

        pvindex:
            list of n1*n2 indices that define the patch

        pbounds:
            list of the neighbor bounding information for the patch
            (2*(n1-1)+2*(n2-1) in length). The first represents the segments at
            the base (j at base and increasing in i), the next is at the right (with i
            at max and j increasing). The third is the top (with j at max and i decreasing)
            and finally the left (i at min and j decreasing).
        """
        n1 = c_int()
        n2 = c_int()
        ppvindex = POINTER(c_int)()
        ppbounds = POINTER(c_int)()

        stat = _egads.EG_getPatch(self._obj, c_int(fIndex), c_int(pIndex),
                                  ctypes.byref(n1), ctypes.byref(n2),
                                  ctypes.byref(ppvindex), ctypes.byref(ppbounds))
        if stat: _raiseStatus(stat)

        n1 = n1.value
        n2 = n2.value

        pvindex = list(ppvindex[:n1*n2])
        pbounds = list(ppbounds[:2*(n1-1)+2*(n2-1)])

        return n1, n2, pvindex, pbounds

#=============================================================================-
    def quadTess(self):
        """
        Takes a triangulation as input and outputs a full quadrilateralization
        of the Body. The algorithm uses the bounds of each Face (the discretization
        of the Loops) and drives the interior towards regularization (4 quad sides
        attached to a vertex) without regard to spacing but maintaining a valid mesh.
        This is the recommended quad approach in that it is robust and does not require
        manual intervention like makeQuads (plus retrieving the quads is much simpler
        and does nor require invoking getQuads and getPatch).

        Returns
        -------
        a triangle-based Tessellation Object, but with pairs of triangles
        (as retrieved by calls to getTessFace) representing each quadrilateral.
        This is marked by the following attributes:
            ".tessType" (ATTRSTRING) is set to "Quad"
            ".mixed" with type ATTRINT and the length is the number of Faces in the Body,
                     where the values are the number of quads (triangle pairs) per Face
        """
        tess = c_ego()
        stat = _egads.EG_quadTess(self._obj, ctypes.byref(tess))
        if stat: _raiseStatus(stat)

        return ego(tess, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def tessMassProps(self):
        """
        Computes and returns the physical and inertial properties of a Tessellation Object.

        Returns
        -------
        volume, surface area (length for EDGE, LOOP or WIREBODY) center of gravity (3)
        inertia matrix at CoG (9)
        """
        lim = None
        data = (c_double*14)()

        stat = _egads.EG_tessMassProps(self._obj, data)

        volume = data[0]
        areaOrLen = data[1]
        CG = tuple(data[2:5])
        I = tuple(data[5:14])

        return volume, areaOrLen, CG, I

#=============================================================================-
    def initEBody(self, angle):
        """
        Takes as input a Body Tessellation Object and returns an Open EBody fully initialized with Effective}
        Objects (that may only contain a single non-effective object). EEdges are automatically merged where
        possible (removing Nodes that touch 2 Edges, unless degenerate or marked as ".Keep"). The tessellation
        is used (and required) in order that single UV space be constructed for EFaces.

        Parameters
        ----------
        angle:
            angle used to determine if Nodes on open
            Edges of Sheet Bodies can be removed. The dot of the tangents at the Node is
            compared to this angle (in degrees). If the dot is greater than the angle,
            the Node does not disappear. The angle is also used to test Edges to see if
            they can be removed. Edges with normals on both trimmed Faces showing deviations
            greater than the input angle will not disappear.

        Returns
        -------
        the resultant Open Effective Topology Body Object
        """
        ebody = c_ego()

        stat = _egads.EG_initEBody(self._obj, c_double(angle), ctypes.byref(ebody))
        if stat: _raiseStatus(stat)

        # the EBody references this tessellation until the EBody is finished
        return ego(ebody, self.context, deleteObject=True, refs=self)

#=============================================================================-
    def finishEBody(self):
        """
        Finish an EBody
        """
        stat = _egads.EG_finishEBody(self._obj)
        if stat: _raiseStatus(stat)

        # the EBody now references the body referenced by the tessellation
        self._refs = self._refs._refs

#=============================================================================-
    def makeAttrEFaces(self, attrName):
        """
        Modifies the EBody by finding "free" Faces (a single Face in an EFace)
        with attrName and the same attribute value(s), thus making a collection of EFaces.
        The single attrName is placed on the EFace(s) unless in "Full Attribute" mode, which
        then performs attribute merging. This function returns the number of EFaces
        possibly constructed (neface). The UVbox can be retrieved via calls to either
        getTopology or getRange with the returned appropriate efaces.


        Parameters
        ----------
        attrName:
            the attribute name used to collect Faces into an EFaces. The attribute value(s)
            are then matched to make the collections.

        Returns
        -------
        the list of EFaces constructed
        """
        attrName = attrName.encode() if isinstance(attrName,str) else attrName

        nefaces = c_int()
        pefaces = POINTER(c_ego)()
        stat = _egads.EG_makeAttrEFaces(self._obj, attrName, ctypes.byref(nefaces), ctypes.byref(pefaces))
        if stat: _raiseStatus(stat)

        efaces = [ego(c_ego(pefaces[i].contents), self.context) for i in range(nefaces.value)]
        _egads.EG_free(pefaces)

        return efaces

#=============================================================================-
    def makeEFace(self, faces):
        """
        Modifies the EBody by removing the nface "free" Faces and replacing them
        with a single eface (returned for convenience, and note that ebody has been updated).
        There are no attributes set on eface unless the "Full Attribution" mode is in place.
        This constructs a single UV for the faces. The UVbox can be retrieved via calls to either
        getTopology or getRange with eface.


        Parameters
        ----------
        faces:
            the list of Face Objects used to make eface.

        Returns
        -------
        the EFace constructed
        """

        nfaces = len(faces)
        pfaces = (c_ego * nfaces)()
        for i in range(nfaces):
            pfaces[i] = faces[i]._obj

        peface = c_ego()
        stat = _egads.EG_makeEFace(self._obj, c_int(nfaces), pfaces, ctypes.byref(peface))
        if stat: _raiseStatus(stat)

        return ego(peface, self.context)

#=============================================================================-
    def effectiveMap(self, eparam):
        """
        Returns the evaluated location in the BRep for the Effective Topology entity.

        Parameters
        ----------
        eparam:
            t for EEdges / the [u, v] for an EFace

        Returns
        -------
        obj:
            the returned source Object in the Body

        param:
            the returned t for an Edge / the returned [u, v] for the Face
        """
        peparam = (c_double * 2)()
        param = (c_double * 2)()
        oclass = self._obj.contents.oclass

        try:
            if (oclass == EEDGE):
                if isinstance(eparam, (int, float)):
                    peparam[0] = eparam
                else:
                    peparam[0] = eparam[0]
            elif (oclass == EFACE):
                peparam[0] = eparam[0]
                peparam[1] = eparam[1]
        except:
            errmsg = 'Failed to convert parameter value: eparam = ' + str(eparam)
            raise ValueError(errmsg)

        pobj = c_ego()
        stat = _egads.EG_effectiveMap(self._obj, peparam, ctypes.byref(pobj), param)
        if stat: _raiseStatus(stat)

        if (oclass == EEDGE):
            param = peparam[0]
        elif (oclass == EFACE):
            param = [peparam[0], peparam[1]]

        return ego(pobj, self.context), param

#=============================================================================-
    def effectiveEdgeList(self):
        """
        Returns the list of Edge entities in the source Body that make up the EEdge.
        A pointer to an integer list of senses for each Edge is returned as well
        as the starting t value in the EEdge (remember that the t will go in the
        opposite direction in the Edge if the sense is SREVERSE).

        Returns
        -------
        edges:
            list of Edges - nedge in length

        senses:
            list of senses - nedge in length

        tstart:
            list of t starting values - nedge in length
       """
        nedges = c_int()
        pedges = POINTER(c_ego)()
        psenses = POINTER(c_int)()
        ptstart = POINTER(c_double)()
        stat = _egads.EG_effectiveEdgeList(self._obj, ctypes.byref(nedges), ctypes.byref(pedges),
                                                      ctypes.byref(psenses), ctypes.byref(ptstart))
        if stat: _raiseStatus(stat)

        edges = [ego(c_ego(pedges[i].contents), self.context) for i in range(nedges.value)]
        senses = [psenses[i] for i in range(nedges.value)]
        tstart = [ptstart[i] for i in range(nedges.value)]

        _egads.EG_free(pedges)
        _egads.EG_free(psenses)
        _egads.EG_free(ptstart)

        return edges, senses, tstart

#=============================================================================-
def makeLoop(edges, geom=None, toler=0.0):
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
        Note: list entries are None-ified when included in LOOPs

    geom:
        SURFACE Object for non-planar LOOPs to be used to bound
        FACEs (can be None)

    toler:
        tolerance used for the operation (0.0 - use EDGE tolerances)

    Returns
    -------
    loop:
        the resultant LOOP Object

    nedges:
        the number of non-None entries in edges when returned
        or error code
    """
    nedges = len(edges)
    edgs = (c_ego * nedges)()
    for i in range(nedges):
        if edges[i] is None:
            edgs[i] = None
        else:
            edgs[i] = edges[i]._obj
            context = edges[i].context

    geo = None if geom is None else geom._obj

    loop = c_ego()
    nloop_edges = _egads.EG_makeLoop(c_int(nedges), edgs, geo, c_double(toler), ctypes.byref(loop))

    if nloop_edges < EGADS_SUCCESS:
        _raiseStatus(nloop_edges)

    refs = []
    for i in range(nedges):
        if not edgs[i]:
            refs.append(edges[i])
            edges[i] = None

    return ego(loop, context, deleteObject=True, refs=refs), edges

#=============================================================================-
def sewFaces(objlist, toler=0.0, manifold=True):
    """
    Creates a MODEL from a collection of Objects. The Objects can
    be either BODYs (not WIREBODY), SHELLs and/or FACEs. After the
    sewing operation, any unconnected Objects are returned as
    BODYs.

    Parameters
    ----------
    objlist:
        List of Objects to sew together

    toler:
        Tolerance used for the operation (0.0 - use Face tolerances)

    manifold:
        Indicates whether to produce manifold/non-manifold geometry

    Returns
    -------
    The resultant MODEL object
    """
    flag = c_int(1)
    if manifold:
        flag = 0

    nobj = len(objlist)
    objs = (c_ego * nobj)()
    for i in range(nobj):
        objs[i] = objlist[i]._obj

    model = c_ego()
    stat = _egads.EG_sewFaces(c_int(nobj), objs, toler, flag, ctypes.byref(model))
    if stat: _raiseStatus(stat)

    return ego(model, objlist[0].context, deleteObject=True, refs=objlist)

#=============================================================================-
def blend(sections, rc1=None, rc2=None):
    """
    Simply lofts the input Objects to create a BODY Object
    (that has the type SOLIDBODY or SHEETBODY). Cubic BSplines are used.
    All sections must have the same number of Edges (except for NODEs)
    and the Edge order in each (defined in a CCW manner) is used to
    specify the loft connectivity.

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

    rc1:
        specifies treatment* at the first section (or None for no treatment)

    rc2:
        specifies treatment* at the first section (or None for no treatment)

    * for NODEs -- elliptical treatment (8 in length): radius of
    curvature1, unit direction, rc2, orthogonal direction;
    nSection must be at least 3 (or 4 for treatments at both ends)
    for other sections -- setting tangency (4 in length):
    magnitude, unit direction for FACEs with 2 or 3 EDGEs -- make
    a Wing Tip-like cap: zero, growthFactor (len of 2)

    Returns
    -------
    the resultant BODY object
    """
    nsec = len(sections)
    secs = (c_ego * nsec)()
    for i in range(nsec):
        secs[i] = sections[i]._obj

    prc1 = None
    if rc1 is not None:
        prc1 = (c_double * 8)(0.0)
        for i in range(min(len(rc1), 8)):
            prc1[i] = rc1[i]

    prc2 = None
    if rc2 is not None:
        prc2 = (c_double * 8)(0.0)
        for i in range(min(len(rc2), 8)):
            prc2[i] = rc2[i]

    body = c_ego()
    stat = _egads.EG_blend(nsec, secs, prc1, prc2, ctypes.byref(body))
    if stat: _raiseStatus(stat)

    return ego(body, sections[0].context, deleteObject=True)

#=============================================================================-
def ruled(sections):
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
    the resultant BODY object
    """
    nsec = len(sections)
    secs = (c_ego * nsec)()
    for i in range(nsec):
        secs[i] = sections[i]._obj

    body = c_ego()
    stat = _egads.EG_ruled(nsec, secs, ctypes.byref(body))
    if stat: _raiseStatus(stat)

    return ego(body, sections[0].context, deleteObject=True)

#=============================================================================-
def skinning(curves, degree=3):
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
    The new BSpline Surface Object
    """
    nCurves = len(curves)
    pcurves = (c_ego * nCurves)()
    for i in range(nCurves):
        pcurves[i] = curves[i]._obj

    geom = c_ego()
    stat = _egads.EG_skinning(c_int(nCurves), pcurves, c_int(degree),
                              ctypes.byref(geom))
    if stat: _raiseStatus(stat)

    return ego(geom, curves[0].context, deleteObject=True)
