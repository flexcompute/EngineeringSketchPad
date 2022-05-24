###########################################################################
#                                                                         #
#   pyCAPS --- Python version of CAPS API                                 #
#                                                                         #
#                                                                         #
#      Copyright 2011-2022, Massachusetts Institute of Technology         #
#      Licensed under The GNU Lesser General Public License, version 2.1  #
#      See http://www.opensource.org/licenses/lgpl-2.1.php                #
#                                                                         #
###########################################################################

import ctypes
from ctypes import POINTER, c_short, c_int, c_ulong, c_ulonglong, c_double, c_void_p, c_char_p
import os
import sys
import json
import copy
import signal
import functools
from sys import version_info
#import xml.etree.ElementTree as xmlElementTree

from pyEGADS import egads

try:
    import numpy
except ImportError:
    numpy = None

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

# load the shared library
if sys.platform.startswith('darwin'):
    _caps = ctypes.CDLL(_ESP_ROOT + "/lib/libcaps.dylib")
elif sys.platform.startswith('linux'):
    _caps = ctypes.CDLL(_ESP_ROOT + "/lib/libcaps.so")
elif sys.platform.startswith('win32'):
    if version_info.major == 3 and version_info.minor < 8:
        _caps = ctypes.CDLL(_ESP_ROOT + "\\lib\\caps.dll")
    else:
        _caps = ctypes.CDLL(_ESP_ROOT + "\\lib\\caps.dll", winmode=0)
else:
    raise IOError("Unknown platform: " + sys.platform)

__all__ = []

# =============================================================================
if sys.platform.startswith('win32'):
    CAPSLONG = c_ulonglong
else:
    CAPSLONG = c_ulong

# =============================================================================
# alias for clarification
enum_capsoType   = c_int
enum_capssType   = c_int
enum_capsvType   = c_int
enum_capsfType   = c_int
enum_capsFixed   = c_int
enum_capsNull    = c_int
enum_capstMethod = c_int
enum_capsdMethod = c_int
enum_capsState   = c_int

# =============================================================================
# define c_capsTuple structure
class c_capsTuple(ctypes.Structure):
    pass

c_capsTuple.__slots__ = [
    'name',
    'value',
]
c_capsTuple._fields_ = [
    ('name', c_char_p),
    ('value', c_char_p),
]

# =============================================================================
# define c_capsOwn structure
class c_capsOwn(ctypes.Structure):
    pass

c_capsOwn.__slots__ = [
    'index',
    'pname',
    'pID',
    'user',
    'datetime',
    'sNum',
]
c_capsOwn._fields_ = [
    ('index', c_int),
    ('pname', c_char_p),
    ('pID', c_char_p),
    ('user', c_char_p),
    ('datetime', c_short * int(6)),
    ('sNum', CAPSLONG),
]

# =============================================================================
# define c_capsObject structure
class c_capsObject(ctypes.Structure):
    pass

c_capsObject.__slots__ = [
    'magicnumber',
    'type',
    'subtype',
    'delMark',
    'name',
    'attrs',
    'blind',
    'flist',
    'nHistory',
    'history',
    'last',
    'parent',
]
c_capsObject._fields_ = [
    ('magicnumber', c_int),
    ('type', c_int),
    ('subtype', c_int),
    ('delMark', c_int),
    ('name', c_char_p),
    ('attrs', POINTER(egads.c_egAttrs)),
    ('blind', c_void_p),
    ('flist', c_void_p),
    ('nHistory', c_int),
    ('history', POINTER(c_capsOwn)),
    ('last', c_capsOwn),
    ('parent', POINTER(c_capsObject)),
]

# define the c_capsObj pointer
c_capsObj = POINTER(c_capsObject)

# =============================================================================
# define c_capsError structure
class c_capsError(ctypes.Structure):
    pass

c_capsError.__slots__ = [
    'errObj',
    'eType',
    'index',
    'nLines',
    'lines',
]
c_capsError._fields_ = [
    ('errObj', POINTER(c_capsObject)),
    ('eType', c_int),
    ('index', c_int),
    ('nLines', c_int),
    ('lines', POINTER(c_char_p)),
]

# =============================================================================
# define c_capsErrs structure
class c_capsErrs(ctypes.Structure):
    pass

c_capsErrs.__slots__ = [
    'nError',
    'errors',
]
c_capsErrs._fields_ = [
    ('nError', c_int),
    ('errors', POINTER(c_capsError)),
]

# =============================================================================
# define c_capsDeriv structure
class c_capsDeriv(ctypes.Structure):
    pass

c_capsDeriv.__slots__ = [
    'name',
    'rank',
    'deriv',
]
c_capsDeriv._fields_ = [
    ('name', c_char_p),
    ('rank', c_int),
    ('deriv', POINTER(c_double)),
]

# =============================================================================
# define c_capsValue structure
class union_capsValue_vals(ctypes.Union):
    pass

union_capsValue_vals.__slots__ = [
    'integer',
    'integers',
    'real',
    'reals',
    'string',
    'tuple',
    'AIMptr',
]
union_capsValue_vals._fields_ = [
    ('integer', c_int),
    ('integers', POINTER(c_int)),
    ('real', c_double),
    ('reals', POINTER(c_double)),
    ('string', c_char_p),
    ('tuple', POINTER(c_capsTuple)),
    ('AIMptr', POINTER(None)),
]

class union_capsValue_limits(ctypes.Union):
    pass

union_capsValue_limits.__slots__ = [
    'ilims',
    'dlims',
]
union_capsValue_limits._fields_ = [
    ('ilims', c_int * int(2)),
    ('dlims', c_double * int(2)),
]

class c_capsValue(ctypes.Structure):
    pass

c_capsValue.__slots__ = [
    'type',
    'length',
    'dim',
    'nrow',
    'ncol',
    'lfixed',
    'sfixed',
    'nullVal',
    'index',
    'pIndex',
    'gInType',
    'vals',
    'limits',
    'units',
    'link',
    'linkMethod',
    'partial',
    'nderiv',
    'derivs',
]
c_capsValue._fields_ = [
    ('type', c_int),
    ('length', c_int),
    ('dim', c_int),
    ('nrow', c_int),
    ('ncol', c_int),
    ('lfixed', c_int),
    ('sfixed', c_int),
    ('nullVal', c_int),
    ('index', c_int),
    ('pIndex', c_int),
    ('gInType', c_int),
    ('vals', union_capsValue_vals),
    ('limits', union_capsValue_limits),
    ('units', c_char_p),
    ('link', POINTER(c_capsObject)),
    ('linkMethod', c_int),
    ('partial', POINTER(c_int)),
    ('nderiv', c_int),
    ('derivs', POINTER(c_capsDeriv)),
]


# =============================================================================
#
# caps.h functions
#
# =============================================================================

# base-level object functions

_caps.caps_revision.argtypes = [POINTER(c_int), POINTER(c_int)]
_caps.caps_revision.restype = None

_caps.caps_info.argtypes = [c_capsObj, POINTER(c_char_p), POINTER(enum_capsoType), POINTER(enum_capssType), POINTER(c_capsObj), POINTER(c_capsObj), POINTER(c_capsOwn)]
_caps.caps_info.restype = c_int

_caps.caps_size.argtypes = [c_capsObj, enum_capsoType, enum_capssType, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_size.restype = c_int

_caps.caps_childByIndex.argtypes = [c_capsObj, enum_capsoType, enum_capssType, c_int, POINTER(c_capsObj)]
_caps.caps_childByIndex.restype = c_int

_caps.caps_childByName.argtypes = [c_capsObj, enum_capsoType, enum_capssType, c_char_p, POINTER(c_capsObj), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_childByName.restype = c_int

_caps.caps_bodyByIndex.argtypes = [c_capsObj, c_int, POINTER(egads.c_ego), POINTER(c_char_p)]
_caps.caps_bodyByIndex.restype = c_int

_caps.caps_ownerInfo.argtypes = [c_capsObj, c_capsOwn, POINTER(c_char_p), 
                                 POINTER(c_char_p), POINTER(c_char_p), POINTER(c_char_p), POINTER(c_int), 
                                 POINTER(POINTER(c_char_p)), POINTER(c_short), POINTER(CAPSLONG)]
_caps.caps_ownerInfo.restype = c_int

# _caps.caps_getHistory.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsOwn))]
# _caps.caps_getHistory.restype = c_int

_caps.caps_delete.argtypes = [c_capsObj]
_caps.caps_delete.restype = c_int

_caps.caps_errorInfo.argtypes = [POINTER(c_capsErrs), c_int, POINTER(POINTER(c_capsObj)), POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_char_p))]
_caps.caps_errorInfo.restype = c_int

_caps.caps_freeError.argtypes = [POINTER(c_capsErrs)]
_caps.caps_freeError.restype = c_int

_caps.caps_freeValue.argtypes = [POINTER(c_capsValue)]
_caps.caps_freeValue.restype = None

# I/O functions 

_caps.caps_writeParameters.argtypes = [c_capsObj, c_char_p]
_caps.caps_writeParameters.restype = c_int

_caps.caps_readParameters.argtypes = [c_capsObj, c_char_p]
_caps.caps_readParameters.restype = c_int

_caps.caps_writeGeometry.argtypes = [c_capsObj, c_int, c_char_p, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_writeGeometry.restype = c_int

# attribute functions 

_caps.caps_attrByName.argtypes = [c_capsObj, c_char_p, POINTER(c_capsObj)]
_caps.caps_attrByName.restype = c_int

_caps.caps_attrByIndex.argtypes = [c_capsObj, c_int, POINTER(c_capsObj)]
_caps.caps_attrByIndex.restype = c_int

_caps.caps_setAttr.argtypes = [c_capsObj, c_char_p, c_capsObj]
_caps.caps_setAttr.restype = c_int

_caps.caps_deleteAttr.argtypes = [c_capsObj, c_char_p]
_caps.caps_deleteAttr.restype = c_int

# problem functions 

_caps.caps_phaseState.argtypes = [c_char_p, c_char_p, POINTER(c_int)]
_caps.caps_phaseState.restype = c_int

_caps.caps_journalState.argtypes = [c_capsObj]
_caps.caps_journalState.restype = c_int

_caps.caps_open.argtypes = [c_char_p, c_char_p, c_int, c_void_p, c_int, POINTER(c_capsObj), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_open.restype = c_int

_caps.caps_close.argtypes = [c_capsObj, c_int, c_char_p]
_caps.caps_close.restype = c_int

_caps.caps_outLevel.argtypes = [c_capsObj, c_int]
_caps.caps_outLevel.restype = c_int

_caps.caps_getRootPath.argtypes = [c_capsObj, POINTER(c_char_p)]
_caps.caps_getRootPath.restype = c_int

_caps.caps_intentPhrase.argtypes = [c_capsObj, c_int, POINTER(c_char_p)]
_caps.caps_intentPhrase.restype = c_int

# analysis functions 

_caps.caps_queryAnalysis.argtypes = [c_capsObj, c_char_p, POINTER(c_int), POINTER(c_int), POINTER(c_int)]
_caps.caps_queryAnalysis.restype = c_int

_caps.caps_getBodies.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(egads.c_ego)), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_getBodies.restype = c_int

_caps.caps_getInput.argtypes = [c_capsObj, c_char_p, c_int, POINTER(c_char_p), POINTER(c_capsValue)]
_caps.caps_getInput.restype = c_int

_caps.caps_getOutput.argtypes = [c_capsObj, c_char_p, c_int, POINTER(c_char_p), POINTER(c_capsValue)]
_caps.caps_getOutput.restype = c_int

_caps.caps_AIMbackdoor.argtypes = [c_capsObj, c_char_p, POINTER(c_char_p)]
_caps.caps_AIMbackdoor.restype = c_int

_caps.caps_makeAnalysis.argtypes = [c_capsObj, c_char_p, c_char_p, c_char_p, c_char_p, POINTER(c_int), POINTER(c_capsObj), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_makeAnalysis.restype = c_int

_caps.caps_dupAnalysis.argtypes = [c_capsObj, c_char_p, POINTER(c_capsObj)]
_caps.caps_dupAnalysis.restype = c_int

_caps.caps_resetAnalysis.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_resetAnalysis.restype = c_int

_caps.caps_dirtyAnalysis.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsObj))]
_caps.caps_dirtyAnalysis.restype = c_int

_caps.caps_analysisInfo.argtypes = [c_capsObj, POINTER(c_char_p), POINTER(c_char_p), POINTER(c_int), POINTER(c_int), POINTER(c_char_p), POINTER(c_int), POINTER(POINTER(c_char_p)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int), POINTER(c_int)]
_caps.caps_analysisInfo.restype = c_int

_caps.caps_preAnalysis.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_preAnalysis.restype = c_int

_caps.caps_system.argtypes = [c_capsObj, c_char_p, c_char_p]
_caps.caps_system.restype = c_int

_caps.caps_execute.argtypes = [c_capsObj, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_execute.restype = c_int

_caps.caps_postAnalysis.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_postAnalysis.restype = c_int

# bound, vertexset and dataset functions

_caps.caps_makeBound.argtypes = [c_capsObj, c_int, c_char_p, POINTER(c_capsObj)]
_caps.caps_makeBound.restype = c_int

_caps.caps_boundInfo.argtypes = [c_capsObj, POINTER(enum_capsState), POINTER(c_int), POINTER(c_double)]
_caps.caps_boundInfo.restype = c_int

_caps.caps_closeBound.argtypes = [c_capsObj]
_caps.caps_closeBound.restype = c_int

_caps.caps_makeVertexSet.argtypes = [c_capsObj, c_capsObj, c_char_p, POINTER(c_capsObj), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_makeVertexSet.restype = c_int

_caps.caps_vertexSetInfo.argtypes = [c_capsObj, POINTER(c_int), POINTER(c_int), POINTER(c_capsObj), POINTER(c_capsObj)]
_caps.caps_vertexSetInfo.restype = c_int

_caps.caps_outputVertexSet.argtypes = [c_capsObj, c_char_p]
_caps.caps_outputVertexSet.restype = c_int

_caps.caps_fillUnVertexSet.argtypes = [c_capsObj, c_int, POINTER(c_double)]
_caps.caps_fillUnVertexSet.restype = c_int

_caps.caps_makeDataSet.argtypes = [c_capsObj, c_char_p, enum_capsfType, c_int, POINTER(c_capsObj), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_makeDataSet.restype = c_int

_caps.caps_dataSetInfo.argtypes = [c_capsObj, POINTER(enum_capsfType), POINTER(c_capsObj), POINTER(enum_capsdMethod)]
_caps.caps_dataSetInfo.restype = c_int

_caps.caps_linkDataSet.argtypes = [c_capsObj, enum_capsdMethod, c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))];
_caps.caps_linkDataSet.restype = c_int

_caps.caps_initDataSet.argtypes = [c_capsObj, c_int, POINTER(c_double), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_initDataSet.restype = c_int

_caps.caps_setData.argtypes = [c_capsObj, c_int, c_int, POINTER(c_double), c_char_p, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_setData.restype = c_int

_caps.caps_getData.argtypes = [c_capsObj, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_char_p), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_getData.restype = c_int

_caps.caps_getDataSets.argtypes = [c_capsObj, c_char_p, POINTER(c_int), POINTER(POINTER(c_capsObj))]
_caps.caps_getDataSets.restype = c_int

_caps.caps_getTriangles.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_int)), 
                                               POINTER(c_int), POINTER(POINTER(c_int)),
                                               POINTER(c_int), POINTER(POINTER(c_int)),
                                               POINTER(c_int), POINTER(POINTER(c_int))]
_caps.caps_getTriangles.restype = c_int

# value functions 

_caps.caps_getValue.argtypes = [c_capsObj, POINTER(enum_capsvType), POINTER(c_int), POINTER(c_int), POINTER(POINTER(None)), POINTER(POINTER(c_int)), POINTER(c_char_p), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_getValue.restype = c_int

_caps.caps_makeValue.argtypes = [c_capsObj, c_char_p, enum_capssType, enum_capsvType, c_int, c_int, POINTER(None), POINTER(c_int), c_char_p, POINTER(c_capsObj)]
_caps.caps_makeValue.restype = c_int

_caps.caps_setValue.argtypes = [c_capsObj, enum_capsvType, c_int, c_int, POINTER(None), POINTER(c_int), c_char_p, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_setValue.restype = c_int

_caps.caps_getLimits.argtypes = [c_capsObj, POINTER(enum_capsvType), POINTER(POINTER(None)), POINTER(c_char_p)]
_caps.caps_getLimits.restype = c_int

_caps.caps_setLimits.argtypes = [c_capsObj, enum_capsvType, POINTER(None), c_char_p, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_setLimits.restype = c_int

_caps.caps_getValueProps.argtypes = [c_capsObj, POINTER(c_int), POINTER(c_int), POINTER(enum_capsFixed), POINTER(enum_capsFixed), POINTER(enum_capsNull)]
_caps.caps_getValueProps.restype = c_int

_caps.caps_setValueProps.argtypes = [c_capsObj, c_int, enum_capsFixed, enum_capsFixed, enum_capsNull, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_setValueProps.restype = c_int

_caps.caps_convertValue.argtypes = [c_capsObj, c_double, c_char_p, POINTER(c_double)]
_caps.caps_convertValue.restype = c_int

_caps.caps_transferValues.argtypes = [c_capsObj, enum_capstMethod, c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_transferValues.restype = c_int

_caps.caps_linkValue.argtypes = [c_capsObj, enum_capstMethod, c_capsObj, POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_linkValue.restype = c_int

_caps.caps_hasDeriv.argtypes = [c_capsObj, POINTER(c_int), POINTER(POINTER(c_char_p)), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_hasDeriv.restype = c_int

_caps.caps_getDeriv.argtypes = [c_capsObj, c_char_p, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_double)), POINTER(c_int), POINTER(POINTER(c_capsErrs))]
_caps.caps_getDeriv.restype = c_int

# units

_caps.caps_convert.argtypes = [c_int, c_char_p, POINTER(c_double), c_char_p, POINTER(c_double)]
_caps.caps_convert.restype = c_int

_caps.caps_unitParse.argtypes = [c_char_p]
_caps.caps_unitParse.restype = c_int

_caps.caps_unitConvertible.argtypes = [c_char_p, c_char_p]
_caps.caps_unitConvertible.restype = c_int

_caps.caps_unitCompare.argtypes = [c_char_p, c_char_p, POINTER(c_int)]
_caps.caps_unitCompare.restype = c_int

_caps.caps_unitMultiply.argtypes = [c_char_p, c_char_p, POINTER(c_char_p)]
_caps.caps_unitMultiply.restype = c_int

_caps.caps_unitDivide.argtypes = [c_char_p, c_char_p, POINTER(c_char_p)]
_caps.caps_unitDivide.restype = c_int

_caps.caps_unitRaise.argtypes = [c_char_p, c_int, POINTER(c_char_p)]
_caps.caps_unitRaise.restype = c_int

_caps.caps_unitOffset.argtypes = [c_char_p, c_double, POINTER(c_char_p)]
_caps.caps_unitOffset.restype = c_int

# others

_caps.caps_externSignal.argtypes = []
_caps.caps_externSignal.restype = None

_caps.caps_rmLock.argtypes = []
_caps.caps_rmLock.restype = None

_caps.caps_printObjects.argtypes = [c_capsObj, c_capsObj, c_int]
_caps.caps_printObjects.restype = None

_caps.caps_outputObjects.argtypes = [c_capsObj, POINTER(c_char_p)]
_caps.caps_outputObjects.restype = c_int


# =============================================================================

# Extract CAPS error codes
_caps_error_codes = {}
globalDict = globals()
minErr = 0
with open(os.path.join(_ESP_ROOT,"include","capsErrors.h")) as fp:
    lines = fp.readlines()
    for line in lines:
        define = line.split()
        if len(define) < 3: continue
        if define[0] != "#define": continue
        _caps_error_codes[int(define[2])] = define[1]
        globalDict[define[1]] = int(define[2])
        __all__.append(define[1])
        minErr = min(int(define[2]), minErr)
globalDict["CAPS_CLOSED"] = minErr
del fp
del lines
del define
del globalDict
del minErr

# Extract CAPS enums
class oFlag:   
    oFileName  = 0
    oMODL      = 1
    oEGO       = 2
    oPhaseName = 3
    oContinue  = 4
    oReload    = 5

__all__.append("oType")
class oType:
     BODIES     = -2
     ATTRIBUTES = -1
#    UNUSED     =  0
     PROBLEM    =  1
     VALUE      =  2
     ANALYSIS   =  3
     BOUND      =  4
     VERTEXSET  =  5
     DATASET    =  6

__all__.append("sType")
class sType:
    NONE         = 0
    STATIC       = 1
    PARAMETRIC   = 2
    GEOMETRYIN   = 3
    GEOMETRYOUT  = 4
    PARAMETER    = 5
    USER         = 6
    ANALYSISIN   = 7
    ANALYSISOUT  = 8
    CONNECTED    = 9
    UNCONNECTED  = 10
    ANALYSISDYNO = 11

__all__.append("eType")
class eType:
    CONTINUATION = -1
    CINFO        =  0
    CWARN        =  1
    CERROR       =  2
    CSTAT        =  3

__all__.append("fType")
class fType:
    FieldIn     = 0
    FieldOut    = 1
    GeomSens    = 2
    TessSens    = 3
    User        = 4
    BuiltIn     = 5

__all__.append("vType")
class vType:
    Boolean   = 0
    Integer   = 1
    Double    = 2
    String    = 3
    Tuple     = 4
    Value     = 5
    DoubleDeriv = 6
    # extra value for Python API
    JsonDict  = 7
    Null      = 8

__all__.append("vDim")
class vDim:
    Scalar  = 0
    Vector  = 1
    Array2D = 2

__all__.append("Fixed")
class Fixed:
    Change = 0
    Fixed  = 1

__all__.append("Null")
class Null:
    NotAllowed = 0
    NotNull    = 1
    IsNull     = 2
    IsPartial  = 3

__all__.append("tMethod")
class tMethod:
    Copy      = 0
    Integrate = 1
    Average   = 2

__all__.append("dMethod")
class dMethod:
    Interpolate = 0
    Conserve    = 1

__all__.append("State")
class State:
    MultipleError = -2
    Open          = -1
    Empty         =  0
    Single        =  1
    Multiple      =  2

# =============================================================================

# Extract OpenCSM error codes
_ocsm_error_codes = {}
with open(os.path.join(_ESP_ROOT,"include","OpenCSM.h")) as fp:
    lines = fp.readlines()
    for line in lines:
        define = line.split()
        if len(define) < 3: continue
        if define[0] != "#define": continue
        try:
            code = int(define[2])
        except:
            continue
        if code > -200: continue
        _ocsm_error_codes[code] = define[1]
del fp
del lines
del define
del code

# =============================================================================

class UnitTupleEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Quantity):
            return self.default( (obj._value, str(obj._units)) )
        return obj

class UnitStringEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, Unit):
            return self.default( str(obj) )
        return obj
# =============================================================================

def signal_handler(sig, frame):
    _caps.caps_rmLock()
    sys.exit(sig)

_caps.caps_externSignal()
signal.signal(signal.SIGSEGV, signal_handler)
signal.signal(signal.SIGINT , signal_handler)
signal.signal(signal.SIGABRT, signal_handler)
if not sys.platform.startswith('win32'):
    signal.signal(signal.SIGHUP , signal_handler)
    signal.signal(signal.SIGBUS , signal_handler)

# =============================================================================

def _decode(data):
    if version_info.major > 2 and isinstance(data, bytes):
        return data.decode()
    return data

# =============================================================================

__all__.append("revision")
def revision():
    """
    Revision - return CAPS revision

    Returns
    -------
        imajor:  major version number
        iminor:  minor version number
    """
    imajor = c_int()
    iminor = c_int()

    stat = _caps.caps_revision(ctypes.byref(imajor), ctypes.byref(iminor))
    if stat: _raiseStatus(stat)

    return (imajor.value, iminor.value)

# =============================================================================

# CAPS error exception class. 
__all__.append("CAPSError")
class CAPSError(Exception):
    
    InternalError = "CAPS Internal Error: "
    
    def __init__(self, status, msg=None, errors=None):
        
        # Combine errors 
        self.errorDict = dict(egads._egads_error_codes)
        self.errorDict.update(_ocsm_error_codes)
        self.errorDict.update(_caps_error_codes)
        
        ## Error code encountered when running the CAPS
        self.status = status
        
        ## Name of error code encountered
        self.errorName = "UNKNOWN_ERROR"
        self.errors = errors

        if status not in self.errorDict:
            msgPre = "CAPS error code " + str(status) + " (undefined error)"
        else:
            msgPre = self.errorDict[status]
            self.errorName = self.errorDict[status]
            
        if msg is not None:
            msg = ": " + msg
        else:
            msg = ""
            
        if errors is not None:
            msg += ":\n"
            msg += "="*80 + "\n"
            for line in errors.info():
                msg += line + "\n"
            msg += "="*80 + "\n"

        super(CAPSError,self).__init__(msgPre + msg)

# =============================================================================

def _raiseStatus(status, msg=None, errors=None):
    # Check status flag returned by CAPS
    
    if status >= CAPS_SUCCESS and status <= CAPS_NOTFOUND: return

    raise CAPSError(status, msg=msg, errors=errors)

# =============================================================================
# __all__.append("c_to_py")
# def c_to_py(c_obj):
#     """
#     Creates a Python class for an exiting c_ego object.
#     
#     Parameters
#     ----------
#     c_obj:
#         the c_capsObj, c_capsOwn, or c_capsErrs instance
#     """
#     if isinstance(c_obj, c_capsObj):
#         return capsObj(c_obj, None, None)
#     elif isinstance(c_obj, c_capsOwn):
#         return capsOwn(c_obj)
#     elif isinstance(c_obj, c_capsErrs):
#         return capsErrs(c_obj)
#     else:
#         raise CAPSError(msg='c_obj must be a c_capsObj, c_capsOwn, or c_capsErrs instance')

#==============================================================================
__all__.append("capsOwn")
class capsOwn:
    """
    Wrapper to represent a capsOwn structure
    """
    def __init__(self, owner, problemObj):
        """
        Constructor for internal use only

        Parameters
        ----------
        owner: 
            c_capsOwn instance

        problemObj:
            capsObj Problem Object
        """
        self._owner = owner
        self._problemObj = problemObj

#==============================================================================
    def info(self):
        """
        Retrieve owner information

        Returns
        -------
        pname:
            the process name
        
        pID:
            the process ID
            
        userID:
            the user ID

        lines:
            Line-by-line list of descriptions
            
        datetime:
            the filled date/time stamp info - 6 in length:
            year, month, day, hour, minute, second
            
        sNum:
            the sequence number (always increasing)
        """
        if self.problemObj._obj is None:
            _raiseStatus(CAPS_CLOSED, "The CAPS Problem object has been closed")

        phase  = c_char_p()
        pname  = c_char_p()
        pID    = c_char_p()
        userID = c_char_p()
        nLines = c_int()
        plines = POINTER(c_char_p)()
        pdatetime = (c_short*6)()
        sNum   = CAPSLONG()
        stat = _caps.caps_ownerInfo(self._problemObj._obj, self._owner, ctypes.byref(phase),
                                    ctypes.byref(pname), ctypes.byref(pID), ctypes.byref(userID), ctypes.byref(nLines), 
                                    ctypes.byref(plines), pdatetime, ctypes.byref(sNum))
        if stat: _raiseStatus(stat)
 
        phase  = _decode(phase.value)
        pname  = _decode(pname.value)
        pID    = _decode(pID.value)
        userID = _decode(userID.value)
        lines  = []
        for i in range(nLines.value):
            lines.append(_decode(plines[i].value))
        datetime = list(pdatetime[:])
        sNum     = sNum.value
 
        return pname, pID, userID, lines, datetime, sNum

# =============================================================================

# Free CAPS errors.
def _caps_freeError(errors):
    
    stat = _caps.caps_freeError(errors)
    if stat: _raiseStatus(stat, msg = "while freeing CAPS capsErrs structure")

#==============================================================================
__all__.append("capsErrs")
class capsErrs:
    """
    Wrapper to represent a list of c_capsErrs structures
    """
    def __init__(self, nErr, errors):
        """
        Constructor for internal use only

        Parameters
        ----------
        nErr: 
            c_int number of errors
            
        errors:
            POINTER(c_capsErrs) instance of errors
        """
        self._finalize = None
        if not isinstance(nErr, c_int): raise CAPSError(msg=CAPSError.InternalError+"nErr must be type c_int")
        if not isinstance(errors, POINTER(c_capsErrs)): raise CAPSError(msg=CAPSError.InternalError+"errors must be type POINTER(c_capsErr)")
        self._nErr = nErr.value
        self._errors = errors
        self._finalize = egads.finalize(self, _caps_freeError, self._errors)

#==============================================================================
    def __del__(self):
        # free the errors
        if self._finalize is not None:
            self._finalize()

#==============================================================================
    def info(self):
        """
        Retrieve error information

        Returns
        -------
            Line-by-line list of error messages 
        """
        lines = []
        for i in range(self._nErr):
            eType = c_int()
            nLines = c_int()
            plines = POINTER(c_char_p)()
            perrObj = POINTER(c_capsObj)()
            stat = _caps.caps_errorInfo(self._errors, c_int(i+1), ctypes.byref(perrObj), ctypes.byref(eType), 
                                        ctypes.byref(nLines), ctypes.byref(plines))
            if stat: _raiseStatus(stat)
            
            for i in range(nLines.value):
                lines.append(_decode(plines[i]))

        return lines

# =============================================================================

# Close a CAPS problem.
def _close_capsProblem(obj, problemState):
    
    stat = _caps.caps_close(obj, 0, None)
    if stat: _raiseStatus(stat, msg = "while closing CAPS Problem")
    problemState.closed = True

#==============================================================================

# free a CAPS Object (only User defined Value and Bound Objects)
def _caps_delete(obj, problemState):
    
    # Check to see if the capsProblem has already been closed and freed this object
    if problemState.closed: return
    
    #stat = _caps.caps_delete(obj)
    #if stat: _raiseStatus(stat, msg = "while deleting CAPS Object")

# =============================================================================
def checkClosed(func):
    """Checks if the problem object has already been closed"""
    @functools.wraps(func)
    def wrapper_checkClosed(*args, **kwargs):
        self = args[0]
        if self.problemObj()._obj is None:
            _raiseStatus(CAPS_CLOSED, "The CAPS Problem object has been closed")
        return func(*args, **kwargs)
    return wrapper_checkClosed

#==============================================================================
__all__.append("phaseState")
def phaseState(prName, phName):
    """
    Check State of CAPS Problem Phase

    Parameters
    ----------
    prName:
        path ending with the CAPS problem name

    phName:
        the current phase name (None is equivalent to 'Scratch')

    Returns
    -------
    bitFlag:
       the returned state (additive): 1 – locked, 2 – closed
    """
    prName = prName.encode() if isinstance(prName, str) else prName
    phName = phName.encode() if isinstance(phName, str) else phName

    bitFlag = c_int()
    stat = _caps.caps_phaseState(prName, phName, ctypes.byref(bitFlag))
    if stat: _raiseStatus(stat)

    return bitFlag.value

#==============================================================================
__all__.append("open")
def open(prName, phName, flag, ptr, outLevel=1):
    """
    Open a CAPS Problem Object

    Parameters
    ----------
    prName:
        path ending with the CAPS problem name
        if exists the stored data initializes the problem, 
        otherwise the directory is created

    phName:
        the current phase name (None is equivalent to 'Scratch')

    ptr:
        path input file name (not needed when restarting) - based on file extension: 
        *.csm initialize the project using the specified OpenCSM file
        *.egads initialize the project based on the static geometry
        - or - 
        pointer to an OpenCSM Model Structure - left open after caps_close

    flag:
        oFlag.oFileName  – ptr is a filename, 
        oFlag.oMODL      – ptr is an OpenCSM Model Structure, 
        oFlag.oEGO       – ptr is a Model ego, 
        oFlag.oPhaseName – ptr is the starting phase name, 
        oFlag.oContinue  – continuation (ptr can be NULL)
        oFlag.oReload    – ptr is the starting phase name with reloading of OpenCSM file

    outLevel:
        0 - minimal, 1 - standard (default), 2 - debug
        
    Returns
    -------
        a capsObj CAPS Problem Object
    """
    prName = prName.encode() if isinstance(prName, str) else prName
    phName = phName.encode() if isinstance(phName, str) else phName
    if isinstance(ptr, (str, bytes)) or \
       (version_info.major <= 2 and isinstance(ptr, unicode)):
        pptr = c_char_p(ptr if isinstance(ptr, bytes) else ptr.encode())
        pptr = ctypes.cast(pptr, c_void_p)
    elif hasattr(ptr, "_modl"):
        flag=oFlag.oMODL
        if not isinstance(ptr._modl, c_void_p):
            raise CAPSError(CAPS_BADVALUE, "ptr has _modl that is not a ctypes.c_void_p. Is it an intance of pyOCSM.Ocsm? ptr = {!r}".format(ptr))
        pptr = ctypes.cast(ptr._modl, c_void_p)
    else:
        raise CAPSError(CAPS_BADVALUE, "ptr type should be str, bytes, or pyOCSM.Ocsm: ptr = {!r}".format(ptr))

    nErr = c_int()
    errs = POINTER(c_capsErrs)()
    problemObj = c_capsObj()
    stat = _caps.caps_open(prName, phName, c_int(flag), pptr, c_int(outLevel), ctypes.byref(problemObj), 
                           ctypes.byref(nErr), ctypes.byref(errs))
    if stat:
        msg = "Failed to open Problem with prName={!r}, phName={!r}, flag={!r}".format(_decode(prName),_decode(phName),flag)
        if flag != oFlag.oMODL:
            msg += ", and ptr={!r}".format(ptr)
        _raiseStatus(stat, msg=msg, errors=capsErrs(nErr, errs))

    return capsObj(problemObj, None, deleteFunction=_close_capsProblem)

#==============================================================================
class capsObj:
    """
    Base class Wrapper to represent a CAPS Object
    """
#==============================================================================
    def __init__(self, obj, problemObj, deleteFunction):
        """
        Constructor only intended for internal use.

        Parameters
        ----------
        obj:
            c_capsObj Object

        problemObj:
            capsObj Problem Object (None if c_capsObj is a Problem Object)

        deleteFunction:
            if not None. the class instance will call deleteFunction
            for the provided c_capsObj during garbage collection
        """
        self._finalize = None
        if not isinstance(obj, c_capsObj): raise CAPSError(msg=CAPSError.InternalError+"obj must be type c_capsObj")
        self._obj = obj
        self._problemObj = problemObj
                
        if problemObj is not None:
            problemState = problemObj._problemState
        else:
            class ProblemState:
                pass
            self._problemState = ProblemState()
            self._problemState.closed = False
            problemState = self._problemState

        if deleteFunction is not None:
            self._finalize = egads.finalize(self, deleteFunction, obj, problemState)

#==============================================================================
    def __del__(self):
        # free the object
        if self._finalize is not None:
            self._finalize()
        self._obj = None
        self._finalize = None
        self._problemObj = None

#==============================================================================
    @checkClosed
    def close(self, complete = 0, phName = None):
        """
        Closes a Problem object
        
        Parameters
        ----------
        self:
           a CAPS Problem Object

        complete:
            0 – the phase is not complete, 1 – the phase is completed and should not be modified

        phName:
            Phase Name of the Scratch phase is closed as complete
        """
        # detach the finalizer
        if self._finalize is not None:
            self._finalize.detach()

        phName   = phName.encode() if isinstance(phName, str) else phName

        stat = _caps.caps_close(self._obj, c_int(complete), phName)
        if stat: _raiseStatus(stat, msg = "while closing CAPS Problem")

        self._obj = None
        self._finalize = None
        self._problemObj = None
        self._problemState.closed = True

#==============================================================================
    @checkClosed
    def __eq__(self, obj):
        """
        checks pointer equality of underlying c_capsObj memory address
        """
        return isinstance(obj, capsObj) and ctypes.addressof(obj._obj.contents) == ctypes.addressof(self._obj.contents)

#==============================================================================
    @checkClosed
    def __ne__(self, obj):
        return not self == obj

#==============================================================================
    def problemObj(self):
        """
        Returns the problemObj associated with this capsObj (self of self is a problemObj)
        """
        return self if self._problemObj is None else self._problemObj

#==============================================================================
    @checkClosed
    def info(self):
        """
        Return information about the object

        Returns
        -------
        name:
            the returned Object name pointer (if any)

        otype:
           the returned Object type: Problem, Value, Analysis, Bound, VertexSet, DataSet

        stype:
            the returned subtype (depending on type)

        link:
            the returned linkage Value Object (None - no link)

        parent:
            the returned parent Object (None for a Problem or an Attribute generated User Value)

        last:
            the returned last owner/history to touch the Object
        """
        name   = c_char_p()
        otype  = c_int()
        stype  = c_int()
        link   = c_capsObj()
        parent = c_capsObj()
        last   = c_capsOwn()
        stat = _caps.caps_info(self._obj, ctypes.byref(name), ctypes.byref(otype),
                               ctypes.byref(stype), ctypes.byref(link), ctypes.byref(parent), ctypes.byref(last))
        if stat: _raiseStatus(stat)
        
        name   = _decode(name.value)
        otype  = otype.value
        stype  = stype.value
        link   = capsObj(link  , self.problemObj(), None) if link   else None
        parent = capsObj(parent, self.problemObj(), None) if parent else None
        last   = capsOwn(last, self.problemObj())         if last   else None
        
        return name, otype, stype, link, parent, last

#==============================================================================
    @checkClosed
    def size(self, type, stype):
        """
        Children sizing information from a Patent Object

        Parameters
        ----------
        type:
           the data type to size: Bodies, Problem, Value, Analysis, Bound, VertexSet, DataSet

        stype:
            the subtype to size (depending on type)

        Returns
        -------
        the size
        """
        size = c_int()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_size(self._obj, c_int(type), c_int(stype), ctypes.byref(size), 
                               ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        return size.value

#==============================================================================
    @checkClosed
    def childByName(self, otype, stype, name):
        """
        Retrieve a CAPS child object by name

        Parameters
        ----------
        otype:
            the Object type to return: Value, Analysis, Bound, VertexSet, DataSet

        stype:
            the subtype to find (depending on type)

        name:
            string name of the object

        Returns
        -------
        child:
            the returned CAPS Object
        """
        name   = name.encode() if isinstance(name,str) else name
        child  = c_capsObj()
        units  = c_char_p()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_childByName(self._obj, c_int(otype), c_int(stype), 
                                      c_char_p(name), ctypes.byref(child),
                                      ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        return capsObj(child, self.problemObj(), None)

#==============================================================================
    @checkClosed
    def childByIndex(self, otype, stype, index):
        """
        Retrieve a CAPS child object by name

        Parameters
        ----------
        otype:
            the Object type to return: Value, Analysis, Bound, VertexSet, DataSet

        stype:
            the subtype to find (depending on type)

        index:
            the index [1 - size]

        Returns
        -------
        child:
            the returned CAPS Object
        """
        child  = c_capsObj()
        units  = c_char_p()
        stat = _caps.caps_childByIndex(self._obj, c_int(otype), c_int(stype), 
                                       c_int(index), ctypes.byref(child))
        if stat: _raiseStatus(stat)
        return capsObj(child, self.problemObj(), None)

#==============================================================================
    # @checkClosed
    # def getHistory(self):
    #     """
    #     Retrieve history of a CAPS Object
    #
    #     Returns
    #     -------
    #         list of capsOwn structures
    #     """
    #     nhist = c_int()
    #     phist = POINTER(c_capsOwn)()
    #     stat = _caps.caps_getHistory(self._obj, ctypes.byref(nhist), ctypes.byref(phist))
    #     if stat: _raiseStatus(stat)
    #
    #     nhist = nhist.value
    #     hist = [capsOwn(phist[i]) for i in range(nhist)]
    #
    #     return hist

#==============================================================================
    @checkClosed
    def attrByName(self, name):
        """
        Retrieve an attribute by name

        Parameters
        ----------
        name:
            string name of the Attribute
            
        Returns
        -------
            User Value Object
        """
        name = name.encode() if isinstance(name,str) else name
        
        attr = c_capsObj()
        
        stat = _caps.caps_attrByName(self._obj, name, ctypes.byref(attr))
        if stat: _raiseStatus(stat)

        deleteFunc = _caps_delete if attr.contents.subtype == sType.USER else None

        return capsObj(attr, self.problemObj(), deleteFunction=deleteFunc)

#==============================================================================
    @checkClosed
    def attrByIndex(self, index):
        """
        Retrieve an Attribute by index

        Parameters
        ----------
        index:
            the index (bias 1) to the list of Attributes
            
        Returns
        -------
            User Value Object
        """
        attr = c_capsObj()
        
        stat = _caps.caps_attrByIndex(self._obj, c_int(index), ctypes.byref(attr))
        if stat: _raiseStatus(stat)
        
        deleteFunc = _caps_delete if attr.contents.subtype == sType.USER else None

        return capsObj(attr, self.problemObj(), deleteFunction=deleteFunc)

#==============================================================================
    @checkClosed
    def setAttr(self, attr, name=None):
        """
        Set an Attribute

        Parameters
        ----------
        attr:
            the Value Object containing the attribute
            The attribute will not maintain the Value Object's shape
            
        name:
            a string referring to the Attribute name - None: use name in attr
            Note: an existing Attribute of this name is overwritten with the new value
        """
        name = name.encode() if isinstance(name,str) else name

        stat = _caps.caps_setAttr(self._obj, name, attr._obj)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def deleteAttr(self, name):
        """
        Set an Attribute

        Parameters
        ----------
        name:
            a string referring to the Attribute to delete
            None deletes all attributes attached to the Object
        """
        name = name.encode() if isinstance(name,str) else name
       
        stat = _caps.caps_deleteAttr(self._obj, name)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def journalState(self):
        """
        Get the journal state for a problem object
        
        Note: Needed to determine if restarting when explicitly executing analyses.

        Parameters
        ----------
        self:
           a CAPS Problem Object

        Returns
        -------
        True if Journaling False otherwise
        """
        stat = _caps.caps_journalState(self._obj)
        if stat == 0: return False
        if stat == 4: return True
        _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def setOutLevel(self, outlevel=1):
        """
        Sets the CAPS verbosity level

        Parameters
        ----------
        outlevel: 0 <= outlevel <= 2
            0 - minimal, 1 - standard (default), 2 - debug

        Returns
        -------
        The old outLevel
        """
        stat = _caps.caps_outLevel(self._obj, c_int(outlevel))
        if stat < 0: _raiseStatus(stat)
        return stat

#==============================================================================
    @checkClosed
    def getRootPath(self):
        """
        Get Problem root

        Parameters
        ----------
        self: 
            a CAPS Problem Object

        Returns
        -------
            the file path to find the root of the Problem/Phase directory structure 
            if on Windows it will contain the drive
        """
        root  = c_char_p()
        stat = _caps.caps_getRootPath(self._obj, ctypes.byref(root))
        if stat: _raiseStatus(stat)
        return _decode(root.value)

#==============================================================================
    # @checkClosed
    def intentPhrase(self, lines):
        """
        Get Problem root
    
        Parameters
        ----------
        self: 
            a CAPS Problem Object
    
        lines:
            List of intent phrases for history tracking
        """
        nlines = len(lines)
        plines = (c_char_p*nlines)()
        for i in range(nlines):
            plines[i] = lines[i].encode() if isinstance(lines[i], str) else lines[i]
    
        stat = _caps.caps_intentPhrase(self._obj, c_int(nlines), plines)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def bodyByIndex(self, index):
        """
        Retrived an EGADS body ego based on index

        Parameters
        ----------
        index:
            the index [1-size] - see caps_size

        Returns
        -------
        body:
            the returned EGADS Body Object
        
        units:
            the Unit declaring the length units - None for unitless values
        """
        body   = egads.c_ego()
        units  = c_char_p()
        stat = _caps.caps_bodyByIndex(self._obj, c_int(index), 
                                      ctypes.byref(body), ctypes.byref(units))
        if stat: _raiseStatus(stat)
        if units: units = Unit(_decode(units.value))
        return egads.c_to_py(body), units

#==============================================================================
    @checkClosed
    def getBodies(self):
        """
        Retrived EGADS ego bodies from an Analysis Object

        Returns
        -------
        bodies:
            the returned list of EGADS ego Body Objects
        """
        nbody  = c_int()
        pbodies = POINTER(egads.c_ego)()
        units  = c_char_p()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_getBodies(self._obj, ctypes.byref(nbody), ctypes.byref(pbodies), ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        bodies = [egads.ego(egads.c_ego(pbodies[i].contents)) for i in range(nbody.value)]
        
        return bodies

#==============================================================================
    @checkClosed
    def convertValue(self, inVal, inUnit):
        """
        Perform unit conversion

        Parameters
        ----------
        inVal:
            the input value to convert
            
        inUint:
            the units of the input value (may be None if self is a Value Object)
        
        Returns
        -------
            inVal in the units of the Value Object units
        """
        inUnit = inUnit.encode() if isinstance(inUnit, str) else inUnit

        outVal = c_double()
        stat = _caps.caps_convertValue(self._obj, c_double(inVal), inUnit, ctypes.byref(outVal))
        if stat: _raiseStatus(stat)
        
        return outVal.value

#==============================================================================
    @checkClosed
    def reset(self):
        """
        Reset the problem. This is equivalent to a clean slate restart.
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_reset(self._obj, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def writeParameters(self, fileName):
        """
        This outputs an OpenCSM Design Parameter file.
        
        Parameters
        ----------
        fileName:
            the name of the parameter file to write
        """
        fileName = fileName.encode() if isinstance(fileName,str) else fileName

        stat = _caps.caps_writeParameters(self._obj, fileName)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def readParameters(self, fileName):
        """
        This reads an OpenCSM Design Parameter file and overwrites (makes dirty) 
        the current state for the GeometryIn Values in the file.
        
        Parameters
        ----------
        fileName:
            the name of the parameter file to read
        """
        fileName = fileName.encode() if isinstance(fileName,str) else fileName

        stat = _caps.caps_readParameters(self._obj, fileName)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def writeGeometry(self, fileName, flag = 1):
        """
        Writes all bodies in the Problem/Analysis Object to a file
        
        Parameters
        ----------
        fileName:
            the name of the file to write - typed by extension (case insensitive): 
            iges/igs - IGES File
            step/stp - STEP File
            brep - OpenCASCADE File
            egads - EGADS file (which includes attribution)
            
        flag:
            the write flag: 0 -- no additional output,  1 -- also write Tessellation Objects for
            EGADS output (only for Analysis Objects)
        """
        fileName = fileName.encode() if isinstance(fileName,str) else fileName

        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_writeGeometry(self._obj, c_int(flag), fileName, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def makeValue(self, vname, stype, data):
        """
        Create a value object
        
        Parameters
        ----------
        vname:
            the Value Object name to be created

        stype:
            the Object subtype: Parameter or User
             
        data:
            the data (may be Unit instances) to store in the Value Object
            if data is string and units is "PATH" - slashes are converted automatically

        Returns
        -------
            the CAPS Value Object
        """
        if isinstance(data, Quantity):
            units = str(data._units)
            data  =     data._value
        else:
            units = None

        vname = vname.encode() if isinstance(vname,str) else vname
        units = units.encode() if isinstance(units,str) else units

        vtype, nrow, ncol, pdata, partial = self._convertData(data)
 
        # Convert integers with units to doubles with units
        if vtype == vType.Integer and units is not None:
            tmp = (c_double*(nrow*ncol))()
            for i in range(nrow*ncol):
                tmp[i] = pdata[i]
            pdata = tmp
            vtype = vType.Double
 
        val = c_capsObj()
        stat = _caps.caps_makeValue(self._obj, vname, c_int(stype), c_int(vtype), 
                                    c_int(nrow), c_int(ncol), pdata, partial, units, ctypes.byref(val))
        if stat: _raiseStatus(stat)

        deleteFunc = _caps_delete if stype == sType.USER else None

        return capsObj(val, self.problemObj(), deleteFunction=deleteFunc)

#==============================================================================
    @checkClosed
    def queryAnalysis(self, aname):
        """
        Query Analysis -- Does not 'load' or create an object
        
        Note: this causes the the DLL/Shared-Object to be loaded (if not already resident)

        Parameters
        ----------
        aname:
            the Value Object name to be created

        Returns
        -------
        nIn:
            the number of inputs
            
        nOut:
            the number of outputs
             
        execute:
            execution flag: 0 - no exec, 1 - AIM Execute exists
        """
        aname = aname.encode() if isinstance(aname,str) else aname
 
        nIn = c_int()
        nOut = c_int()
        execute = c_int()
        stat = _caps.caps_queryAnalysis(self._obj, aname, 
                                        ctypes.byref(nIn), ctypes.byref(nOut), 
                                        ctypes.byref(execute))
        if stat: _raiseStatus(stat)
        
        nIn = nIn.value
        nOut = nOut.value
        execute = execute.value

        return nIn, nOut, execute

#==============================================================================
    @checkClosed
    def getInput(self):
        raise CAPSError(msg="Implement")

#==============================================================================
    @checkClosed
    def getOutput(self):
        raise CAPSError(msg="Implement")

#==============================================================================
    @checkClosed
    def makeAnalysis(self, aname, name, unitSys=None, intent=None, execute=1):
        """
        Load Analysis into the Problem

        Notes: 
            This causes the the DLL/Shared-Object to be loaded (if not already resident)

            If execute is 1 and the AIM has aimExecute, aimExecute automatically runs after preAnalysis
            and if the execution is not asynchronous aimPostAnalysis is automatically run. 

        Parameters
        ----------
        aname:
            the Analysis (and AIM plugin) name

        name:
            the unique supplied name for this instance (can be  None)

        unitSys: 
            string describing the unit system to be used by the AIM (can be None)
            see specific AIM documentation for a list of strings for which the AIM will respond

        intent: 
            the intent string used to pass Bodies to the AIM, None - no filtering

        execute: 
            the execution flag:  0 - No auto execution, 1 - Auto execute if available

        Returns
        -------
            the Analysis Object
        """
        aname   = aname.encode()   if isinstance(aname  ,str) else aname
        name    = name.encode()    if isinstance(name   ,str) else name
        unitSys = unitSys.encode() if isinstance(unitSys,str) else unitSys
        
        if unitSys:
            if isinstance(unitSys,dict):
                unitSys = json.dumps(unitSys, cls=UnitStringEncoder).encode()
            else:
                _raiseStatus(CAPS_UNITERR, "unitSys must be a dictionary: unitSys = {!r}".format(unitSys))
                
        if isinstance(intent, (str, bytes)):
            intent = intent.encode() if isinstance(intent, str) else intent
        elif isinstance(intent, list):
            intent = b";".join( [i.encode() if isinstance(i ,str) else i for i in intent] )

        if intent == b"": intent = None 
        analysis = c_capsObj()
        execute = c_int(execute)
        nErr = c_int()
        errs = POINTER(c_capsErrs)()

        stat = _caps.caps_makeAnalysis(self._obj, aname, name, unitSys, intent,
                                       ctypes.byref(execute), ctypes.byref(analysis), 
                                       ctypes.byref(nErr), ctypes.byref(errs))
        if stat:
            msg = "Failed to make AIM '" + _decode(aname) + "'"
            if name is not None:
                msg += " with name '" + _decode(name) + "'"
            _raiseStatus(stat, msg=msg, errors=capsErrs(nErr, errs))

        return capsObj(analysis, self.problemObj(), None)

#==============================================================================
    @checkClosed
    def makeBound(self, dim, bname):
        """
        Creates an Bound Object

        Parameters
        ----------
        dim:
              the dimensionality of the Bound (1 - 3)
              
        bname:
             string associated with "capsBound" attribute on bodies

        Returns
        -------
            the Bound Object
        """
        bname = bname.encode() if isinstance(bname ,str) else bname

        bound = c_capsObj()
        stat = _caps.caps_makeBound(self._obj, c_int(dim), bname, ctypes.byref(bound))
        if stat: _raiseStatus(stat)

        return capsObj(bound, self.problemObj(), deleteFunction=_caps_delete)

#==============================================================================
    @checkClosed
    def getValue(self):
        """
        Returns the data in the Value Object
        
        Returns
        ----------
        data:
            the data stored in the Value Object
        """
        vtype   = c_int()
        nrow    = c_int()
        ncol    = c_int()
        partial = POINTER(c_int)()
        pdata   = c_void_p()
        units   = c_char_p()
        nErr    = c_int()
        errs    = POINTER(c_capsErrs)()

        stat = _caps.caps_getValue(self._obj, ctypes.byref(vtype), 
                                   ctypes.byref(nrow), ctypes.byref(ncol), 
                                   ctypes.byref(pdata), ctypes.byref(partial), 
                                   ctypes.byref(units), ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

        vtype = vtype.value
        nrow = nrow.value
        ncol = ncol.value
        units = _decode(units.value) if units else None

        if not pdata:
            return None

        # Transpose for column vector
        if ncol == 1 and nrow > 1:
            ncol = nrow
            nrow = 1

        if vtype == vType.String:
            size = ctypes.sizeof(ctypes.c_char)
            slen = 0

            data = [None]*nrow
            for i in range(nrow):
                data[i] = [None]*ncol
                for j in range(ncol):
                    chars = ctypes.cast(c_void_p(pdata.value + slen*size), c_char_p)
                    slen += len(chars.value)+1
                    data[i][j] = _decode(chars.value)

        if vtype == vType.Boolean or \
           vtype == vType.Integer or \
           vtype == vType.Double  or \
           vtype == vType.DoubleDeriv:
            
            if vtype == vType.Double or vtype == vType.DoubleDeriv:
                pdata = ctypes.cast(pdata, POINTER(c_double))
            else:
                pdata = ctypes.cast(pdata, POINTER(c_int))
                
            data = [[pdata[i*ncol+j] for j in range(ncol)] for i in range(nrow)]

        elif vtype == vType.Tuple:
            pdata = ctypes.cast(pdata, POINTER(c_capsTuple))
            data = {}
            for i in range(nrow):
                for j in range(ncol):
                    
                    name = _decode(pdata[i*ncol+j].name)
                    valu = _decode(pdata[i*ncol+j].value)
                    try:
                        valu = json.loads( valu )
                    except:
                        pass

                    if partial and partial[i*ncol+j] == Null.IsNull:
                        data[name] = None
                    else:
                        data[name] = valu

            return data

        if partial:
            data = [[data[i][j] if partial[i*ncol+j] == Null.NotNull else None for j in range(ncol)] for i in range(nrow)]

        if nrow == 1 and ncol == 1:
            data = data[0][0]
        elif nrow == 1 or ncol == 1:
            data = data[0]

        return _withUnits(data, units)

#==============================================================================
    @checkClosed
    def getValueUnit(self):
        """
        Returns the unit in the Value Object
        
        Returns
        ----------
        units:
            Unit instance of the units or None
        """
        vtype   = c_int()
        nrow    = c_int()
        ncol    = c_int()
        partial = POINTER(c_int)()
        pdata   = c_void_p()
        units   = c_char_p()
        nErr    = c_int()
        errs    = POINTER(c_capsErrs)()

        stat = _caps.caps_getValue(self._obj, ctypes.byref(vtype), 
                                   ctypes.byref(nrow), ctypes.byref(ncol), 
                                   ctypes.byref(pdata), ctypes.byref(partial), 
                                   ctypes.byref(units), ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        units = _decode(units.value) if units else None
        
        if units is None:
            return None
        
        return Unit(units)

#==============================================================================
    # Determine the equivalent CAPS value object type for a Python object     
    def _getType(self, data, level = 0):
    
        # If we have a list loop back in to see what type of value the element is 
        if isinstance(data, list) or \
           isinstance(data, tuple):
            vtype = self._getType(data[0], level+1)
            for i in range(1,len(data)):
                # Allow convertsion from Integer to Double
                vtypei = self._getType(data[i])
                if vtype == vType.Null:
                    vtype = vtypei
                if vtypei == vType.Null:
                    vtypei = vtype
                if vtype == vType.Integer and vtypei == vType.Double:
                    vtype = vType.Double
                if vtypei == vType.Integer and vtype == vType.Double:
                    vtypei = vType.Double
                if vtype != vtypei:
                    raise CAPSError(CAPS_BADVALUE, "List entries must all be same type! data: {!r}".format(data))
            return vtype

        if isinstance(data, bool) and (str(data) == 'True' or str(data) == 'False'):
            return vType.Boolean
        
        if isinstance(data, int):
            return vType.Integer

        if numpy and (isinstance(data, numpy.int32) or isinstance(data, numpy.int64)):
            return vType.Integer

        if isinstance(data, float):
            return vType.Double

        if numpy and (isinstance(data, numpy.float32) or isinstance(data, numpy.float64)):
            return vType.Double

        if numpy and isinstance(data, numpy.ndarray):
            return self._getType(data[0], level+1)

        if isinstance(data, str):
            return vType.String

        # dictionaries are converted to json strings
        if isinstance(data, dict):
            if level == 0:
                return vType.Tuple
            else:
                return vType.JsonDict

        if data is None:
            return vType.Null

        raise CAPSError(CAPS_BADVALUE, "Unsupported data type! type: {!r}, data: {!r}".format(type(data), data))

#==============================================================================
    def _getShape(self, vtype, data):
        
        nrow = len(data)
        ncol = 1
        
        if isinstance(data, dict):
            return nrow, ncol
        elif isinstance(data[0], list) or \
           isinstance(data[0], tuple):
            ncol = len(data[0])
        elif numpy and isinstance(data, numpy.ndarray):
            if len(data.shape) == 1:
                return data.shape[0], 1
            else:
                return data.shape
        else:
            ncol = 1
            
        # Make sure shape is consistent
        if ncol != 1:
            for i in range(nrow):
                if len(data[i]) != ncol:
                    raise CAPSError(CAPS_BADVALUE, "Inconsistent column sizes!")
        
        return nrow, ncol

#==============================================================================
    def _convertData(self, data):

        vtype = self._getType(data)
        
        # Convert data to list 
        if  not (isinstance(data, list) or 
                 isinstance(data, tuple) or 
                (numpy and isinstance(data, numpy.ndarray)) or 
                 isinstance(data, dict)):
            data = [data] # Convert to list

        nrow, ncol = self._getShape(vtype, data)

        length = nrow*ncol
        partial = None

        # Boolean, Integer, or Double
        if vtype == vType.Boolean or vtype == vType.Integer or vtype == vType.Double:
            
            partial = (enum_capsNull*(length))()
            for i in range(length):
                partial[i] = Null.NotNull
            isPartial = False

            if vtype == vType.Double:
                val = (c_double*length)()
            else:
                val = (c_int*length)()
        
            for i in range(nrow):
                for j in range(ncol):
                    if not hasattr(data[i], "__len__"):
                        if data[i] is None:
                            partial[i*ncol+j] = Null.IsNull
                            val[i*ncol+j] = 0.0
                            isPartial = True
                            continue
                        else:
                            val[i*ncol+j] = data[i]
                    else:
                        if data[i][j] is None:
                            partial[i*ncol+j] = Null.IsNull
                            val[i*ncol+j] = 0.0
                            isPartial = True
                            continue
                        else:
                            val[i*ncol+j] = data[i][j]

            if not isPartial: partial = None

        # String
        elif vtype == vType.String:
            if ncol == 1:
                val = b"\0".join([data[i].encode() if isinstance(data[i], str) else data[i] for i in range(nrow)])
            else:
                val = b"\0".join([data[i][j].encode() if isinstance(data[i][j], str) else data[i][j] for j in range(ncol) for i in range(nrow)])

        # Tuple e.g. dict
        elif vtype == vType.Tuple:

            val = (c_capsTuple*length)()

            for i, key in enumerate(sorted(data.keys())):
               
                temp = json.dumps(data[key], cls=UnitTupleEncoder)
                
                byteName   = key.encode() if isinstance(key, str) else key
                byteString = temp.encode() if isinstance(temp, str) else temp

                val[i].name  = byteName
                val[i].value = byteString
                
        # Dict converted to Json string
        elif vtype == vType.JsonDict:
            val = ';'.join( [ json.dumps(i) for i in data ] )
            val = val.encode() if isinstance(val, str) else val
            vtype = vType.String

        # Null (None)
        elif vtype == vType.Null:
            val = None

        # Unknown type
        else:
            raise CAPSError(CAPS_BADVALUE, "Can not convert Python object to type: " + str(vtype))

        return vtype, nrow, ncol, val, partial

#==============================================================================
    @checkClosed
    def setValue(self, data):
        """
        Set data in the Value object
        
        Parameters
        ----------
        data:
            data to store in the Value Object
            
        units:
            pointer to the string declaring the units - None for unitless values
            if data is string and units is "PATH" - slashes are converted automatically
        """
        if isinstance(data, Quantity):
            units = str(data._units)
            data  =     data._value
        else:
            units = None

        units = units.encode() if isinstance(units,str) else units
        if units is not None and not isinstance(units, bytes):
            raise CAPSError(CAPS_UNITERR, "units must be a string: units = {!r}".format(units))

        vtype, nrow, ncol, pdata, partial = self._convertData(data)

        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_setValue(self._obj, c_int(vtype), 
                                   c_int(nrow), c_int(ncol), 
                                   pdata, partial, units, ctypes.byref(nErr), ctypes.byref(errs))
        if stat:
            name, otype, stype, link, parent, last = self.info()
            _raiseStatus(stat, msg="Trying to set value for: {!r}".format(name), errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def getLimits(self):
        """
        Returns the value limits
        
        Returns
        -------
        limits:
            list of length 2 or None if not set 

        units:
            units of the limits or None
        """
        vtype = c_int()
        plimits = c_void_p()
        punits = c_char_p()

        stat = _caps.caps_getLimits(self._obj, ctypes.byref(vtype), ctypes.byref(plimits), ctypes.byref(punits))
        if stat: _raiseStatus(stat)
        
        if not plimits:
            return None
        
        vtype = vtype.value
        units = _decode(punits.value) if punits else None

        if vtype == vType.Integer:
            limits = list(ctypes.cast(plimits, POINTER(c_int))[0:2])
        elif vtype == vType.Double:
            limits = list(ctypes.cast(plimits, POINTER(c_double))[0:2])
        else:
            raise CAPSError(CAPS_BADVALUE, CAPSError.InternalError+"Unknown vtype: " + str(vtype))

        return _withUnits(limits, units)

#==============================================================================
    @checkClosed
    def setLimits(self, limits):
        """
        Set the value limits
        
        Paramters
        -------
        limits:
            list (or Unit) of length 2 of int or float
        """
        if isinstance(limits, Quantity):
            units  = str(limits._units)
            limits =     limits._value
        else:
            units = None
        
        if limits is None:
            vtype = 0
            plimits = None
        else:
            if not (isinstance(limits, list) or \
                    isinstance(limits, tuple)) or \
               len(limits) != 2:
                raise CAPSError(CAPS_BADVALUE, "limits should be 2 element list - [min value, max value]!")
            
            for i in limits:
                if not isinstance(i, (int, float)):
                    raise CAPSError(CAPS_BADVALUE, "Invalid element type for limits value, only int or float values are valid!")

            vtype = self._getType(limits)
           
            if vtype == vType.Integer:
                plimits = (c_int*2)()
            elif vtype == vType.Double:
                plimits = (c_double*2)()
            else:
                raise CAPSError(CAPS_BADVALUE, CAPSError.InternalError+"Unknown vtype: " + str(vtype))

            plimits[:2] = limits[:]
            plimits = ctypes.cast(plimits, ctypes.c_void_p)
            
        units = units.encode() if isinstance(units, str) else units
        if units is not None and not isinstance(units, bytes):
            raise CAPSError(CAPS_UNITERR, "units must be a string: units = {!r}".format(units))

        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_setLimits(self._obj, c_int(vtype), plimits, units, 
                                    ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def getValueProps(self):
        """
        Returns the value properties
        
        Returns
        -------
        dim:
            the returned dimensionality:
            0 - scalar only, 1 - vector or scalar, 2 - scalar, vector, or 2D array

        pmtr:
            the returned flag: 0 - normal, 1 - GeometryIn type -> OCSM_CFGPMTR
            
        lfix:
            0 - the length(s) can change, 1 -- the length is fixed
            
        sfix:
            0 - the Shape can change, 1 - Shape is fixed

        ntype:
            0 - None invalid, 1 - not None, 2 - is None, 3 - partial None
        """
        dim   = c_int()
        pmtr  = c_int()
        lfix  = c_int()
        sfix  = c_int()
        ntype = c_int()

        stat = _caps.caps_getValueProps(self._obj, ctypes.byref(dim), ctypes.byref(pmtr), 
                                        ctypes.byref(lfix), ctypes.byref(sfix), ctypes.byref(ntype))
        if stat: _raiseStatus(stat)
        
        dim   = dim.value
        pmtr  = pmtr.value
        lfix  = lfix.value
        sfix  = sfix.value
        ntype = ntype.value

        return dim, pmtr, lfix, sfix, ntype

#==============================================================================
    @checkClosed
    def setValueProps(self, dim, lfix, sfix, ntype):
        """
        Sets the value properties (only for the User & Parameter subtypes) 
        
        Parameters
        -------
        dim:
            the dimensionality:
            0 - scalar only, 1 - vector or scalar, 2 - scalar, vector, or 2D array
            
        lfix:
            0 - the length(s) can change, 1 -- the length is fixed
            
        sfix:
            0 - the Shape can change, 1 - Shape is fixed

        ntype:
            0 - None invalid, 1 - not None, 2 - is None, 3 - partial None
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_setValueProps(self._obj, c_int(dim),
                                        c_int(lfix), c_int(sfix), c_int(ntype),
                                        ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def transferValues(self, tmethod, src):
        """
        Transfers values from src Value Object (not for AIM pointer or Tuple vtypes) - or - DataSet Object
        to self Value Object
        
        Parameters
        ----------
        self:
            the self Value Object to receive the data
            Notes: 
                Must not be GeometryOut or AnalysisOut
                Shapes must be compatible
                Overwrites any Linkage

        tmethod:
            0 - copy, 1 - integrate, 2 - weighted average  -- (1 & 2 only for DataSet src)
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()

        stat = _caps.caps_transferValues(src._obj, c_int(tmethod), self._obj,
                                         ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def linkValue(self, link, tmethod):
        """
        Links values from Value Object (not for Tuple vtype or Value subtype User)  - or - DataSet Object
        to a target Value Object
        
        self is the target Value Object which will recived the data
        Notes: 
                Must not be GeometryOut or AnalysisOut
                (Sub)shapes must be compatible

        Note: circular linkages are not allowed!
        
        Parameters
        ----------
        link:
            the source Value or DataSet Object

        tmethod:
            0 - copy, 1 - integrate, 2 - weighted average  -- (1 & 2 only for DataSet self)
        """
        if not isinstance(link, capsObj):
            raise CAPSError(CAPS_BADVALUE, "link must be of type capsObj")

        if link._obj.contents.type != oType.VALUE:
            raise CAPSError(CAPS_BADVALUE, "link must be a Value Object")

        if self._obj.contents.type != oType.VALUE:
            raise CAPSError(CAPS_BADVALUE, "self must be a Value Object")

        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_linkValue(link._obj, c_int(tmethod), self._obj, 
                                    ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def removeLink(self):
        """
        Remove a link to this Value Object
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_linkValue(None, c_int(0), self._obj, 
                                    ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def hasDeriv(self):
        """
        Get a list of the Derivatives available
        
        Only for GeometryOut/AnalysisOut or Double /w Dot Value Objects
        
        Returns
        -------
            the list of derivative names
            derivatives derived from vectors/arrays will have "[n]" or "[n,m]" appended
        """
        
        nderiv = c_int()
        pnames = POINTER(c_char_p)()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()

        stat = _caps.caps_hasDeriv(self._obj, ctypes.byref(nderiv), ctypes.byref(pnames),
                                   ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        if not pnames:
            return None
        
        nderiv = nderiv.value
        names = [None]*nderiv
        for i in range(nderiv):
            names[i] = _decode(pnames[i])
        
        egads.free(pnames)
        
        return names

#==============================================================================
    @checkClosed
    def getDeriv(self, name):
        """
        Get Derivative values
        
        Only for GeometryOut/AnalysisOut or Double /w Dot Value Objects
        
        Parameters
        ----------
        name:
            the name of the derivative

        Returns
        -------
            list of values (tuples if rank > 1) of sensitvities
        """
        
        name    = name.encode() if isinstance(name,str) else name
        len     = c_int()
        len_wrt = c_int()
        pderiv  = POINTER(c_double)()
        nErr    = c_int()
        errs    = POINTER(c_capsErrs)()

        stat = _caps.caps_getDeriv(self._obj, name, ctypes.byref(len), ctypes.byref(len_wrt), ctypes.byref(pderiv),
                                   ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        len     = len.value
        len_wrt = len_wrt.value
        
        if len == 0:
            return None
        
        if len == 1 and len_wrt == 1:
            deriv = pderiv[0]
        elif len == 1 or len_wrt == 1:
            deriv = list(pderiv[0:len*len_wrt])
        else:
            deriv = [None]*len
            for i in range(len):
                deriv[i] = list(pderiv[i*len_wrt:(i+1)*len_wrt])
        
        return deriv

#==============================================================================
    @checkClosed
    def dupAnalysis(self, name):
        """
        Initializeds a new analysis based on this Analysis Object
        
        Parameters
        ----------
        self:
            CAPS Analysis Object
            
        name:
            the name of the duplicated analysis object
            
        Returns
        -------
            a copy of the CAPS Analsysis Object
        """
        analysis = c_capsObj()
        name = name.encode() if isinstance(name,str) else name
        stat = _caps.caps_dupAnalysis(self._obj, name, ctypes.byref(analysis))
        if stat: _raiseStatus(stat)

        return capsObj(analysis, self.problemObj(), None)

#==============================================================================
    @checkClosed
    def resetAnalysis(self):
        """
        Reset the Analysis
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_resetAnalysis(self._obj, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def dirtyAnalysis(self):
        """
        Resturns the list of dirty Analysis Objects
        """
        nAobj = c_int()
        aobjs = POINTER(c_capsObj)()
        stat = _caps.caps_dirtyAnalysis(self._obj, ctypes.byref(nAobj), ctypes.byref(aobjs))
        if stat: _raiseStatus(stat)

        anlysis = [capsObj(c_capsObj(aobjs[i].contents), self.problemObj(), None) for i in range(nAobj.value)]
 
        egads.free(aobjs)

        return anlysis

#==============================================================================
    @checkClosed
    def analysisInfo(self):
        """
        Get information about the Anaysis Object

        Returns
        -------
        dir:
            string specifying the directory for file I/O

        unitSys:
            string describing the unit system used by the AIM (can be None)
            
        major:
            the returned AIM major version number
            
        minor:
            the returned AIM minor version number
            
        intent:
            the intent string used to pass Bodies to the AIM
            
        fnames:
            a list of strings with the field/DataSet names
            
        franks:
            a list of ranks associated with each field

        fInOut:
            a list of FIELDIN/FieldOut flags associated with each field
            
        execute:
            execution flag: 0 - no exec, 1 - AIM Execute runs analysis, 2 - Auto Exec

        status:
            0 - up to date, 1 - dirty Analysis inputs, 2 - dirty Geometry inputs
            3 - both Geometry & Analysis inputs are dirty, 4 - new geometry,
            5 - post Analysis required, 6 -- Execution & post Analysis required
        """ 
        dir     = c_char_p()
        unitSys = c_char_p()
        major   = c_int()
        minor   = c_int()
        intent  = c_char_p()
        nfields = c_int()
        pfnames = POINTER(c_char_p)()
        pfranks = POINTER(c_int)()
        pfInOut = POINTER(c_int)()
        execute = c_int()
        status  = c_int()

        stat = _caps.caps_analysisInfo(self._obj, 
                                       ctypes.byref(dir), 
                                       ctypes.byref(unitSys), 
                                       ctypes.byref(major), ctypes.byref(minor), 
                                       ctypes.byref(intent), 
                                       ctypes.byref(nfields), ctypes.byref(pfnames), ctypes.byref(pfranks), ctypes.byref(pfInOut), 
                                       ctypes.byref(execute), ctypes.byref(status))
        if stat: _raiseStatus(stat)

        dir     = _decode(dir.value)     if dir     else None
        unitSys = _decode(unitSys.value) if unitSys else None
        major   = major.value
        minor   = minor.value
        intent  = _decode(intent.value) if intent and intent.value != b"" else None
        nfields = nfields.value
        fnames  = [_decode(pfnames[i]) for i in range(nfields)]
        franks  = [pfranks[i]          for i in range(nfields)]
        fInOut  = [pfInOut[i]          for i in range(nfields)]
        execute = execute.value
        status  = status.value

        if intent is not None and ";" in intent: intent = intent.split(";")

        return dir, unitSys, major, minor, intent, fnames, franks, fInOut, execute, status

#==============================================================================
    @checkClosed
    def preAnalysis(self):
        """
        Generate Analysis Inputs
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_preAnalysis(self._obj, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def system(self, cmd, rpath=None):
        """
        Execute the Command Line String
        
        Notes: 
            1. only needed when explicitly executing the appropriate analysis solver (i.e., not using the AIM)
            2. should be invoked after caps_preAnalysis and before caps_postAnalysis
            3. this must be used instead of the OS system call to ensure that journaling properly functions
                
        Parameters
        ----------
        self:
            CAPS Analysis Object

        cmd:
            the command line string to execute

        rpath:
            the relative path from the Analysis' directory or None (in the Analysis path)
        """
        cmd = cmd.encode() if isinstance(cmd,str) else cmd
        rpath = rpath.encode() if isinstance(rpath,str) else rpath
        stat = _caps.caps_system(self._obj, rpath, cmd)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def execute(self):
        """
        Execute Analysis if AIM does execution or AutoExec
        """
        status = c_int()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_execute(self._obj, ctypes.byref(status), ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def postAnalysis(self):
        """
        Mark Analysis as run
        
        Note: this clears all Analysis Output Objects to force reloads/recomputes
        """
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_postAnalysis(self._obj, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def boundInfo(self):
        """
        Get information about the Bound Object

        Returns
        -------
        state:
            the Bound state: 
            -1 - Open, 0 - Empty & Closed, 1 - single BRep entity, 
            2 - multiple BRep entities, -2 - multiple BRep entities with Errors in reparameterization

        dim:
            the dimensiounalty of the Bound (1 - 3
            
        plims:
            the parametrization limits  (2 values when dim is 1, 4 when dim is 2)
        """ 
        state   = c_int()
        dim     = c_int()
        plims  = POINTER(c_double)()
        stat = _caps.caps_boundInfo(self._obj, ctypes.byref(state), 
                                    ctypes.byref(dim), ctypes.byref(plims))
        if stat: _raiseStatus(stat)
        
        state   = state.value
        dim     = dim.value
        plims   = plims[:2*dim]
        
        return state, dim, plims

#==============================================================================
    @checkClosed
    def closeBound(self):
        """
        Completes the Bound Object
        """ 
        stat = _caps.caps_closeBound(self._obj)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def makeVertexSet(self, analysis, vname=None):
        """
        Create a VertexSet on the Bound Object

        Parameters
        ----------
        analysis:
            the Analysis Object (None - Unconnected) 

        vname:
            string naming the VertexSet (can be None for a Connected VertexSet)
        
        Returns
        -------
        the new VertexSet Object
        """ 
        vname = vname.encode() if isinstance(vname,str) else vname

        vset  = c_capsObj()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_makeVertexSet(self._obj, analysis._obj, 
                                        vname, ctypes.byref(vset),
                                        ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        return capsObj(vset, self.problemObj(), None)

#==============================================================================
    @checkClosed
    def getDataSets(self, dname):
        """
        Get DataSet Objects by name
        
        Parameters
        ----------
        dname:
            string containing the name of the DataSet 
        """
        dname = dname.encode() if isinstance(dname, str) else dname
        nobj  = c_int()
        dsets = POINTER(c_capsObj())
        stat = _caps.caps_getDataSets(self._obj, dname, ctypes.byref(nobj), ctypes.byref(dsets))
        if stat: _raiseStatus(stat)
        
        datasets = [DataSet(c_capsObj(dsets[i].contents)) for i in range(nobj.value)]
        egads.free(dsets)
        
        return datasets

#==============================================================================
    @checkClosed
    def getTriangles(self):
        """
        Get Triangulations for a 2D VertexSet
        
        Returns
        -------
        Gtris: 
            element mesh triangles: a list of 3-tuple indices (bias 1) referencing Geometry-based points

        Gsegs: 
            element mesh segments: a list of 2-tuple indices (bias 1) referencing Geometry-based points
    
        Dtris: 
            element data triangles: a list of 3-tuple indices (bias 1) referencing Data-based points

        Dsegs: 
            element data segments: a list of 2-tuple indices (bias 1) referencing Data-based points
        """
        nGtris  = c_int()
        pGtris  = POINTER(c_int)()
        nGsegs  = c_int()
        pGsegs  = POINTER(c_int)()
        nDtris  = c_int()
        pDtris  = POINTER(c_int)()
        nDsegs  = c_int()
        pDsegs  = POINTER(c_int)()
        stat = _caps.caps_getTriangles(self._obj, ctypes.byref(nGtris), ctypes.byref(pGtris), 
                                                  ctypes.byref(nGsegs), ctypes.byref(pGsegs), 
                                                  ctypes.byref(nDtris), ctypes.byref(pDtris),
                                                  ctypes.byref(nDsegs), ctypes.byref(pDsegs))
        if stat: _raiseStatus(stat)
        
        Gtris = None
        if pGtris:
            Gtris = [tuple(pGtris[3*i:3*(i+1)]) for i in range(nGtris.value)]

        Gsegs = None
        if pGsegs:
            Gsegs = [tuple(pGsegs[2*i:2*(i+1)]) for i in range(nGsegs.value)]

        Dtris = None
        if pDtris:
            Dtris = [tuple(pGtris[3*i:3*(i+1)]) for i in range(nDtris.value)]

        Dsegs = None
        if pDsegs:
            Dsegs = [tuple(pDsegs[2*i:2*(i+1)]) for i in range(nDsegs.value)]

        egads.free(pGtris)
        egads.free(pGsegs)
        egads.free(pDtris)
        egads.free(pDsegs)
        
        return Gtris, Gsegs, Dtris, Dsegs

#==============================================================================
    @checkClosed
    def vertexSetInfo(self):
        """
        Get information about the VertexSet Object

        Returns
        -------
        nGpts: 
            the returned number of Geometry points in the VertexSet
            
        nDpts:
            the returned number of point Data positions in the VertexSet
            
        bound:
            the returned associated Bound Object
            
        analysis:
            the returned associated Analysis Object (None -- Unconnected)
        """ 
        nGpts    = c_int()
        nDpts    = c_int()
        bound    = c_capsObj()
        analysis = c_capsObj()
        stat = _caps.caps_vertexSetInfo(self._obj, ctypes.byref(nGpts), ctypes.byref(nDpts), 
                                        ctypes.byref(bound), ctypes.byref(analysis))
        if stat: _raiseStatus(stat)
        
        nGpts    = nGpts.value
        nDpts    = nDpts.value
        bound    = Bound(bound, deleteObject=False)
        analysis = Analysis(analysis, deleteObject=False) if analysis else None
        
        return nGpts, nDpts, bound, analysis

#==============================================================================
    @checkClosed
    def fillUnVertexSet(self, xyzs):
        """
        Fill an Unconnected VertexSet
        
        Parameters
        ----------
        xyzs:
            list of 3 valued tuple point coorediantes
        """
        npts = len(xyzs)
        pxyzs = (c_double*(3*npts))()
        for i in range(npts):
            pxyzs[3*i:3*(i+1)] = xyzs[i][:]
        
        stat = _caps.caps_fillUnVertexSets(self._obj, c_int(npts), pxyzs)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def makeDataSet(self, dname, dmethod, rank=0):
        """
        Create a DataSet Object -- associated Bound must be Open
        
        Parameters
        ----------
        dname:
            a string containing the name of the DataSet (i.e., "pressure")
        
        dmethod:
            the method used for data transfers: (FieldIn, FieldOut, GeomSens, TessSens, User)
            
        rank:
            the rank of the data for a User field (e.g., 1 -- scalar, 3 -- vector), ignored otherwise
        """
        dname = dname.encode() if isinstance(dname,str) else dname

        dset = c_capsObj()
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_makeDataSet(self._obj, dname, c_int(dmethod), c_int(rank), ctypes.byref(dset),
                                      ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

        return capsObj(dset, self.problemObj(), None)

#==============================================================================
    @checkClosed
    def dataSetInfo(self):
        """
        Returns info about a DataSet Object
        
        Returns
        ----------
        ftype:
            the Field type of the DataSet: (FieldIn, FieldOut, GeomSens, TessSens, User)
        
        link:
            a linked DataSet Object or None

        dmethod:
            the method used for data transfers: 
        """
        dname = dname.encode() if isinstance(dname,str) else dname

        ftype = c_int()
        link = c_capsObj()
        dmethod = ()
        stat = _caps.caps_makeDataSet(self._obj, c_int(rank), ctypes.byref(link),
                                      ctypes.byref(dmethod))
        if stat: _raiseStatus(stat)

        link = capsObj(link, self.problemObj(), None) if link else None

        return ftype.value, link, dmethod.value

#==============================================================================
    @checkClosed
    def linkDataSet(self, link, dmethod):
        """
        Links data from DataSet Object
        
        self is the target DataSet Object which must be a FieldIn DataSet
        
        Parameters
        ----------
        link:
            the source DataSet Object

        dmethod:
            0 - Interpolate, 1 - Conserve
        """
        if not isinstance(link, capsObj):
            raise CAPSError(CAPS_BADVALUE, "link must be of type capsObj")

        if link._obj.contents.type != oType.DATASET:
            raise CAPSError(CAPS_BADVALUE, "link must be a DataSet Object")

        if self._obj.contents.type != oType.DATASET:
            raise CAPSError(CAPS_BADVALUE, "self must be a DataSet Object")

        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_linkDataSet(link._obj, c_int(dmethod), self._obj,
                                      ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def initDataSet(self, startup):
        """
        Initialize DataSet for cyclic/incremental startup

        Note: invocations of caps_getData and aim_getDataSet will 
        return this data (and a length of 1) until properly filled.
        
        Parameters
        ----------
        startup:
            the constant startup data (scalar or tuple)
        """
        rank = 1 if isinstance(startup,float) else len(startup)
        pstartup = (c_double*rank)()
        pstartup[:] = startup
        
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_initDataSet(self._obj, c_int(rank), pstartup,
                                      ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def getData(self):
        """
        Get data from the DataSet

        Returns
        -------
        data:
            the data (rank long tuples list, npts in legnth)
            
        units:
            string declaring the units
        """ 
        npt   = c_int()
        rank  = c_int()
        pdata = POINTER(c_double)()
        units = c_char_p()
        
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_getData(self._obj, ctypes.byref(npt), ctypes.byref(rank), 
                                  ctypes.byref(pdata), ctypes.byref(units), 
                                  ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))
        
        npt    = npt.value
        rank   = rank.value
        if rank == 1:
            data = [pdata[i] for i in range(npt)]
        else:
            data = [tuple(pdata[rank*i:rank*(i+1)]) for i in range(npt)]
        units  = _decode(units.value) if units else None
        
        return _withUnits(data, units)

#==============================================================================
    @checkClosed
    def setData(self, data):
        """
        Set data for the DataSet

        Parameters
        ----------
        data:
            the data (or Unit) (list of rank long tuples)
        """ 
        if isinstance(data, Quantity):
            units = data._units
            data  = data._value
        else:
            units = None
        
        npt   = len(data)
        rank  = 1 if isinstance(data[0],float) else len(data[0])
        pdata = POINTER(c_double*(rank*npt))()
        units = units.encode() if isinstance(units, str) else units
        
        for i in range(npt):
            pdata[rank*i:rank*(i+1)] = data[i]
        
        nErr = c_int()
        errs = POINTER(c_capsErrs)()
        stat = _caps.caps_setData(self._obj, c_int(npt), c_int(rank), 
                                  pdata, units, ctypes.byref(nErr), ctypes.byref(errs))
        if stat: _raiseStatus(stat, errors=capsErrs(nErr, errs))

#==============================================================================
    @checkClosed
    def fillUnVertexSet(self, xyzs):
        """
        Fill an Unconnected VertexSet
        
        Parameters
        ----------
        xyzs:
            list of 3 valued tuple point coorediantes
        """
        npts = len(xyzs)
        pxyzs = (c_double*(3*npts))()
        for i in range(npts):
            pxyzs[3*i:3*(i+1)] = xyzs[i][:]
        
        stat = _caps.caps_fillUnVertexSets(self._obj, c_int(npts), pxyzs)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def outputObjects(self):
        stat = _caps.caps_outputObjects(self._obj, None)
        if stat: _raiseStatus(stat)

#==============================================================================
    @checkClosed
    def AIMbackdoor(self, JSONin):
        
        JSONin = JSONin.encode() if isinstance(JSONin, str) else JSONin
        
        pJSONout = c_char_p()
        stat = _caps.caps_AIMbackdoor(self._obj, JSONin, ctypes.byref(pJSONout))
        if stat: _raiseStatus(stat)
        
        JSONout = _decode(pJSONout.value)
        egads.free(pJSONout)

        return JSONout


#==============================================================================
# Helper function to return Unit class or values
def _withUnits(value, units):
    if units is None:
        return value
    else:
        return Quantity(value, Unit(units))

#==============================================================================
__all__.append("Unit")
class Unit(object):
    """
    Represents a unit which can be manipulated into other units
    """

#==============================================================================
    def __init__(self, units):
        """
        Parameters
        ----------
        units:
            a udunits string representing the units
        """
        if isinstance(units, str):
            self._units = units
        elif isinstance(units, (int, float)):
            self._units = str(units)
        elif isinstance(units, Unit):
            self._units = units._units
        elif isinstance(units, Quantity):
            self._units = str(units)
        else:
            _raiseStatus(CAPS_UNITERR, "Unit.__init__ accepts str, Unit, or Quantity")

#==============================================================================
    def __hash__(self):
        return hash( (self.__class__.__name__, self._units) )

#==============================================================================
    def __repr__(self):
        '''x.__repr__() <==> repr(x)
        '''
        return '<{0} {1}>'.format(self.__class__.__name__, self)

#==============================================================================
    def __str__(self):
        '''x.__str__() <==> str(x)
        '''
        return self._units

#==============================================================================
    def __copy__(self):
        return self.__class__(self._unit)

    def __deepcopy__(self, memo):
        '''Used if copy.deepcopy is called on the variable.
        '''
        other = self.__copy__()
        memo[id(self)] = other
        return other

#==============================================================================
    def __eq__(self, other):
        '''Comparison operator: x.__eq__(y) <==> x == y
        '''
        if isinstance(other, self.__class__):
            compare = c_int()
            stat = _caps.caps_unitCompare(self._units.encode(), other._units.encode(), ctypes.byref(compare))
            if stat: _raiseStatus(stat, "Cannot {!r} == {!r}".format(self, other))
            return compare.value == 0
        else:
            return False

#==============================================================================
    def __ne__(self, other):
        '''Comparison operator: x.__ne__(y) <==> x != y
        '''
        return not self.__eq__(other)

#==============================================================================
    def __gt__(self, other):
        '''Comparison operator: x.__gt__(y) <==> x > y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} > {!r}".format(self, other))

        compare = c_int()
        stat = _caps.caps_unitCompare(self._units.encode(), other._units.encode(), ctypes.byref(compare))
        if stat: _raiseStatus(stat, "Cannot {!r} > {!r}".format(self, other))
        return compare.value > 0

#==============================================================================
    def __ge__(self, other):
        '''Comparison operator: x.__ge__(y) <==> x >= y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} >= {!r}".format(self, other))

        compare = c_int()
        stat = _caps.caps_unitCompare(self._units.encode(), other._units.encode(), ctypes.byref(compare))
        if stat: _raiseStatus(stat, "Cannot {!r} >= {!r}".format(self, other))
        return compare.value >= 0

#==============================================================================
    def __lt__(self, other):
        '''Comparison operator: x.__lt__(y) <==> x < y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} < {!r}".format(self, other))

        compare = c_int()
        stat = _caps.caps_unitCompare(self._units.encode(), other._units.encode(), ctypes.byref(compare))
        if stat: _raiseStatus(stat, "Cannot {!r} < {!r}".format(self, other))
        return compare.value < 0

#==============================================================================
    def __le__(self, other):
        '''Comparison operator: x.__le__(y) <==> x <= y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} <= {!r}".format(self, other))

        compare = c_int()
        stat = _caps.caps_unitCompare(self._units.encode(), other._units.encode(), ctypes.byref(compare))
        if stat: _raiseStatus(stat, "Cannot {!r} <= {!r}".format(self, other))
        return compare.value <= 0

#==============================================================================
    def __sub__(self, other):
        '''Binary operator: x.__sub__(y) <==> x - y
        '''
        if not isinstance(other, (int, float)):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} - {!r}".format(self, other))

        newUnits = c_char_p()
        stat = _caps.caps_unitOffset(self._units.encode(), c_double(-other), 
                                     ctypes.byref(newUnits))
        if stat: _raiseStatus(stat, "Cannot {!r} - {!r}".format(self, other))
        units = _decode(newUnits.value)
        egads.free(newUnits)

        return Unit(units)

#==============================================================================
    def __add__(self, other):
        '''Binary operator: x.__add__(y) <==> x + y
        '''
        if not isinstance(other, (int, float)):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} + {!r}".format(self, other))

        newUnits = c_char_p()
        stat = _caps.caps_unitOffset(self._units.encode(), c_double(other), 
                                     ctypes.byref(newUnits))
        if stat: _raiseStatus(stat, "Cannot {!r} + {!r}".format(self, other))
        units = _decode(newUnits.value)
        egads.free(newUnits)

        return Unit(units)

#==============================================================================
    def __mul__(self, other):
        '''Binary operator: x.__mul__(y) <==> x * y
        '''
        if isinstance(other, self.__class__):
            
            newUnits = c_char_p()
            stat = _caps.caps_unitMultiply(self._units.encode(), other._units.encode(), 
                                           ctypes.byref(newUnits))
            if stat: _raiseStatus(stat, "Cannot {!r} * {!r}".format(self, other))
            units = _decode(newUnits.value)
            egads.free(newUnits)

            return Unit(units)
        else:
            return Quantity(other, self)

        #_raiseStatus(CAPS_BADVALUE, "Cannot {!r} * {!r}".format(self, other))

#==============================================================================
    def __div__(self, other):
        '''Binary operator: x.__div__(y) <==> x/y
        '''
        if isinstance(other, self.__class__):

            newUnits = c_char_p()
            stat = _caps.caps_unitDivide(self._units.encode(), other._units.encode(), 
                                         ctypes.byref(newUnits))
            if stat: _raiseStatus(stat, "Cannot {!r} / {!r}".format(self, other))
            units = _decode(newUnits.value)
            egads.free(newUnits)

            return Unit(units)
        else:
            return Quantity(1/other, 1/self)

        #_raiseStatus(CAPS_BADVALUE, "Cannot {!r} / {!r}".format(self, other))

#==============================================================================
    def __pow__(self, other, modulo=None):
        '''Binary operator: x.__pow__(y) <==> x ** y
        '''
        if modulo is not None:
            raise CAPSError(CAPS_UNITERR, "3-argument power not supported")

        if not isinstance(other, int):
            raise CAPSError(CAPS_UNITERR, "Unit can only be raised by integer powers! Power = {!r}".format(other))

        newUnits = c_char_p()
        stat = _caps.caps_unitRaise(self._units.encode(), c_int(other), 
                                     ctypes.byref(newUnits))
        if stat: _raiseStatus(stat, "Cannot {!r} ** {!r}".format(self, other))
        units = _decode(newUnits.value)
        egads.free(newUnits)

        return Unit(units)

#==============================================================================
    def __isub__(self, other):
        '''In-place operator: x.__isub__(y) <==> x -= y
        '''
        try:
            newUnit = self - other
            self._units = newUnit._units
            return self
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} -= {!r}".format(self, other))

#==============================================================================
    def __iadd__(self, other):
        '''In-place operator: x.__iadd__(y) <==> x += y
        '''
        try:
            newUnit = self + other
            self._units = newUnit._units
            return self
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} += {!r}".format(self, other))

#==============================================================================
    def __imul__(self, other):
        '''In-place operator: x.__imul__(y) <==> x *= y
        '''
        try:
            newUnit = self * other
            self._units = newUnit._units
            return self
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} *= {!r}".format(self, other))

#==============================================================================
    def __idiv__(self, other):
        '''In-place operator: x.__idiv__(y) <==> x /= y
        '''
        try:
            newUnit = self / other
            self._units = newUnit._units
            return self
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} /= {!r}".format(self, other))

#==============================================================================
    def __ipow__(self, other):
        '''In-place operator: x.__ipow__(y) <==> x **= y
        '''
        try:
            newUnit = self ** other
            self._units = newUnit._units
            return self
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} **= {!r}".format(self, other))

#==============================================================================
    def __rsub__(self, other):
        '''Reflected binary operator: x.__rsub__(y) <==> y - x
        '''
        try:
            return -self + other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} - {!r}".format(other, self))

#==============================================================================
    def __radd__(self, other):
        '''Reflected binary operator: x.__radd__(y) <==> y + x
        '''
        return self + other

#==============================================================================
    def __rmul__(self, other):
        '''Reflected binary operator: x.__rmul__(y) <==> y * x
        '''
        return self * other

#==============================================================================
    def __rdiv__(self, other):
        '''Reflected binary operator: x.__rdiv__(y) <==> y / x
        '''
        try:
            return (self ** -1) * other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} / {!r}".format(other, self))

#==============================================================================
    def __floordiv__(self, other):
        '''x.__floordiv__(y) <==> x//y <==> x / y
        '''
        return self / other

#==============================================================================
    def __ifloordiv__(self, other):
        '''x.__ifloordiv__(y) <==> x//=y <==> x /= y
        '''
        return self / other

#==============================================================================
    def __rfloordiv__(self, other):
        '''x.__rfloordiv__(y) <==> y//x <==> y / x
        '''
        try:
            return (self ** -1) * other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} // {!r}".format(other, self))

#==============================================================================
    def __truediv__(self, other):
        '''x.__truediv__(y) <==> x / y
        '''
        return self.__div__(other)

#==============================================================================
    def __itruediv__(self, other):
        '''x.__itruediv__(y) <==> x /= y
        '''
        return self.__idiv__(other)

#==============================================================================
    def __rtruediv__(self, other):
        '''x.__rtruediv__(y) <==> y / x
        '''
        return self.__rdiv__(other)

#==============================================================================
    def __mod__(self, other):
        '''x.__mod__(y) <==> y % x
        '''
        _raiseStatus(CAPS_UNITERR, "Cannot {!r} % {!r}".format(self, other))

#==============================================================================
    def __neg__(self):
        '''Unary operator: x.__neg__() <==> -x
        '''
        return self * Unit(-1)

#==============================================================================
    def __pos__(self):
        '''Unary operator: x.__pos__() <==> +x
        '''
        return self


#==============================================================================
__all__.append("Quantity")
class Quantity(object):
    """
    Represents a value with units
    """

    #
    # Make Quantity dominant when multiplying with a numpy.array for left or right
    #
    __array_priority__ = 20.0

#==============================================================================
    def __init__(self, value, units):
        """
        Parameters
        ----------
        value:
            the scalar value

        units:
            the units of value
        """
        self._value = value
        if isinstance(units,str):
            self._units = Unit(units)
        elif isinstance(units,Unit):
            self._units = units
        else:
            raise CAPSError(CAPS_UNITERROR, "units must be str or Unit")

#==============================================================================
    def value(self):
        """
        Returns the Value in the Quantity
        """
        return self._value

#==============================================================================
    def unit(self):
        """
        Returns the Unit in the Quantity
        """
        return self._unit

#==============================================================================
    def convert(self, toUnits):
        """
        Convertes to other units
        
        Parameters
        ----------
        toUnits:
            Unit to convert to
        
        Returns
        -------
        quantity converted to toUnits
        """
        if not isinstance(toUnits,(Unit,str)):
            if stat: _raiseStatus(CAPS_UNITERR, "toUnits mist be instance of Unit or str. toUnits = {!r}".format(toUnits))
        
        if isinstance(toUnits,str):
            toUnits = Unit(toUnits)
        
        tounits = toUnits._units.encode()

        fromunits = self._units._units.encode()
        val = c_double()

        if hasattr(self._value, '__len__'):
            count = len(self._value)
            fromVal = (c_double*count)()
            toVal   = (c_double*count)()
            
            fromVal[:] = self._value[:]

            stat = _caps.caps_convert(c_int(count), fromunits, fromVal,
                                                    tounits  , toVal  )
            if stat: _raiseStatus(stat, "Cannot convert {!r} to {!r}".format(self, toUnits))
   
            value    = copy.copy(self._value)
            value[:] = toVal[:]
        else:
            fromVal = c_double(self._value)
            toVal   = c_double()

            stat = _caps.caps_convert(c_int(1), fromunits, ctypes.byref(fromVal),
                                                tounits  , ctypes.byref(toVal)  )
            if stat: _raiseStatus(stat, "Cannot convert {!r} to {!r}".format(self, toUnits))
            
            value = toVal.value
   
        return Quantity(value, toUnits)

#==============================================================================
    def __hash__(self):
        return hash( (self.__class__.__name__, self._value, self._units) )

#==============================================================================
    def __repr__(self):
        '''x.__repr__() <==> repr(x)
        '''
        return '<{0} {1}>'.format(self.__class__.__name__, self)

#==============================================================================
    def __str__(self):
        '''x.__str__() <==> str(x)
        '''
        return str(self._value) + ' ' + self._units._units

#==============================================================================
    def __copy__(self):
        return self.__class__(self._value, self._unit)

    def __deepcopy__(self, memo):
        '''Used if copy.deepcopy is called on the variable.
        '''
        other = self.__copy__()
        memo[id(self)] = other
        return other

#==============================================================================
    def __bool__(self):
        '''Truth value testing and the built-in operation ``bool``

            x.__bool__() <==> x != 0

        '''
        return bool(self._value)

#==============================================================================
    def __float__(self): # conversion to float (makes trig functions work)
        return self / self.__class__(1, "rad")

#==============================================================================
    def __len__(self):
        return len(self._value)

#==============================================================================
    def __getitem__(self,index):
        ''' returns self's value sliced to index with self's units
        '''
        if hasattr(self._value, '__len__'):
            return self.__class__(self._value[index], self._units)
        elif index == 0 or index == -1:
            return self.__class__(self._value, self._units)
        else:
            raise IndexError

#==============================================================================
    def __setitem__(self,index,other):
        ''' makes a slice assignment on self's value based on index
        '''
        other = other.convert(self._units)
        if hasattr(self._value, '__len__'):
            self._value[index] = other._value
        elif index == 0 or index == -1:
           self._value = other._value
        else:
           raise IndexError

#==============================================================================
    def __eq__(self, other):
        '''Comparison operator: x.__eq__(y) <==> x == y
        '''
        if isinstance(other, self.__class__):
            try:
                other = other.convert(self._units)
            except CAPSError:
                return False

            return self._value == other._value and \
                   self._units == other._units
        else:
            if other == 0:
                return self._value == 0
            return False

#==============================================================================
    def __ne__(self, other):
        '''Comparison operator: x.__ne__(y) <==> x != y
        '''
        return not self.__eq__(other)

#==============================================================================
    def __gt__(self, other):
        '''Comparison operator: x.__gt__(y) <==> x > y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} > {!r}".format(self, other))

        other = other.convert(self._units)

        return self._value > other._value

#==============================================================================
    def __ge__(self, other):
        '''Comparison operator: x.__ge__(y) <==> x >= y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} >= {!r}".format(self, other))

        other = other.convert(self._units)

        return self._value >= other._value

#==============================================================================
    def __lt__(self, other):
        '''Comparison operator: x.__lt__(y) <==> x < y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} < {!r}".format(self, other))

        other = other.convert(self._units)

        return self._value < other._value

#==============================================================================
    def __le__(self, other):
        '''Comparison operator: x.__le__(y) <==> x <= y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} <= {!r}".format(self, other))

        other = other.convert(self._units)

        return self._value <= other._value

#==============================================================================
    def __sub__(self, other):
        '''Binary operator: x.__sub__(y) <==> x - y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} - {!r}".format(self, other))

        other = other.convert(self._units)

        return self.__class__(self._value - other._value, self._units)

#==============================================================================
    def __add__(self, other):
        '''Binary operator: x.__add__(y) <==> x + y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} + {!r}".format(self, other))

        other = other.convert(self._units)

        return self.__class__(self._value + other._value, self._units)

#==============================================================================
    def __mul__(self, other):
        '''Binary operator: x.__mul__(y) <==> x * y
        '''
        if isinstance(other, Unit):
            other = 1 * other

        if isinstance(other, self.__class__):
            
            units = self._units * other._units
        
            if hasattr(self._value, '__len__'):
                value = copy.copy(self._value)
                for i in range(len(self._value)):
                    value[i] = self._value[i] * other._value
            elif hasattr(other._value, '__len__'):
                value = copy.copy(other._value)
                for i in range(len(other._value)):
                    value[i] = self._value * other._value[i]
            else:
                value = self._value*other._value

            if " " in units._units:
                split = units._units.split(" ")
                if split[1] == "1":
                    scale = float(split[0])
                    if hasattr(value, '__len__'):
                        for i in range(len(value)):
                            value[i] *= scale
                    else:
                        value *= scale
                    return value
            elif units._units == "1":
                return value

            return self.__class__(value, units)

        elif isinstance(other, (float,int)):

            return self.__class__(self._value*other, self._units)

        elif hasattr(other, '__len__'):
            
            value = copy.copy(other)
            for i in range(len(other)):
                value[i] = self._value * other[i]

            return self.__class__(value, self._units)

        _raiseStatus(CAPS_BADVALUE, "Cannot {!r} * {!r}".format(self, other))

#==============================================================================
    def __div__(self, other):
        '''Binary operator: x.__div__(y) <==> x/y
        '''
        if isinstance(other, Unit):
            other = 1 * other

        if isinstance(other, self.__class__):

            units = self._units / other._units

            if hasattr(self._value, '__len__'):
                value = copy.copy(self._value)
                for i in range(len(self._value)):
                    value[i] = self._value[i] / other._value
            else:
                value = self._value/other._value
                
            if " " in units._units:
                split = units._units.split(" ")
                if split[1] == "1":
                    scale = float(split[0])
                    if hasattr(value, '__len__'):
                        for i in range(len(value)):
                            value[i] *= scale
                    else:
                        value *= scale
                    return value
            elif units._units == "1":
                return value

            return self.__class__(value, units)

        elif isinstance(other, (float,int)):

            return self.__class__(self._value/other, self._units)

        _raiseStatus(CAPS_UNITERR, "Cannot {!r} / {!r}".format(self, other))

#==============================================================================
    def __pow__(self, other, modulo=None):
        '''Binary operator: x.__pow__(y) <==> x ** y
        '''
        if modulo is not None:
            raise CAPSError(CAPS_UNITERR, "3-argument power not supported")

        if not isinstance(other, int):
            raise CAPSError(CAPS_UNITERR, "Unit can only be raised by integer powers! Power = {!r}".format(other))

        return self.__class__(self._value**other, self._units**other)

#==============================================================================
    def __isub__(self, other):
        '''In-place operator: x.__isub__(y) <==> x -= y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} -= {!r}".format(self, other))

        other = other.convert(self._units)
        self._value -= other._value

        return self

#==============================================================================
    def __iadd__(self, other):
        '''In-place operator: x.__iadd__(y) <==> x += y
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} += {!r}".format(self, other))

        other = other.convert(self._units)
        self._value += other._value

        return self

#==============================================================================
    def __imul__(self, other):
        '''In-place operator: x.__imul__(y) <==> x *= y
        '''
        try:
            newval = self * other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} *= {!r}".format(self, other))

        self._value = newval._value
        self._units = newval._units
        return self

#==============================================================================
    def __idiv__(self, other):
        '''In-place operator: x.__idiv__(y) <==> x /= y
        '''
        try:
            newval = self / other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} /= {!r}".format(self, other))

        self._value = newval._value
        self._units = newval._units
        return self

#==============================================================================
    def __ipow__(self, other):
        '''In-place operator: x.__ipow__(y) <==> x **= y
        '''
        try:
            newval = self ** other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} **= {!r}".format(self, other))

        self._value = newval._value
        self._units = newval._units
        return self

#==============================================================================
    def __rsub__(self, other):
        '''Reflected binary operator: x.__rsub__(y) <==> y - x
        '''
        try:
            return -self + other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} - {!r}".format(other, self))

#==============================================================================
    def __radd__(self, other):
        '''Reflected binary operator: x.__radd__(y) <==> y + x
        '''
        return self + other

#==============================================================================
    def __rmul__(self, other):
        '''Reflected binary operator: x.__rmul__(y) <==> y * x
        '''
        return self * other

#==============================================================================
    def __rdiv__(self, other):
        '''Reflected binary operator: x.__rdiv__(y) <==> y / x
        '''
        try:
            return (self ** -1) * other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} / {!r}".format(other, self))

#==============================================================================
    def __floordiv__(self, other):
        '''x.__floordiv__(y) <==> x//y <==> x / y
        '''
        return self / other

#==============================================================================
    def __ifloordiv__(self, other):
        '''x.__ifloordiv__(y) <==> x//=y <==> x /= y
        '''
        return self / other

#==============================================================================
    def __rfloordiv__(self, other):
        '''x.__rfloordiv__(y) <==> y//x <==> y / x
        '''
        try:
            return (self ** -1) * other
        except CAPSError as e:
            raise CAPSError(e.errorCode, "Cannot {!r} // {!r}".format(other, self))

#==============================================================================
    def __truediv__(self, other):
        '''x.__truediv__(y) <==> x / y
        '''
        return self.__div__(other)

#==============================================================================
    def __itruediv__(self, other):
        '''x.__itruediv__(y) <==> x /= y
        '''
        return self.__idiv__(other)

#==============================================================================
    def __rtruediv__(self, other):
        '''x.__rtruediv__(y) <==> y / x
        '''
        return self.__rdiv__(other)

#==============================================================================
    def __mod__(self, other):
        '''x.__mod__(y) <==> y % x
        '''
        if not isinstance(other, self.__class__):
            _raiseStatus(CAPS_UNITERR, "Cannot {!r} % {!r}".format(self, other))

        other = other.convert(self._units)

        return  self.__class__(self._value % other._value, self._units)

#==============================================================================
    def __round__(self, n=None):
        '''Unary operator: x.__round__(n) <==> round(x, n)
        '''
        return self.__class__(round(self._value, n), self._units)

#==============================================================================
    def __abs__(self):
        '''Unary operator: x.__abs__() <==> abs(x)
        '''
        return self.__class__(abs(self._value), self._units)

#==============================================================================
    def __neg__(self):
        '''Unary operator: x.__neg__() <==> -x
        '''
        return self * -1

#==============================================================================
    def __pos__(self):
        '''Unary operator: x.__pos__() <==> +x
        '''
        return self

#==============================================================================
# class units:
#     # Collection of all available units
#     pass
# 
# UDUNITS2_XML_PATH = os.environ['UDUNITS2_XML_PATH']
# tree = xmlElementTree.parse(UDUNITS2_XML_PATH)
# root = tree.getroot()
# prefix = []
# 
# for imp in root:
#     tree = xmlElementTree.parse(os.path.dirname(UDUNITS2_XML_PATH) + os.sep + imp.text)
#     for pfx in tree.findall('prefix'):
#         for sym in pfx.iter('symbol'):
#             try:
#                 prefix.append(_decode(sym.text.encode('ascii')))
#             except UnicodeEncodeError:
#                 continue
# 
# for imp in root:
#     tree = xmlElementTree.parse(os.path.dirname(UDUNITS2_XML_PATH) + os.sep + imp.text)
#     for ut in tree.findall('unit'):
#         for sym in ut.iter('symbol'):
#             try:
#                 symbol = _decode(sym.text.encode('ascii'))
#             except UnicodeEncodeError:
#                 continue
#             setattr(units, symbol, unit(symbol))
#             for pfx in prefix:
#                 setattr(units, pfx+symbol, unit(pfx+symbol))
