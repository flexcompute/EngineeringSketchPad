###########################################################################
#                                                                         #
# pyEGADS --- Python version of EGADS API                                 #
#                                                                         #
#                                                                         #
#      Copyright 2011-2024, Massachusetts Institute of Technology         #
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
from sys import version_info

from . import egads_common
from .egads_common import *
from .egads_common import _ESP_ROOT, _egads_error_codes

if hasattr(egads_common, "_egads"):
    raise ImportError("egads has already been loaded! egads and egadslite cannot be loaded at the same time.")

__all__ = egads_common._all_

# load the shared library
if sys.platform.startswith('darwin'):
    _egadslite = ctypes.CDLL(_ESP_ROOT + "/lib/libegadslite.dylib")
elif sys.platform.startswith('linux'):
    _egadslite = ctypes.CDLL(_ESP_ROOT + "/lib/libegadslite.so")
elif sys.platform.startswith('win32'):
    if version_info.major == 3 and version_info.minor < 8:
        _egadslite = ctypes.CDLL(_ESP_ROOT + "\\lib\\egadslite.dll")
    else:
        _egadslite = ctypes.CDLL(_ESP_ROOT + "\\lib\\egadslite.dll", winmode=0)
else:
    raise IOError("Unknown platform: " + sys.platform)

# =============================================================================
#
# egads(lite).h functions
#
# Functions removed from egadslite are commented out for documentation purposes
#
# =============================================================================

# egads functions declarations
_egadslite.EG_free.argtypes = [c_void_p]
_egadslite.EG_free.restype = None

_egadslite.EG_revision.argtypes = [POINTER(c_int), POINTER(c_int), POINTER(c_char_p)]
_egadslite.EG_revision.restype = None

_egadslite.EG_open.argtypes = [POINTER(c_ego)]
_egadslite.EG_open.restype = c_int

_egadslite.EG_importModel.argtypes = [c_ego, c_size_t, c_char_p, POINTER(c_ego)]
_egadslite.EG_importModel.restype = c_int

_egadslite.EG_deleteObject.argtypes = [c_ego]
_egadslite.EG_deleteObject.restype = c_int

#_egadslite.EG_makeTransform.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_makeTransform.restype = c_int
delattr(Context, "makeTransform")

#_egadslite.EG_getTransformation.argtypes = [c_ego, POINTER(c_double)]
#_egadslite.EG_getTransformation.restype = c_int
delattr(ego, "getTransform")

_egadslite.EG_getContext.argtypes = [c_ego, POINTER(c_ego)]
_egadslite.EG_getContext.restype = c_int

_egadslite.EG_setOutLevel.argtypes = [c_ego, c_int]
_egadslite.EG_setOutLevel.restype = c_int

#_egadslite.EG_updateThread.argtypes = [c_ego]
#_egadslite.EG_updateThread.restype = c_int
delattr(Context, "updateThread")

_egadslite.EG_getInfo.argtypes = [c_ego, POINTER(c_int), POINTER(c_int), POINTER(c_ego), POINTER(c_ego), POINTER(c_ego)]
_egadslite.EG_getInfo.restype = c_int

#_egadslite.EG_copyObject.argtypes = [c_ego, POINTER(None), POINTER(c_ego)]
#_egadslite.EG_copyObject.restype = c_int
delattr(ego, "copyObject")

#_egadslite.EG_flipObject.argtypes = [c_ego, POINTER(c_ego)]
#_egadslite.EG_flipObject.restype = c_int
delattr(ego, "flipObject")

_egadslite.EG_close.argtypes = [c_ego]
_egadslite.EG_close.restype = c_int

#_egadslite.EG_setUserPointer.argtypes = [c_ego, POINTER(None)]
#_egadslite.EG_setUserPointer.restype = c_int
#delattr(ego, "setUserPointer")

#_egadslite.EG_getUserPointer.argtypes = [c_ego, POINTER(POINTER(None))]
#_egadslite.EG_getUserPointer.restype = c_int
#delattr(ego, "getUserPointer")

#_egadslite.EG_setFullAttrs.argtypes = [c_ego, c_int]
#_egadslite.EG_setFullAttrs.restype = c_int
delattr(Context, "setFullAttrs")

#_egadslite.EG_attributeAdd.argtypes = [c_ego, c_char_p, c_int, c_int, POINTER(c_int), POINTER(c_double), c_char_p]
#_egadslite.EG_attributeAdd.restype = c_int
delattr(ego, "attributeAdd")

#_egadslite.EG_attributeDel.argtypes = [c_ego, c_char_p]
#_egadslite.EG_attributeDel.restype = c_int
delattr(ego, "attributeDel")

_egadslite.EG_attributeNum.argtypes = [c_ego, POINTER(c_int)]
_egadslite.EG_attributeNum.restype = c_int

_egadslite.EG_attributeGet.argtypes = [c_ego, c_int, POINTER(c_char_p), POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egadslite.EG_attributeGet.restype = c_int

_egadslite.EG_attributeRet.argtypes = [c_ego, c_char_p, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egadslite.EG_attributeRet.restype = c_int

#_egadslite.EG_attributeDup.argtypes = [c_ego, c_ego]
#_egadslite.EG_attributeDup.restype = c_int
delattr(ego, "attributeDup")

#_egadslite.EG_attributeAddSeq.argtypes = [c_ego, c_char_p, c_int, c_int, POINTER(c_int), POINTER(c_double), c_char_p]
#_egadslite.EG_attributeAddSeq.restype = c_int
delattr(ego, "attributeAddSeq")

_egadslite.EG_attributeNumSeq.argtypes = [c_ego, c_char_p, POINTER(c_int)]
_egadslite.EG_attributeNumSeq.restype = c_int

_egadslite.EG_attributeRetSeq.argtypes = [c_ego, c_char_p, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double)), POINTER(c_char_p)]
_egadslite.EG_attributeRetSeq.restype = c_int

_egadslite.EG_getGeometry.argtypes = [c_ego, POINTER(c_int), POINTER(c_int), POINTER(c_ego), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
_egadslite.EG_getGeometry.restype = c_int

_egadslite.EG_getGeometryLen.argtypes = [c_ego, POINTER(c_int), POINTER(c_int)]
_egadslite.EG_getGeometryLen.restype = None

#_egadslite.EG_makeGeometry.argtypes = [c_ego, c_int, c_int, c_ego, POINTER(c_int), POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_makeGeometry.restype = c_int
delattr(Context, "makeGeometry")

_egadslite.EG_getRange.argtypes = [c_ego, POINTER(c_double), POINTER(c_int)]
_egadslite.EG_getRange.restype = c_int

_egadslite.EG_evaluate.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
_egadslite.EG_evaluate.restype = c_int

_egadslite.EG_invEvaluate.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egadslite.EG_invEvaluate.restype = c_int

_egadslite.EG_invEvaluateGuess.argtypes = [c_ego, POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egadslite.EG_invEvaluateGuess.restype = c_int

_egadslite.EG_arcLength.argtypes = [c_ego, c_double, c_double, POINTER(c_double)]
_egadslite.EG_arcLength.restype = c_int

_egadslite.EG_curvature.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
_egadslite.EG_curvature.restype = c_int

#_egadslite.EG_approximate.argtypes = [c_ego, c_int, c_double, POINTER(c_int), POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_approximate.restype = c_int
delattr(Context, "approximate")

#_egadslite.EG_fitTriangles.argtypes = [c_ego, c_int, POINTER(c_double), c_int, POINTER(c_int), POINTER(c_int), c_double, POINTER(c_ego)]
#_egadslite.EG_fitTriangles.restype = c_int
delattr(Context, "fitTriangles")

#_egadslite.EG_addKnots.argtypes = [c_ego, c_int, c_int, POINTER(c_double), c_int, c_int, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_addKnots.restype = c_int
delattr(ego, "addKnots")

#_egadslite.EG_otherCurve.argtypes = [c_ego, c_ego, c_double, POINTER(c_ego)]
#_egadslite.EG_otherCurve.restype = c_int
delattr(ego, "otherCurve")

#_egadslite.EG_isSame.argtypes = [c_ego, c_ego]
#_egadslite.EG_isSame.restype = c_int
delattr(ego, "isSame")

#_egadslite.EG_isoCline.argtypes = [c_ego, c_int, c_double, POINTER(c_ego)]
#_egadslite.EG_isoCline.restype = c_int
delattr(ego, "isoCline")

#_egadslite.EG_isIsoPCurve.argtypes = [c_ego, POINTER(c_int), POINTER(c_double), POINTER(c_int)]
#_egadslite.EG_isIsoPCurve.restype = c_int
delattr(ego, "isIsoPCurve")

#_egads.EG_mergeBSplineCurves.argtypes = [c_ego, c_ego, POINTER(c_ego)]
#_egads.EG_mergeBSplineCurves.restype = c_int
delattr(ego, "mergeBSplineCurves")

#_egadslite.EG_convertToBSpline.argtypes = [c_ego, POINTER(c_ego)]
#_egadslite.EG_convertToBSpline.restype = c_int
delattr(ego, "convertToBSpline")

#_egadslite.EG_convertToBSplineRange.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_convertToBSplineRange.restype = c_int
delattr(ego, "convertToBSplineRange")

#_egadslite.EG_skinning.argtypes = [c_int, POINTER(c_ego), c_int, POINTER(c_ego)]
#_egadslite.EG_skinning.restype = c_int
del globals()["skinning"]

_egadslite.EG_tolerance.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_tolerance.restype = c_int

_egadslite.EG_getTolerance.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_getTolerance.restype = c_int

#_egadslite.EG_setTolerance.argtypes = [c_ego, c_double]
#_egadslite.EG_setTolerance.restype = c_int

_egadslite.EG_getTopology.argtypes = [c_ego, POINTER(c_ego), POINTER(c_int), POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(POINTER(c_int))]
_egadslite.EG_getTopology.restype = c_int

#_egadslite.EG_makeTopology.argtypes = [c_ego, c_ego, c_int, c_int, POINTER(c_double), c_int, POINTER(c_ego), POINTER(c_int), POINTER(c_ego)]
#_egadslite.EG_makeTopology.restype = c_int
delattr(Context, "makeTopology")

#_egadslite.EG_makeLoop.argtypes = [c_int, POINTER(c_ego), c_ego, c_double, POINTER(c_ego)]
#_egadslite.EG_makeLoop.restype = c_int
del globals()["makeLoop"]

#_egadslite.EG_getArea.argtypes = [c_ego, POINTER(c_double), POINTER(c_double)]
#_egadslite.EG_getArea.restype = c_int
delattr(ego, "getArea")

#_egadslite.EG_makeFace.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_makeFace.restype = c_int
delattr(ego, "makeFace")

_egadslite.EG_getBodyTopos.argtypes = [c_ego, c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_ego))]
_egadslite.EG_getBodyTopos.restype = c_int

_egadslite.EG_indexBodyTopo.argtypes = [c_ego, c_ego]
_egadslite.EG_indexBodyTopo.restype = c_int

_egadslite.EG_objectBodyTopo.argtypes = [c_ego, c_int, c_int, POINTER(c_ego)]
_egadslite.EG_objectBodyTopo.restype = c_int

_egadslite.EG_inTopology.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_inTopology.restype = c_int

_egadslite.EG_inFace.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_inFace.restype = c_int

_egadslite.EG_getEdgeUV.argtypes = [c_ego, c_ego, c_int, c_double, POINTER(c_double)]
_egadslite.EG_getEdgeUV.restype = c_int

_egadslite.EG_getEdgeUVs.argtypes = [c_ego, c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double)]
_egadslite.EG_getEdgeUVs.restype = c_int

#_egadslite.EG_getPCurve.argtypes = [c_ego, c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
#_egadslite.EG_getPCurve.restype = c_int
#delattr(ego, "getPCurve")

_egadslite.EG_getWindingAngle.argtypes = [c_ego, c_double, POINTER(c_double)]
_egadslite.EG_getWindingAngle.restype = c_int

_egadslite.EG_getBody.argtypes = [c_ego, POINTER(c_ego)]
_egadslite.EG_getBody.restype = c_int

#_egadslite.EG_makeSolidBody.argtypes = [c_ego, c_int, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_makeSolidBody.restype = c_int
delattr(Context, "makeSolidBody")

_egadslite.EG_getBoundingBox.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_getBoundingBox.restype = c_int

#_egadslite.EG_getMassProperties.argtypes = [c_ego, POINTER(c_double)]
#_egadslite.EG_getMassProperties.restype = c_int
delattr(ego, "getMassProperties")

#_egadslite.EG_isEquivalent.argtypes = [c_ego, c_ego]
#_egadslite.EG_isEquivalent.restype = c_int
delattr(ego, "isEquivalent")

#_egadslite.EG_sewFaces.argtypes = [c_int, POINTER(c_ego), c_double, c_int, POINTER(c_ego)]
#_egadslite.EG_sewFaces.restype = c_int
del globals()["sewFaces"]

#_egadslite.EG_makeNmWireBody.argtypes = [ c_int, POINTER(c_ego), c_double, POINTER(c_ego)]
#_egadslite.EG_makeNmWireBody.restype = c_int
del globals()["makeNmWireBody"]

#_egadslite.EG_replaceFaces.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
#_egadslite.EG_replaceFaces.restype = c_int
delattr(ego, "replaceFaces")

#_egads.EG_removeNodes.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
#_egads.EG_removeNodes.restype = c_int
delattr(ego, "removeNodes")

#_egadslite.EG_mapBody.argtypes = [c_ego, c_ego, c_char_p, POINTER(c_ego)]
#_egadslite.EG_mapBody.restype = c_int
delattr(ego, "mapBody")

#_egadslite.EG_mapBody2.argtypes = [c_ego, c_char_p, c_char_p, POINTER(c_ego)]
#_egadslite.EG_mapBody2.restype = c_int
delattr(ego, "mapBody2")

#_egadslite.EG_matchBodyEdges.argtypes = [c_ego, c_ego, c_double, POINTER(c_int), POINTER(POINTER(c_int))]
#_egadslite.EG_matchBodyEdges.restype = c_int
delattr(ego, "matchBodyEdges")

#_egadslite.EG_matchBodyFaces.argtypes = [c_ego, c_ego, c_double, POINTER(c_int), POINTER(POINTER(c_int))]
#_egadslite.EG_matchBodyFaces.restype = c_int
delattr(ego, "matchBodyFaces")

_egadslite.EG_setTessParam.argtypes = [c_ego, c_int, c_double, POINTER(c_double)]
_egadslite.EG_setTessParam.restype = c_int

_egadslite.EG_makeTessGeom.argtypes = [c_ego, POINTER(c_double), POINTER(c_int), POINTER(c_ego)]
_egadslite.EG_makeTessGeom.restype = c_int

_egadslite.EG_getTessGeom.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_double))]
_egadslite.EG_getTessGeom.restype = c_int

_egadslite.EG_makeTessBody.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego)]
_egadslite.EG_makeTessBody.restype = c_int

_egadslite.EG_remakeTess.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_double)]
_egadslite.EG_remakeTess.restype = c_int

_egadslite.EG_finishTess.argtypes = [c_ego, POINTER(c_double)]
_egadslite.EG_finishTess.restype = c_int

#_egadslite.EG_mapTessBody.argtypes = [c_ego, c_ego, POINTER(c_ego)]
#_egadslite.EG_mapTessBody.restype = c_int
delattr(ego, "mapTessBody")

_egadslite.EG_locateTessBody.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(c_double), POINTER(c_int), POINTER(c_double)]
_egadslite.EG_locateTessBody.restype = c_int

_egadslite.EG_getTessEdge.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double))]
_egadslite.EG_getTessEdge.restype = c_int

_egadslite.EG_getTessFace.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_int))]
_egadslite.EG_getTessFace.restype = c_int

_egadslite.EG_getTessLoops.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_int))]
_egadslite.EG_getTessLoops.restype = c_int

_egadslite.EG_getTessQuads.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_int))]
_egadslite.EG_getTessQuads.restype = c_int

_egadslite.EG_makeQuads.argtypes = [c_ego, POINTER(c_double), c_int]
_egadslite.EG_makeQuads.restype = c_int

_egadslite.EG_getQuads.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(POINTER(c_double)), POINTER(POINTER(c_double)), POINTER(POINTER(c_int)), POINTER(POINTER(c_int)), POINTER(c_int)]
_egadslite.EG_getQuads.restype = c_int

_egadslite.EG_getPatch.argtypes = [c_ego, c_int, c_int, POINTER(c_int), POINTER(c_int), POINTER(POINTER(c_int)), POINTER(POINTER(c_int))]
_egadslite.EG_getPatch.restype = c_int

_egadslite.EG_quadTess.argtypes = [c_ego, POINTER(c_ego)]
_egadslite.EG_quadTess.restype = c_int

_egadslite.EG_insertEdgeVerts.argtypes = [c_ego, c_int, c_int, c_int, POINTER(c_double)]
_egadslite.EG_insertEdgeVerts.restype = c_int

_egadslite.EG_deleteEdgeVert.argtypes = [c_ego, c_int, c_int, c_int]
_egadslite.EG_deleteEdgeVert.restype = c_int

_egadslite.EG_moveEdgeVert.argtypes = [c_ego, c_int, c_int, c_double]
_egadslite.EG_moveEdgeVert.restype = c_int

_egadslite.EG_openTessBody.argtypes = [c_ego]
_egadslite.EG_openTessBody.restype = c_int

_egadslite.EG_initTessBody.argtypes = [c_ego, POINTER(c_ego)]
_egadslite.EG_initTessBody.restype = c_int

_egadslite.EG_statusTessBody.argtypes = [c_ego, POINTER(c_ego), POINTER(c_int), POINTER(c_int)]
_egadslite.EG_statusTessBody.restype = c_int

_egadslite.EG_setTessEdge.argtypes = [c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double)]
_egadslite.EG_setTessEdge.restype = c_int

_egadslite.EG_setTessFace.argtypes = [c_ego, c_int, c_int, POINTER(c_double), POINTER(c_double), c_int, POINTER(c_int)]
_egadslite.EG_setTessFace.restype = c_int

_egadslite.EG_localToGlobal.argtypes = [c_ego, c_int, c_int, POINTER(c_int)]
_egadslite.EG_localToGlobal.restype = c_int

_egadslite.EG_getGlobal.argtypes = [c_ego, c_int, POINTER(c_int), POINTER(c_int), POINTER(c_double)]
_egadslite.EG_getGlobal.restype = c_int

#_egadslite.EG_saveTess.argtypes = [c_ego, c_char_p]
#_egadslite.EG_saveTess.restype = c_int
#delattr(ego, "saveTess")

#_egadslite.EG_loadTess.argtypes = [c_ego, c_char_p, POINTER(c_ego)]
#_egadslite.EG_loadTess.restype = c_int
#delattr(ego, "loadTess")

#_egadslite.EG_tessMassProps.argtypes = [c_ego, POINTER(c_double)]
#_egadslite.EG_tessMassProps.restype = c_int
delattr(ego, "tessMassProps")

#_egadslite.EG_fuseSheets.argtypes = [c_ego, c_ego, POINTER(c_ego)]
#_egadslite.EG_fuseSheets.restype = c_int
delattr(ego, "fuseSheets")

#_egadslite.EG_generalBoolean.argtypes = [c_ego, c_ego, c_int, c_double, POINTER(c_ego)]
#_egadslite.EG_generalBoolean.restype = c_int
delattr(ego, "generalBoolean")

#_egadslite.EG_solidBoolean.argtypes = [c_ego, c_ego, c_int, POINTER(c_ego)]
#_egadslite.EG_solidBoolean.restype = c_int

#_egadslite.EG_intersection.argtypes = [c_ego, c_ego, POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(c_ego)]
#_egadslite.EG_intersection.restype = c_int
delattr(ego, "intersection")

#_egadslite.EG_imprintBody.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
#_egadslite.EG_imprintBody.restype = c_int
delattr(ego, "imprintBody")

#_egadslite.EG_filletBody.argtypes = [c_ego, c_int, POINTER(c_ego), c_double, POINTER(c_ego), POINTER(POINTER(c_int))]
#_egadslite.EG_filletBody.restype = c_int
delattr(ego, "filletBody")

#_egadslite.EG_chamferBody.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego), c_double, c_double, POINTER(c_ego), POINTER(POINTER(c_int))]
#_egadslite.EG_chamferBody.restype = c_int
delattr(ego, "chamferBody")

#_egadslite.EG_hollowBody.argtypes = [c_ego, c_int, POINTER(c_ego), c_double, c_int, POINTER(c_ego), POINTER(POINTER(c_int))]
#_egadslite.EG_hollowBody.restype = c_int
delattr(ego, "hollowBody")

#_egadslite.EG_extrude.argtypes = [c_ego, c_double, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_extrude.restype = c_int
delattr(ego, "extrude")

#_egadslite.EG_rotate.argtypes = [c_ego, c_double, POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_rotate.restype = c_int
delattr(ego, "rotate")

#_egadslite.EG_sweep.argtypes = [c_ego, c_ego, c_int, POINTER(c_ego)]
#_egadslite.EG_sweep.restype = c_int
delattr(ego, "sweep")

#_egadslite.EG_loft.argtypes = [c_int, POINTER(c_ego), c_int, POINTER(c_ego)]
#_egadslite.EG_loft.restype = c_int
#del globals()["loft"]

#_egadslite.EG_blend.argtypes = [c_int, POINTER(c_ego), POINTER(c_double), POINTER(c_double), POINTER(c_ego)]
#_egadslite.EG_blend.restype = c_int
del globals()["blend"]

#_egadslite.EG_ruled.argtypes = [c_int, POINTER(c_ego), POINTER(c_ego)]
#_egadslite.EG_ruled.restype = c_int
del globals()["ruled"]

#_egadslite.EG_initEBody.argtypes = [c_ego, c_double, POINTER(c_ego)]
#_egadslite.EG_initEBody.restype = c_int
delattr(ego, "initEBody")

#_egadslite.EG_finishEBody.argtypes = [c_ego]
#_egadslite.EG_finishEBody.restype = c_int
delattr(ego, "finishEBody")

#_egadslite.EG_makeEFace.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
#_egadslite.EG_makeEFace.restype = c_int
delattr(ego, "makeEFace")

#_egadslite.EG_makeAttrEFaces.argtypes = [c_ego, c_char_p, POINTER(c_int), POINTER(POINTER(c_ego))]
#_egadslite.EG_makeAttrEFaces.restype = c_int
delattr(ego, "makeAttrEFaces")

_egadslite.EG_effectiveMap.argtypes = [c_ego, POINTER(c_double), POINTER(c_ego), POINTER(c_double)]
_egadslite.EG_effectiveMap.restype = c_int

_egadslite.EG_effectiveEdgeList.argtypes = [c_ego, POINTER(c_int), POINTER(POINTER(c_ego)), POINTER(POINTER(c_int)), POINTER(POINTER(c_double))]
_egadslite.EG_effectiveEdgeList.restype = c_int

_egadslite.EG_inTriExact.argtypes = [POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)]
_egadslite.EG_inTriExact.restype = c_int

# =============================================================================
# Register the egadslite shared library with common
egads_common._egads = _egadslite
