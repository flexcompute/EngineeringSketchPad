###########################################################################
#                                                                         #
# pyEGADS --- Python version of EGADS API                                 #
#                                                                         #
#                                                                         #
#      Copyright 2011-2025, Massachusetts Institute of Technology         #
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
    raise ImportError("egadslite has already been loaded! egads and egadslite cannot be loaded at the same time.")

__all__ = egads_common._all_

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

_egads.EG_contextCopy.argtypes = [c_ego, c_ego, POINTER(c_ego)]
_egads.EG_contextCopy.restype = c_int

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

_egads.EG_addKnots.argtypes = [c_ego, c_int, c_int, POINTER(c_double), c_int, c_int, POINTER(c_double), POINTER(c_ego)]
_egads.EG_addKnots.restype = c_int

_egads.EG_otherCurve.argtypes = [c_ego, c_ego, c_double, POINTER(c_ego)]
_egads.EG_otherCurve.restype = c_int

_egads.EG_isSame.argtypes = [c_ego, c_ego]
_egads.EG_isSame.restype = c_int

_egads.EG_isoCline.argtypes = [c_ego, c_int, c_double, POINTER(c_ego)]
_egads.EG_isoCline.restype = c_int

_egads.EG_isIsoPCurve.argtypes = [c_ego, POINTER(c_int), POINTER(c_double), POINTER(c_int)]
_egads.EG_isIsoPCurve.restype = c_int

_egads.EG_mergeBSplineCurves.argtypes = [c_ego, c_ego, POINTER(c_ego)]
_egads.EG_mergeBSplineCurves.restype = c_int

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

_egads.EG_makeNmWireBody.argtypes = [ c_int, POINTER(c_ego), c_double, POINTER(c_ego)]
_egads.EG_makeNmWireBody.restype = c_int

_egads.EG_replaceFaces.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_replaceFaces.restype = c_int

_egads.EG_removeNodes.argtypes = [c_ego, c_int, POINTER(c_ego), POINTER(c_ego)]
_egads.EG_removeNodes.restype = c_int

_egads.EG_mapBody.argtypes = [c_ego, c_ego, c_char_p, POINTER(c_ego)]
_egads.EG_mapBody.restype = c_int

_egads.EG_mapBody2.argtypes = [c_ego, c_char_p, c_char_p, POINTER(c_ego)]
_egads.EG_mapBody2.restype = c_int

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

_egads.EG_zeroGeometry_dot.argtypes = [c_ego]
_egads.EG_zeroGeometry_dot.restype = c_int

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
# Register the egads shared library with common
egads_common._egads = _egads

delattr(Context, "importModel")
