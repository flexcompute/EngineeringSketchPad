function EG_alloc(nbytes)
    ccall((:EG_alloc, C_egadslib), Ptr{Cvoid}, (Csize_t,), nbytes)
end

function EG_calloc(nele, size)
    ccall((:EG_calloc, C_egadslib), Ptr{Cvoid}, (Csize_t, Csize_t), nele, size)
end

function EG_reall(ptr, nbytes)
    ccall((:EG_reall, C_egadslib), Ptr{Cvoid}, (Ptr{Cvoid}, Csize_t), ptr, nbytes)
end

function EG_strdup(str)
    ccall((:EG_strdup, C_egadslib), Ptr{Cchar}, (Ptr{Cchar},), str)
end

function EG_free(pointer)
    ccall((:EG_free, C_egadslib), Cvoid, (Ptr{Cvoid},), pointer)
end

function EG_revision(major, minor, OCCrev)
    ccall((:EG_revision, C_egadslib), Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cchar}}), major, minor, OCCrev)
end

function EG_open(context)
    ccall((:EG_open, C_egadslib), Cint, (Ptr{ego},), context)
end

function EG_loadModel(context, bflg, name, model)
    ccall((:EG_loadModel, C_egadslib), Cint, (ego, Cint, Ptr{Cchar}, Ptr{ego}), context, bflg, name, model)
end

function EG_saveModel(model, name)
    ccall((:EG_saveModel, C_egadslib), Cint, (ego, Ptr{Cchar}), model, name)
end

function EG_exportModel(model, nbytes, stream)
    ccall((:EG_exportModel, C_egadslib), Cint, (ego, Ptr{Csize_t}, Ptr{Ptr{Cchar}}), model, nbytes, stream)
end

function EG_importModel(context, nbytes, stream, model)
    ccall((:EG_importModel, C_egadslib), Cint, (ego, Csize_t, Ptr{Cchar}, Ptr{ego}), context, nbytes, stream, model)
end

function EG_deleteObject(object)
    ccall((:EG_deleteObject, C_egadslib), Cint, (ego,), object)
end

function EG_makeTransform(context, xform, oform)
    ccall((:EG_makeTransform, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), context, xform, oform)
end

function EG_getTransformation(oform, xform)
    ccall((:EG_getTransformation, C_egadslib), Cint, (ego, Ptr{Cdouble}), oform, xform)
end

function EG_getContext(object, context)
    ccall((:EG_getContext, C_egadslib), Cint, (ego, Ptr{ego}), object, context)
end

function EG_setOutLevel(context, outLevel)
    ccall((:EG_setOutLevel, C_egadslib), Cint, (ego, Cint), context, outLevel)
end

function EG_updateThread(context)
    ccall((:EG_updateThread, C_egadslib), Cint, (ego,), context)
end

function EG_getInfo(object, oclass, mtype, topObj, prev, next)
    ccall((:EG_getInfo, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{ego}, Ptr{ego}), object, oclass, mtype, topObj, prev, next)
end

function EG_copyObject(object, oform, copy)
    ccall((:EG_copyObject, C_egadslib), Cint, (ego, Ptr{Cvoid}, Ptr{ego}), object, oform, copy)
end

function EG_flipObject(object, flippedCopy)
    ccall((:EG_flipObject, C_egadslib), Cint, (ego, Ptr{ego}), object, flippedCopy)
end

function EG_close(context)
    ccall((:EG_close, C_egadslib), Cint, (ego,), context)
end

function EG_setUserPointer(context, ptr)
    ccall((:EG_setUserPointer, C_egadslib), Cint, (ego, Ptr{Cvoid}), context, ptr)
end

function EG_getUserPointer(context, ptr)
    ccall((:EG_getUserPointer, C_egadslib), Cint, (ego, Ptr{Ptr{Cvoid}}), context, ptr)
end

function EG_setFullAttrs(context, flag)
    ccall((:EG_setFullAttrs, C_egadslib), Cint, (ego, Cint), context, flag)
end

function EG_attributeAdd(obj, name, type, len, ints, reals, str)
    ccall((:EG_attributeAdd, C_egadslib), Cint, (ego, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}), obj, name, type, len, ints, reals, str)
end

function EG_attributeDel(object, name)
    ccall((:EG_attributeDel, C_egadslib), Cint, (ego, Ptr{Cchar}), object, name)
end

function EG_attributeNum(obj, num)
    ccall((:EG_attributeNum, C_egadslib), Cint, (ego, Ptr{Cint}), obj, num)
end

function EG_attributeGet(obj, index, name, atype, len, ints, reals, str)
    ccall((:EG_attributeGet, C_egadslib), Cint, (ego, Cint, Ptr{Ptr{Cchar}}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, index, name, atype, len, ints, reals, str)
end

function EG_attributeRet(obj, name, atype, len, ints, reals, str)
    ccall((:EG_attributeRet, C_egadslib), Cint, (ego, Ptr{Cchar}, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, name, atype, len, ints, reals, str)
end

function EG_attributeDup(src, dst)
    ccall((:EG_attributeDup, C_egadslib), Cint, (ego, ego), src, dst)
end

function EG_attributeAddSeq(obj, name, type, len, ints, reals, str)
    ccall((:EG_attributeAddSeq, C_egadslib), Cint, (ego, Ptr{Cchar}, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cchar}), obj, name, type, len, ints, reals, str)
end

function EG_attributeNumSeq(obj, name, num)
    ccall((:EG_attributeNumSeq, C_egadslib), Cint, (ego, Ptr{Cchar}, Ptr{Cint}), obj, name, num)
end

function EG_attributeRetSeq(obj, name, index, atype, len, ints, reals, str)
    ccall((:EG_attributeRetSeq, C_egadslib), Cint, (ego, Ptr{Cchar}, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cchar}}), obj, name, index, atype, len, ints, reals, str)
end

function EG_getGeometry(geom, oclass, mtype, refGeom, ivec, rvec)
    ccall((:EG_getGeometry, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Cint}, Ptr{ego}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), geom, oclass, mtype, refGeom, ivec, rvec)
end

function EG_makeGeometry(context, oclass, mtype, refGeom, ivec, rvec, geom)
    ccall((:EG_makeGeometry, C_egadslib), Cint, (ego, Cint, Cint, ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}), context, oclass, mtype, refGeom, ivec, rvec, geom)
end

function EG_getRange(geom, range, periodic)
    ccall((:EG_getRange, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}), geom, range, periodic)
end

function EG_evaluate(geom, param, results)
    ccall((:EG_evaluate, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), geom, param, results)
end

function EG_invEvaluate(geom, xyz, param, results)
    ccall((:EG_invEvaluate, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), geom, xyz, param, results)
end

function EG_invEvaluateGuess(geom, xyz, param, results)
    ccall((:EG_invEvaluateGuess, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), geom, xyz, param, results)
end

function EG_arcLength(geom, t1, t2, alen)
    ccall((:EG_arcLength, C_egadslib), Cint, (ego, Cdouble, Cdouble, Ptr{Cdouble}), geom, t1, t2, alen)
end

function EG_curvature(geom, param, results)
    ccall((:EG_curvature, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), geom, param, results)
end

function EG_approximate(context, maxdeg, tol, sizes, xyzs, bspline)
    ccall((:EG_approximate, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Ptr{ego}), context, maxdeg, tol, sizes, xyzs, bspline)
end

function EG_fitTriangles(context, npts, xyzs, ntris, tris, tric, tol, bspline)
    ccall((:EG_fitTriangles, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Cint, Ptr{Cint}, Ptr{Cint}, Cdouble, Ptr{ego}), context, npts, xyzs, ntris, tris, tric, tol, bspline)
end

function EG_otherCurve(surface, curve, tol, newcurve)
    ccall((:EG_otherCurve, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{ego}), surface, curve, tol, newcurve)
end

function EG_isSame(geom1, geom2)
    ccall((:EG_isSame, C_egadslib), Cint, (ego, ego), geom1, geom2)
end

function EG_isoCline(surface, UV, val, newcurve)
    ccall((:EG_isoCline, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{ego}), surface, UV, val, newcurve)
end

function EG_convertToBSpline(geom, bspline)
    ccall((:EG_convertToBSpline, C_egadslib), Cint, (ego, Ptr{ego}), geom, bspline)
end

function EG_convertToBSplineRange(geom, range, bspline)
    ccall((:EG_convertToBSplineRange, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), geom, range, bspline)
end

function EG_skinning(nCurves, curves, skinning_degree, surface)
    ccall((:EG_skinning, C_egadslib), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}), nCurves, curves, skinning_degree, surface)
end

function EG_tolerance(topo, tol)
    ccall((:EG_tolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}), topo, tol)
end

function EG_getTolerance(topo, tol)
    ccall((:EG_getTolerance, C_egadslib), Cint, (ego, Ptr{Cdouble}), topo, tol)
end

function EG_setTolerance(topo, tol)
    ccall((:EG_setTolerance, C_egadslib), Cint, (ego, Cdouble), topo, tol)
end

function EG_getTopology(topo, geom, oclass, type, limits, nChildren, children, sense)
    ccall((:EG_getTopology, C_egadslib), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}), topo, geom, oclass, type, limits, nChildren, children, sense)
end

function EG_makeTopology(context, geom, oclass, mtype, limits, nChildren, children, senses, topo)
    ccall((:EG_makeTopology, C_egadslib), Cint, (ego, ego, Cint, Cint, Ptr{Cdouble}, Cint, Ptr{ego}, Ptr{Cint}, Ptr{ego}), context, geom, oclass, mtype, limits, nChildren, children, senses, topo)
end

function EG_makeLoop(nedge, edges, geom, toler, result)
    ccall((:EG_makeLoop, C_egadslib), Cint, (Cint, Ptr{ego}, ego, Cdouble, Ptr{ego}), nedge, edges, geom, toler, result)
end

function EG_getArea(object, limits, area)
    ccall((:EG_getArea, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}), object, limits, area)
end

function EG_makeFace(object, mtype, limits, face)
    ccall((:EG_makeFace, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}), object, mtype, limits, face)
end

function EG_getBodyTopos(body, src, oclass, ntopo, topos)
    ccall((:EG_getBodyTopos, C_egadslib), Cint, (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{ego}}), body, src, oclass, ntopo, topos)
end

function EG_indexBodyTopo(body, src)
    ccall((:EG_indexBodyTopo, C_egadslib), Cint, (ego, ego), body, src)
end

function EG_objectBodyTopo(body, oclass, index, obj)
    ccall((:EG_objectBodyTopo, C_egadslib), Cint, (ego, Cint, Cint, Ptr{ego}), body, oclass, index, obj)
end

function EG_inTopology(topo, xyz)
    ccall((:EG_inTopology, C_egadslib), Cint, (ego, Ptr{Cdouble}), topo, xyz)
end

function EG_inFace(face, uv)
    ccall((:EG_inFace, C_egadslib), Cint, (ego, Ptr{Cdouble}), face, uv)
end

function EG_getEdgeUV(face, edge, sense, t, UV)
    ccall((:EG_getEdgeUV, C_egadslib), Cint, (ego, ego, Cint, Cdouble, Ptr{Cdouble}), face, edge, sense, t, UV)
end

function EG_getEdgeUVs(face, edge, sense, nts, ts, UV)
    ccall((:EG_getEdgeUVs, C_egadslib), Cint, (ego, ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), face, edge, sense, nts, ts, UV)
end

function EG_getPCurve(face, edge, sense, mtype, ivec, rvec)
    ccall((:EG_getPCurve, C_egadslib), Cint, (ego, ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), face, edge, sense, mtype, ivec, rvec)
end

function EG_getWindingAngle(edge, t, angle)
    ccall((:EG_getWindingAngle, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}), edge, t, angle)
end

function EG_getBody(object, body)
    ccall((:EG_getBody, C_egadslib), Cint, (ego, Ptr{ego}), object, body)
end

function EG_makeSolidBody(context, stype, rvec, body)
    ccall((:EG_makeSolidBody, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{ego}), context, stype, rvec, body)
end

function EG_getBoundingBox(topo, bbox)
    ccall((:EG_getBoundingBox, C_egadslib), Cint, (ego, Ptr{Cdouble}), topo, bbox)
end

function EG_getMassProperties(topo, result)
    ccall((:EG_getMassProperties, C_egadslib), Cint, (ego, Ptr{Cdouble}), topo, result)
end

function EG_isEquivalent(topo1, topo2)
    ccall((:EG_isEquivalent, C_egadslib), Cint, (ego, ego), topo1, topo2)
end

function EG_sewFaces(nobj, objs, toler, flag, result)
    ccall((:EG_sewFaces, C_egadslib), Cint, (Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}), nobj, objs, toler, flag, result)
end

function EG_replaceFaces(body, nobj, objs, result)
    ccall((:EG_replaceFaces, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), body, nobj, objs, result)
end

function EG_mapBody(sBody, dBody, fAttr, mapBody)
    ccall((:EG_mapBody, C_egadslib), Cint, (ego, ego, Ptr{Cchar}, Ptr{ego}), sBody, dBody, fAttr, mapBody)
end

function EG_matchBodyEdges(body1, body2, toler, nmatch, match)
    ccall((:EG_matchBodyEdges, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}), body1, body2, toler, nmatch, match)
end

function EG_matchBodyFaces(body1, body2, toler, nmatch, match)
    ccall((:EG_matchBodyFaces, C_egadslib), Cint, (ego, ego, Cdouble, Ptr{Cint}, Ptr{Ptr{Cint}}), body1, body2, toler, nmatch, match)
end

function EG_setTessParam(context, iparam, value, oldvalue)
    ccall((:EG_setTessParam, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{Cdouble}), context, iparam, value, oldvalue)
end

function EG_makeTessGeom(obj, params, sizes, tess)
    ccall((:EG_makeTessGeom, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{ego}), obj, params, sizes, tess)
end

function EG_getTessGeom(tess, sizes, xyz)
    ccall((:EG_getTessGeom, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cdouble}}), tess, sizes, xyz)
end

function EG_makeTessBody(object, params, tess)
    ccall((:EG_makeTessBody, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}), object, params, tess)
end

function EG_remakeTess(tess, nobj, objs, params)
    ccall((:EG_remakeTess, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{Cdouble}), tess, nobj, objs, params)
end

function EG_finishTess(tess, params)
    ccall((:EG_finishTess, C_egadslib), Cint, (ego, Ptr{Cdouble}), tess, params)
end

function EG_mapTessBody(tess, body, mapTess)
    ccall((:EG_mapTessBody, C_egadslib), Cint, (ego, ego, Ptr{ego}), tess, body, mapTess)
end

function EG_locateTessBody(tess, npt, ifaces, uv, itri, results)
    ccall((:EG_locateTessBody, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cdouble}), tess, npt, ifaces, uv, itri, results)
end

function EG_getTessEdge(tess, eIndex, len, xyz, t)
    ccall((:EG_getTessEdge, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}), tess, eIndex, len, xyz, t)
end

function EG_getTessFace(tess, fIndex, len, xyz, uv, ptype, pindex, ntri, tris, tric)
    ccall((:EG_getTessFace, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}), tess, fIndex, len, xyz, uv, ptype, pindex, ntri, tris, tric)
end

function EG_getTessLoops(tess, fIndex, nloop, lIndices)
    ccall((:EG_getTessLoops, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cint}}), tess, fIndex, nloop, lIndices)
end

function EG_getTessQuads(tess, nquad, fIndices)
    ccall((:EG_getTessQuads, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{Cint}}), tess, nquad, fIndices)
end

function EG_makeQuads(tess, params, fIndex)
    ccall((:EG_makeQuads, C_egadslib), Cint, (ego, Ptr{Cdouble}, Cint), tess, params, fIndex)
end

function EG_getQuads(tess, fIndex, len, xyz, uv, ptype, pindex, npatch)
    ccall((:EG_getQuads, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}, Ptr{Cint}), tess, fIndex, len, xyz, uv, ptype, pindex, npatch)
end

function EG_getPatch(tess, fIndex, patch, nu, nv, ipts, bounds)
    ccall((:EG_getPatch, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cint}}), tess, fIndex, patch, nu, nv, ipts, bounds)
end

function EG_quadTess(tess, quadTess)
    ccall((:EG_quadTess, C_egadslib), Cint, (ego, Ptr{ego}), tess, quadTess)
end

function EG_insertEdgeVerts(tess, eIndex, vIndex, npts, t)
    ccall((:EG_insertEdgeVerts, C_egadslib), Cint, (ego, Cint, Cint, Cint, Ptr{Cdouble}), tess, eIndex, vIndex, npts, t)
end

function EG_deleteEdgeVert(tess, eIndex, vIndex, dir)
    ccall((:EG_deleteEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cint), tess, eIndex, vIndex, dir)
end

function EG_moveEdgeVert(tess, eIndex, vIndex, t)
    ccall((:EG_moveEdgeVert, C_egadslib), Cint, (ego, Cint, Cint, Cdouble), tess, eIndex, vIndex, t)
end

function EG_openTessBody(tess)
    ccall((:EG_openTessBody, C_egadslib), Cint, (ego,), tess)
end

function EG_initTessBody(object, tess)
    ccall((:EG_initTessBody, C_egadslib), Cint, (ego, Ptr{ego}), object, tess)
end

function EG_statusTessBody(tess, body, state, np)
    ccall((:EG_statusTessBody, C_egadslib), Cint, (ego, Ptr{ego}, Ptr{Cint}, Ptr{Cint}), tess, body, state, np)
end

function EG_setTessEdge(tess, eIndex, len, xyz, t)
    ccall((:EG_setTessEdge, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}), tess, eIndex, len, xyz, t)
end

function EG_setTessFace(tess, fIndex, len, xyz, uv, ntri, tris)
    ccall((:EG_setTessFace, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Cint}), tess, fIndex, len, xyz, uv, ntri, tris)
end

function EG_localToGlobal(tess, index, _local, _global)
    ccall((:EG_localToGlobal, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}), tess, index, _local, _global)
end

function EG_getGlobal(tess, _global, ptype, pindex, xyz)
    ccall((:EG_getGlobal, C_egadslib), Cint, (ego, Cint, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), tess, _global, ptype, pindex, xyz)
end

function EG_saveTess(tess, name)
    ccall((:EG_saveTess, C_egadslib), Cint, (ego, Ptr{Cchar}), tess, name)
end

function EG_loadTess(body, name, tess)
    ccall((:EG_loadTess, C_egadslib), Cint, (ego, Ptr{Cchar}, Ptr{ego}), body, name, tess)
end

function EG_tessMassProps(tess, props)
    ccall((:EG_tessMassProps, C_egadslib), Cint, (ego, Ptr{Cdouble}), tess, props)
end

function EG_fuseSheets(src, tool, sheet)
    ccall((:EG_fuseSheets, C_egadslib), Cint, (ego, ego, Ptr{ego}), src, tool, sheet)
end

function EG_generalBoolean(src, tool, oper, tol, model)
    ccall((:EG_generalBoolean, C_egadslib), Cint, (ego, ego, Cint, Cdouble, Ptr{ego}), src, tool, oper, tol, model)
end

function EG_solidBoolean(src, tool, oper, model)
    ccall((:EG_solidBoolean, C_egadslib), Cint, (ego, ego, Cint, Ptr{ego}), src, tool, oper, model)
end

function EG_intersection(src, tool, nedge, facEdg, model)
    ccall((:EG_intersection, C_egadslib), Cint, (ego, ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{ego}), src, tool, nedge, facEdg, model)
end

function EG_imprintBody(src, nedge, facEdg, result)
    ccall((:EG_imprintBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), src, nedge, facEdg, result)
end

function EG_filletBody(src, nedge, edges, radius, result, facemap)
    ccall((:EG_filletBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}), src, nedge, edges, radius, result, facemap)
end

function EG_chamferBody(src, nedge, edges, faces, dis1, dis2, result, facemap)
    ccall((:EG_chamferBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}, Cdouble, Cdouble, Ptr{ego}, Ptr{Ptr{Cint}}), src, nedge, edges, faces, dis1, dis2, result, facemap)
end

function EG_hollowBody(src, nface, faces, offset, join, result, facemap)
    ccall((:EG_hollowBody, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Cdouble, Cint, Ptr{ego}, Ptr{Ptr{Cint}}), src, nface, faces, offset, join, result, facemap)
end

function EG_extrude(src, dist, dir, result)
    ccall((:EG_extrude, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}), src, dist, dir, result)
end

function EG_rotate(src, angle, axis, result)
    ccall((:EG_rotate, C_egadslib), Cint, (ego, Cdouble, Ptr{Cdouble}, Ptr{ego}), src, angle, axis, result)
end

function EG_sweep(src, spine, mode, result)
    ccall((:EG_sweep, C_egadslib), Cint, (ego, ego, Cint, Ptr{ego}), src, spine, mode, result)
end

function EG_loft(nsec, secs, opt, result)
    ccall((:EG_loft, C_egadslib), Cint, (Cint, Ptr{ego}, Cint, Ptr{ego}), nsec, secs, opt, result)
end

function EG_blend(nsec, secs, rc1, rcN, result)
    ccall((:EG_blend, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{ego}), nsec, secs, rc1, rcN, result)
end

function EG_ruled(nsec, secs, result)
    ccall((:EG_ruled, C_egadslib), Cint, (Cint, Ptr{ego}, Ptr{ego}), nsec, secs, result)
end

function EG_initEBody(tess, angle, ebody)
    ccall((:EG_initEBody, C_egadslib), Cint, (ego, Cdouble, Ptr{ego}), tess, angle, ebody)
end

function EG_finishEBody(EBody)
    ccall((:EG_finishEBody, C_egadslib), Cint, (ego,), EBody)
end

function EG_makeEFace(EBody, nFace, Faces, EFace)
    ccall((:EG_makeEFace, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{ego}), EBody, nFace, Faces, EFace)
end

function EG_makeAttrEFaces(EBody, attrName, nEFace, EFaces)
    ccall((:EG_makeAttrEFaces, C_egadslib), Cint, (ego, Ptr{Cchar}, Ptr{Cint}, Ptr{Ptr{ego}}), EBody, attrName, nEFace, EFaces)
end

function EG_effectiveMap(EObject, eparam, Object, param)
    ccall((:EG_effectiveMap, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{ego}, Ptr{Cdouble}), EObject, eparam, Object, param)
end

function EG_effectiveEdgeList(EEdge, nedge, edges, senses, tstart)
    ccall((:EG_effectiveEdgeList, C_egadslib), Cint, (ego, Ptr{Cint}, Ptr{Ptr{ego}}, Ptr{Ptr{Cint}}, Ptr{Ptr{Cdouble}}), EEdge, nedge, edges, senses, tstart)
end

function EG_effectiveTri(EObj, uv, fIndex, itri, w)
    ccall((:EG_effectiveTri, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cint}, Ptr{Cint}, Ptr{Cdouble}), EObj, uv, fIndex, itri, w)
end

function EG_inTriExact(t1, t2, t3, p, w)
    ccall((:EG_inTriExact, C_egadslib), Cint, (Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), t1, t2, t3, p, w)
end

function EG_makeGeometry_dot(context, oclass, mtype, refGeom, ints, data, data_dot, geom)
    ccall((:EG_makeGeometry_dot, C_egadslib), Cint, (ego, Cint, Cint, ego, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Ptr{egObject}}), context, oclass, mtype, refGeom, ints, data, data_dot, geom)
end

function EG_setGeometry_dot(geom, oclass, mtype, ints, reals, reals_dot)
    ccall((:EG_setGeometry_dot, C_egadslib), Cint, (ego, Cint, Cint, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), geom, oclass, mtype, ints, reals, reals_dot)
end

function EG_getGeometry_dot(geom, reals, reals_dot)
    ccall((:EG_getGeometry_dot, C_egadslib), Cint, (ego, Ptr{Ptr{Cdouble}}, Ptr{Ptr{Cdouble}}), geom, reals, reals_dot)
end

function EG_hasGeometry_dot(geom)
    ccall((:EG_hasGeometry_dot, C_egadslib), Cint, (ego,), geom)
end

function EG_copyGeometry_dot(obj, mat, mat_dot, copy)
    ccall((:EG_copyGeometry_dot, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, ego), obj, mat, mat_dot, copy)
end

function EG_evaluate_dot(geom, params, params_dot, results, results_dot)
    ccall((:EG_evaluate_dot, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), geom, params, params_dot, results, results_dot)
end

function EG_approximate_dot(bspline, mDeg, tol, sizes, xyzs, xyzs_dot)
    ccall((:EG_approximate_dot, C_egadslib), Cint, (ego, Cint, Cdouble, Ptr{Cint}, Ptr{Cdouble}, Ptr{Cdouble}), bspline, mDeg, tol, sizes, xyzs, xyzs_dot)
end

function EG_skinning_dot(surface, nCurves, curves)
    ccall((:EG_skinning_dot, C_egadslib), Cint, (ego, Cint, Ptr{ego}), surface, nCurves, curves)
end

function EG_makeTopology_dot(context, geom, oclass, mtype, limits, limits_dot, nChildren, children, senses, topo)
    ccall((:EG_makeTopology_dot, C_egadslib), Cint, (Ptr{egObject}, Ptr{egObject}, Cint, Cint, Ptr{Cdouble}, Ptr{Cdouble}, Cint, Ptr{Ptr{egObject}}, Ptr{Cint}, Ptr{Ptr{egObject}}), context, geom, oclass, mtype, limits, limits_dot, nChildren, children, senses, topo)
end

function EG_makeSolidBody_dot(body, stype, data, data_dot)
    ccall((:EG_makeSolidBody_dot, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{Cdouble}), body, stype, data, data_dot)
end

function EG_setRange_dot(object, oclass, range, range_dot)
    ccall((:EG_setRange_dot, C_egadslib), Cint, (ego, Cint, Ptr{Cdouble}, Ptr{Cdouble}), object, oclass, range, range_dot)
end

function EG_getRange_dot(geom, range, range_dot, periodic)
    ccall((:EG_getRange_dot, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cint}), geom, range, range_dot, periodic)
end

function EG_makeFace_dot(face, object, limits, limits_dot)
    ccall((:EG_makeFace_dot, C_egadslib), Cint, (ego, ego, Ptr{Cdouble}, Ptr{Cdouble}), face, object, limits, limits_dot)
end

function EG_tessMassProps_dot(tess, xyz_dot, props, props_dot)
    ccall((:EG_tessMassProps_dot, C_egadslib), Cint, (ego, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), tess, xyz_dot, props, props_dot)
end

function EG_extrude_dot(body, src, dist, dist_dot, dir, dir_dot)
    ccall((:EG_extrude_dot, C_egadslib), Cint, (ego, ego, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), body, src, dist, dist_dot, dir, dir_dot)
end

function EG_rotate_dot(body, src, angle, angle_dot, axis, axis_dot)
    ccall((:EG_rotate_dot, C_egadslib), Cint, (ego, ego, Cdouble, Cdouble, Ptr{Cdouble}, Ptr{Cdouble}), body, src, angle, angle_dot, axis, axis_dot)
end

function EG_ruled_dot(body, nSection, sections)
    ccall((:EG_ruled_dot, C_egadslib), Cint, (ego, Cint, Ptr{ego}), body, nSection, sections)
end

function EG_blend_dot(body, nSection, sections, rc1, rc1_dot, rcN, rcN_dot)
    ccall((:EG_blend_dot, C_egadslib), Cint, (ego, Cint, Ptr{ego}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}, Ptr{Cdouble}), body, nSection, sections, rc1, rc1_dot, rcN, rcN_dot)
end

