using Test
using egadslite, DelimitedFiles
const egads = egadslite

#------------------------------------------------------------------------------
@time @testset " Julia egadslite API: unit testing " begin
    imajor, iminor, OCCrev = egadslite.revision()
    @info "********* EGADSLITE version $imajor.$iminor OCC $OCCrev *********"
    @test typeof(imajor) == Cint
    @test typeof(iminor) == Cint
    @test typeof(OCCrev) == String
    context = egads.Context()
    level   = egads.setOutLevel(context,0)
    data    = read(open("facebody.lite", "r"))
    model   = egads.importModel(context, data)
    @test typeof(model) == egads.Ego

    topo = egads.getTopology(model)
    @test length(topo.children) == 1

    box  = topo.children[1]
    info = egads.getInfo(box)
    @test info.oclass == egads.BODY && info.mtype == egads.SOLIDBODY

    btopo = egads.getTopology(box)
    @test length(btopo.children) == 1

    faces = egads.getBodyTopos(box, egads.FACE)
    @test length(faces) == 6  # egads file is a box

    ftopo = egads.getTopology(faces[1])
    face1 = ftopo.children[1]
    surf  = ftopo.geom
    info  = egads.getGeometry(surf)
    @test info.oclass == egads.SURFACE && info.mtype == egads.PLANE

    rInfo = egads.getRange(surf)
    @test rInfo.periodic == false
    t = [0.5,0.0]
    eval  = egads.evaluate(surf,t)
    @test sum(abs.(eval.X     .- [0.0,0.0,t[1]])) < 1.e-14
    @test sum(abs.(eval.dX[1]  .- [0.0,0.0,1.0])) < 1.e-14
    @test sum(abs.(eval.ddX[1] .- [0.0,0.0,0.0])) < 1.e-14
    inve  = egads.invEvaluate(surf, eval.X)
    @test inve.params == t

    edges = egads.getBodyTopos(box, egads.EDGE)
    rInfo = egads.getRange(edges[1])
    arcL  = egads.arcLength(edges[1], rInfo.range[1], rInfo.range[2])
    @test abs(arcL .- 1.0) < 1.e-10
    info = egads.getBodyTopos(box, egads.FACE, ref = edges[1])

    @test length(info) == 2 && info[1]._ego == faces[1]._ego

    curv = egads.curvature(info[1], t)
    @test iszero(curv.curvature)

    @test egads.indexBodyTopo(box, edges[1]) == 1
    obj = egads.objectBodyTopo(box,egads.EDGE, 1)
    @test obj._ego == edges[1]._ego

    bbox = egads.getBoundingBox(box)
    @test bbox == [0.0 0.0 1.0; 0.0 1.0 1.0]

    @test egads.inTopology(faces[1], [0.0,0.0,1.0])
    @test egads.inFace(faces[1], [0.0,0.0,1.0])
    @test egads.getEdgeUV(faces[1], edges[1],1,0.5) == [0.5,0.0]

    body = egads.getBody(model)
    @test typeof(body) == egads.Ego

    @test egads.getTolerance(faces[1]) == egads.tolerance(faces[1])
    angle = egads.getWindingAngle(edges[1], 0.5)

    tess  = egads.makeTessBody(box, [0.2,0.1,10.0])
    qtess = egads.quadTess(tess)
    egads.remakeTess!(qtess, vcat(faces[1],edges[1]), [0.1,0.1,10.0])

    te   = egads.getTessEdge(qtess,1)
    info = egads.getTessFace(qtess,1)
    lI   = egads.getTessLoops(qtess,1)
    ste  = size(te.ts,1)
    @test lI[1] == ste * 2 + (ste-2) * 2

    egads.makeQuads!(tess, [0.1,0.1,10.0],1)
    @test egads.getTessQuads(tess)[1] == 1
    qeles = egads.getQuads(tess, 1)
    pinfo = egads.getPatch(tess, 1, 1)
    newV = 0.5 * (te.ts[end-2] + te.ts[end-1])
    egads.moveEdgeVert!(qtess,1, ste - 2, newV)

    te2 = egads.getTessEdge(qtess,1)

    @test te2.ts[end-2] == 0.5 * (te.ts[end-2] + te.ts[end-1])
    egads.deleteEdgeVert!(qtess,1,ste-2,1)
    te3 = egads.getTessEdge(qtess,1)

    @test length(te3.ts) == ste-1
    tp = te3.ts[length(te3.ts)-2]

    egads.insertEdgeVerts!(qtess,1,length(te3.ts)-2,newV)
    te3 = egads.getTessEdge(qtess,1)
    @test te3.ts == te2.ts

    egads.openTessBody(qtess)
    nTess = egads.initTessBody(box)
    info  = egads.statusTessBody(nTess)
    for k = 1 : 4
        to = egads.getTessEdge(qtess,k)
        egads.setTessEdge!(nTess, k, to.xyzs, to.ts)
    end
    ten = egads.getTessEdge(nTess,1)
    @test ten == te3

    info = egads.getTessFace(qtess,1)
    egads.setTessFace!(nTess, 1, info.xyz, info.uv, info.tris)
    egads.finishTess!(nTess, [0.2,0.3,40])

    gbl = egads.localToGlobal(nTess, 1,2)
    res = egads.getGlobal(nTess, gbl)
    @test res.ptype > 0
    gbl = egads.localToGlobal(nTess, 1,1)
    res = egads.getGlobal(nTess, gbl)
    @test res.ptype == 0

    @test egads.attributeNum(faces[1]) == 2
    @test egads.attributeRet(faces[1], "faceMap")[1] == 2
    @test egads.attributeRet(faces[1], "faceDouble") == [0.1,0.2]
    nS = egads.attributeNumSeq(box,"faceDouble")
    aS = egads.attributeRetSeq(box,"faceDouble", 1)
    @test nS == 0 && aS == [0.1,0.2]

    @test egads.attributeNum(box) == 3
    info = egads.attributeGet(box, 1)
    @test info.name == "faceDouble" && info.attrVal == [0.1,0.2]
    info = egads.attributeGet(box, 2)
    @test info.name == "faceString" && info.attrVal == "myAttr"
    info = egads.attributeGet(box, 3)
    csys =egads.CSys([0,0, 0, 1,0, 0,0,1,0])
    @test info.name == "faceCSYS" && info.attrVal.val == csys.val && info.attrVal.ortho == csys.ortho

end
