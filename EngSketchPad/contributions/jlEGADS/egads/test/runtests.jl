@__DIR__###########################################################################
#                                                                         #
#                         unit testing for jlEGADS                        #
#                                                                         #
#                                                                         #
#      Copyright 2011-2024, Massachusetts Institute of Technology         #
#      Licensed under The GNU Lesser General Public License, version 2.1  #
#      See http://www.opensource.org/licenses/lgpl-2.1.php                #
#                                                                         #
#                                                                         #
###########################################################################

using Test

push!(LOAD_PATH, ".")
using egads


#------------------------------------------------------------------------------
function make_wire(context, xyz, class = egads.LOOP)
  d,n = size(xyz)
  @test d == 3
  xyz = Cdouble.(xyz)
  c1n = vcat([1:n;],1)
  nodes = [egads.makeTopology!(context, egads.NODE ; reals=xyz[:,j]) for j = 1:n]

  dir  = hcat([ xyz[:,c1n[j+1]] .- xyz[:,c1n[j]] for j = 1:n]...)
  lines = [egads.makeGeometry(context,egads.CURVE, egads.LINE, vcat(xyz[:,j],dir[:,j]) ) for j = 1:n]
  range = zeros(2,n)
  for j in 1:n
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j]])
    range[1,j] = aux.params
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j+1]])
    range[2,j] = aux.params
  end
  edges = [egads.makeTopology!(context,egads.EDGE ; mtype = egads.TWONODE, geom=lines[j], children = [nodes[c1n[j]], nodes[c1n[j+1]]], reals = range[:,j]) for j = 1:n]

  loop =  egads.makeTopology!(context,egads.LOOP; mtype=egads.CLOSED, children = edges, senses=fill(egads.SFORWARD,n))
  class == egads.LOOP && return loop
  face = egads.makeFace(loop, egads.SFORWARD)
  class == egads.SURFACE && return face
  return egads.makeTopology!(context,egads.BODY; mtype = egads.FACEBODY, children = face)
end

#------------------------------------------------------------------------------
function test_revision(context)
  maxR, minR, occ = egads.revision()
  @info "********* EGADS version $maxR.$minR OCC $occ *********"
  return true
end

#------------------------------------------------------------------------------
function test_makeLoopFace(context)
  xyz   = Cdouble.([[0,0,0] [0,1,0] [1,1,0] [1,0,0]])
  d, n  = size(xyz)
  c1n   = vcat([1:1:n;],1)
  nodes = [egads.makeTopology!(context, egads.NODE ; reals=xyz[:,j]) for j = 1:n]

  dir   = hcat([ xyz[:,c1n[j+1]] .- xyz[:,c1n[j]] for j = 1:n]...)
  lines = [egads.makeGeometry(context,egads.CURVE, egads.LINE, vcat(xyz[:,j],dir[:,j])) for j = 1:n]
  range = zeros(2,n)

  for j in 1:n
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j]])
    range[1,j] = aux.params
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j+1]])
    range[2,j] = aux.params
  end

  edges = [egads.makeTopology!(context,egads.EDGE; mtype = egads.TWONODE, geom=lines[j],
           children = [nodes[c1n[j]], nodes[c1n[j+1]]], reals = range[:,j]) for j = 1:n]

  loop1 = egads.makeTopology!(context,egads.LOOP; mtype=egads.CLOSED,
          children = edges, senses=fill(egads.SFORWARD,n))

  loop2, edges = egads.makeLoop!(edges)

  a = egads.getArea(loop1) ; b = egads.getArea(loop2)
  @test a - b < 1.e-10

  faceL1 = egads.makeFace(loop1, egads.SFORWARD)
  faceL2 = egads.makeFace(loop2, egads.SFORWARD)

  data   = egads.getInfo(faceL1)
  @test data.oclass == egads.FACE && data.mtype == egads.SFORWARD

  a1 = egads.getArea(faceL1) ; a2 = egads.getArea(faceL2)
  @test  a1 - a2 < 1.e-10

  plane = egads.makeGeometry(context, egads.SURFACE, egads.PLANE,[0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0])
  facePlane = egads.makeFace(plane, egads.SFORWARD ; rdata = [0.0,1.0,0.0,1.0])

  data = egads.getInfo(facePlane)
  @test data.oclass == egads.FACE && data.mtype == egads.SFORWARD
  return true
end

#------------------------------------------------------------------------------
function test_makeTransform(context)
  mat1  = [ [1.0,0.0,0.0,4.0] [0.0,1.0,0.0,8.0] [0.0,0.0,1.0,12.0]]
  trans = egads.makeTransform(context,mat1)
  mat2  = egads.getTransform(trans)
  @test mat1 == mat2
  egads.updateThread!(context)
  return true
end

#------------------------------------------------------------------------------
function test_makeGeometry_curves(context)

  xyz    = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  cent   =  [0.0,0.0,0.0]
  xyzAX  = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]

  radius = 0.5
  rM     = 1.5
  rm     = 1.0

  #CIRCLE
  circle = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, vcat(cent, xyzAX[1], xyzAX[2], radius))
  geoC   = egads.getGeometry(circle)
  @test geoC.oclass == egads.CURVE
  @test geoC.mtype  == egads.CIRCLE
  #ELLIPSE
  ellipse = egads.makeGeometry(context, egads.CURVE, egads.ELLIPSE,  vcat(cent, xyzAX[1], xyzAX[2], rM, rm))
  geoE    = egads.getGeometry(ellipse)
  @test geoE.oclass == egads.CURVE
  @test geoE.mtype  == egads.ELLIPSE
  @test geoE.reals[1:3] == xyz[1]

  #PARABOLA
  parabola = egads.makeGeometry(context, egads.CURVE, egads.PARABOLA,  vcat(cent, xyzAX[1], xyzAX[2], rM))
  geoP     = egads.getGeometry(parabola)
  @test geoP.oclass == egads.CURVE
  @test geoP.mtype  == egads.PARABOLA
  @test geoP.reals[end] == rM

  #HIPERBOLA
  hiperbola = egads.makeGeometry(context, egads.CURVE, egads.HYPERBOLA,  vcat(cent, xyzAX[1], xyzAX[2], rm, rM))
  geoH      = egads.getGeometry(hiperbola)
  @test geoH.oclass == egads.CURVE
  @test geoH.mtype  == egads.HYPERBOLA
  @test geoH.reals[1:3] == xyz[1]

  #OFFSET CURVE
  offset  = egads.makeGeometry(context, egads.CURVE, egads.OFFSET,  vcat(xyzAX[1], 2.0), rGeom = circle)

  geoO    = egads.getGeometry(offset)
  @test geoO.oclass == egads.CURVE
  @test geoO.mtype  == egads.OFFSET
  @test geoO.reals  == vcat(xyzAX[1],2.0)
  return true
end

#------------------------------------------------------------------------------
function test_makeGeometry_surfaces(context)

  xyz    = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  cent   =  [0.0,0.0,0.0]
  xyzAX  = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  radius = 0.5
  rM     = 1.5
  rm     = 1.0
  angle  = 30.0 * 3.14/180.
  #SPHERE
  sphere = egads.makeGeometry(context, egads.SURFACE, egads.SPHERICAL, vcat(cent, xyzAX[1], xyzAX[2], radius))
  geoS   = egads.getGeometry(sphere)
  @test geoS.oclass == egads.SURFACE
  @test geoS.mtype  == egads.SPHERICAL

  #CONICAL
  conical = egads.makeGeometry(context, egads.SURFACE, egads.CONICAL,  vcat(cent, xyzAX[1], xyzAX[2], xyzAX[3], angle, rM))
  geoC    = egads.getGeometry(conical)
  @test geoC.oclass == egads.SURFACE
  @test geoC.mtype  == egads.CONICAL
  @test geoC.reals[4:6] == xyzAX[1]

  #CYLINDER
  cylinder = egads.makeGeometry(context, egads.SURFACE, egads.CYLINDRICAL,  vcat(cent, xyzAX[1], xyzAX[2], xyzAX[3], rM))
  geoC     = egads.getGeometry(cylinder)
  @test geoC.oclass == egads.SURFACE
  @test geoC.mtype  == egads.CYLINDRICAL

  #TOROIDAL
  toroidal = egads.makeGeometry(context, egads.SURFACE, egads.TOROIDAL,  vcat(cent, xyzAX[1], xyzAX[2], xyzAX[3], rm, rM))
  geoT     = egads.getGeometry(toroidal)
  @test geoT.oclass == egads.SURFACE
  @test geoT.mtype  == egads.TOROIDAL
  @test geoT.reals[4:12] == vcat(xyzAX...)

  #SURFACE REVOLUTION
  line       = egads.makeGeometry(context, egads.CURVE, egads.LINE, vcat(xyz[1], xyz[2]))
  revolution = egads.makeGeometry(context, egads.SURFACE, egads.REVOLUTION,  vcat(cent, xyzAX[3]); rGeom = line)
  geoR       = egads.getGeometry(revolution)
  @test geoR.oclass == egads.SURFACE
  @test geoR.mtype  == egads.REVOLUTION
  @test geoR.reals[1:3] == cent

  #EXTRUSION SURFACE
  circle     = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, vcat(cent, xyzAX[1], xyzAX[2], radius))
  extrusion  = egads.makeGeometry(context, egads.SURFACE, egads.EXTRUSION, xyzAX[3]; rGeom = circle)
  geoE       = egads.getGeometry(extrusion)
  @test geoE.oclass == egads.SURFACE
  @test geoE.mtype  == egads.EXTRUSION
  @test geoE.reals  == xyzAX[3]

  #OFFSET SURFACE
  offset  = egads.makeGeometry(context, egads.SURFACE, egads.OFFSET, 2.0; rGeom = cylinder)
  geoO    = egads.getGeometry(offset)
  @test geoO.oclass == egads.SURFACE
  @test geoO.mtype  == egads.OFFSET
  @test geoO.reals  == 2.0
  return true
end

#------------------------------------------------------------------------------
function test_makeBezier(context)
  #BEZIER SURFACE
  nCPu = 4
  nCPv = 4

  rdata  = [[0.0,0.0,0.0] [1.0,0.0,0.1] [1.5,1.0,0.7] [0.25,0.75,0.6]  #=
            =# [0.0,0.0,1.0] [1.0,0.0,1.1] [1.5,1.0,1.7] [0.25,0.75,1.6] #=
            =# [0.0,0.0,2.0] [1.0,0.0,2.1] [1.5,1.0,2.7] [0.25,0.75,2.6] #=
            =# [0.0,0.0,3.0] [1.0,0.0,3.1] [1.5,1.0,3.7] [0.25,0.75,3.6]]

  idata  = [0, nCPu-1, nCPu, nCPv-1,nCPv]
  beS    = egads.makeGeometry(context, egads.SURFACE, egads.BEZIER,rdata ; ints = idata)
  geoB   = egads.getGeometry(beS)
  @test geoB.oclass == egads.SURFACE
  @test geoB.mtype  == egads.BEZIER
  @test geoB.ints   == idata

  # make Nodes
  rn    = hcat(rdata[:,1], rdata[:,4], rdata[:,16] , rdata[:,12] )
  nodes = [egads.makeTopology!(context, egads.NODE ; reals=rn[:,j]) for j = 1:4]
  #BEZIER CURVE
  idata = [0, nCPu-1,nCPu]
  idx   = [[1,2,3,4], [4,8,12,16], [16,15,14,13],[13,9,5,1]]
  beC   = [egads.makeGeometry(context, egads.CURVE, egads.BEZIER, rdata[:,idx[j]] ; ints = idata) for j =1:4]
  geoB   = egads.getGeometry(beC[1])
  @test geoB.oclass == egads.CURVE
  @test geoB.mtype  == egads.BEZIER
  @test geoB.ints   == idata
  c1n   = vcat([1:4;],1)
  edges = [egads.makeTopology!(context, egads.EDGE; mtype = egads.TWONODE, geom = beC[j], children=vcat(nodes[c1n[j]],nodes[c1n[j+1]]), reals=[0.0,1.0]) for j = 1:4]
  pcurves = [egads.otherCurve(beS, e) for e in edges]
  edges   = vcat(edges, pcurves)

  loopF   = egads.makeTopology!(context, egads.LOOP ; mtype = egads.CLOSED, geom=beS, children= edges, senses = fill(egads.SFORWARD,4) )
  @test egads.getArea(loopF) > 0.0

  faceF = egads.makeTopology!(context, egads.FACE; mtype = egads.SFORWARD, geom=beS, children= loopF)
  # making a face with the loop reversed
  piv = [3,2,1,4,7,6,5,8]


  loopR = egads.makeTopology!(context, egads.LOOP ; mtype = egads.CLOSED, geom=beS, children=edges[piv], senses=fill(egads.SREVERSE,4))

  @test egads.getArea(loopR) < 0.0
  faceR = egads.makeTopology!(context, egads.FACE; mtype = egads.SREVERSE, geom=beS, children= loopR)
  @test egads.getArea(faceF) - egads.getArea(faceR) < 1.e-10

  plane = egads.makeGeometry(context, egads.SURFACE , egads.PLANE, [0.0,0.0,0.0,1.0,0.0,0.0, 0.0,1.0,0.0])

  planeFace = egads.makeFace(plane, egads.SREVERSE ; rdata = [0.0 1.0 0.0 1.0])
  info      = egads.getInfo(planeFace)

  @test info.oclass == egads.FACE && info.mtype == egads.SREVERSE

  a1 = egads.getArea(planeFace)
  a2 = egads.getArea(faceR)
  @test a1 - a2 < 1.e-10


  flip  = egads.flipObject(planeFace)

  a1   = egads.getArea(planeFace)

  @test a1 - a2 < 1.e-10

  # Outer Loops
  # RIGHT HANDED
  rOut = make_wire(context, [[0.0,0,0] [1.0,0,0] [1.0,1.0,0] [0,1.0,0]], egads.LOOP)

  flip_rOut = egads.flipObject(rOut)

  #LEFT HANDED
  lOut = make_wire(context, [[0,0,0] [0,1,0] [1,1,0] [1,0,0]], egads.LOOP)
  flip_lOut = egads.flipObject(lOut)
  # Inner Loops
  rIn = make_wire(context, [[0.25,0.25,0] [0.75,0.25,0] [0.75,0.75,0] [0.25,0.75,0]], egads.LOOP)


  flip_rIn = egads.flipObject(rIn)
  lIn      = make_wire(context, [[0.25,0.25,0] [0.25,0.75,0] [0.75,0.75,0] [0.75,0.25,0]], egads.LOOP)
  flip_lIn = egads.flipObject(lIn)


  # FORWARDegads.FACES
  fFace1 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SFORWARD, geom=plane, children=[rOut, lIn], senses=[egads.SOUTER, egads.SINNER])
  fFace2 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SFORWARD, geom=plane, children=[rOut, flip_rIn], senses=[egads.SOUTER, egads.SINNER])
  fFace3 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SFORWARD, geom=plane, children=[flip_lOut, lIn], senses=[egads.SOUTER, egads.SINNER])

  # REVERSE FACES
  rFace1 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SREVERSE, geom=plane, children=[lOut, rIn], senses=[egads.SOUTER, egads.SINNER])
  rFace2 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SREVERSE, geom=plane, children=[lOut, flip_lIn], senses=[egads.SOUTER, egads.SINNER])
  rFace3 = egads.makeTopology!(context, egads.FACE ; mtype = egads.SREVERSE, geom=plane, children=[flip_rOut, lIn], senses=[egads.SOUTER, egads.SINNER])
  return true
end

#------------------------------------------------------------------------------
function test_convertToBSpline(context)

  dataE = [1.0,0.0,0.0,1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4.0, 0.5]
  geoE  = egads.makeGeometry(context, egads.CURVE, egads.ELLIPSE, dataE)

  nCP   = 20
  rInfo = egads.getRange(geoE)

  geoSpline = egads.convertToBSplineRange(geoE, rInfo.range)

  geoInfo = egads.getGeometry(geoSpline)
  @test geoInfo.oclass == egads.CURVE && geoInfo.mtype == egads.BSPLINE
  h      = (rInfo.range[2] - rInfo.range[1]) / (nCP + 1)
  jeval  = [egads.evaluate(geoE, rInfo.range[1] + h * j) for j in 1:nCP]
  data3D = hcat(getfield.(jeval,:X)...)

  # Fit spline through data points
  bspline = egads.approximate(context, [nCP,0], data3D)

  # LOOK BACKWARDS STUFF
  geoInfo = egads.getGeometry(bspline)

  @test geoInfo.oclass == egads.CURVE && geoInfo.mtype == egads.BSPLINE

  for j in 1:nCP
    invJ  = egads.invEvaluate(bspline, data3D[:,j])
    e_xyz = sum((invJ.result .- data3D[:,j]).^2)
    @test e_xyz < 1.e-5
  end
  dupSpline = egads.convertToBSpline(bspline)
  stat      = egads.isSame(bspline, dupSpline)
  return true
end

#------------------------------------------------------------------------------
function test_mergeBSplineCurves(context)

  pts = zeros(3,2)
  pts[:,1] = [0.0, 0.0, 0.0]
  pts[:,2] = [1.0, 0.0, 0.0]

  bspline0 = egads.approximate(context, [2, 0], pts)

  pts[:,1] = [1.0, 0.0, 0.0]
  pts[:,2] = [2.0, 0.0, 0.0]

  bspline1 = egads.approximate(context, [2, 0], pts)

  bspline3 = egads.mergeBSplineCurves(bspline0, bspline1)

  data = egads.getGeometry(bspline3)

  @test data.oclass == egads.CURVE
  @test data.mtype == egads.BSPLINE
end

#------------------------------------------------------------------------------
function test_skinning(context)

  npts   = 20
  radius = 1.0
  points = zeros(3,npts)
  for j = 0:npts-1
    points[:,j+1] = [radius .* cos(pi * j / npts), radius .* sin(pi * j / npts), 0.0]
  end

  curve1       = egads.approximate(context, [npts, 0], points)
  points[3,:] .= 1.0
  curve2       = egads.approximate(context, [npts, 0], points)
  bspline = egads.skinning([curve1,curve2], 1)

  egadsfile = joinpath(@__DIR__, "skinning.egads")

  range, per = egads.getRange(bspline)
  face       = egads.makeFace(bspline, egads.SFORWARD ; rdata =  range)
  splineBody = egads.makeTopology!(context, egads.BODY;mtype=egads.FACEBODY,children=face)

  faces      = egads.getBodyTopos(splineBody,egads.FACE)

  body2 = egads.getBody(faces[1])

  stat = egads.isEquivalent(body2, splineBody)
  egads.saveModel!(body2,egadsfile,true)
  rm(egadsfile)
  return true
end

#------------------------------------------------------------------------------
function test_attribution(context)

  egads.setFullAttrs!(context, 1)

  box1 = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context,egads.BOX, [1.0,0.0,0.0,1.0,1.0,1.0])

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  egads.attributeAdd!(faces1[2], "faceMap", 2)
  egads.attributeAdd!(faces2[1], "faceMap", 1)
 
  attrVal = egads.attributeRet(faces1[2], "faceMap")
  nAttr   = egads.attributeNum(faces1[2])

  @test attrVal[1] == 2 && nAttr == 1

  attrVal = egads.attributeRet(faces2[1], "faceMap")
  nAttr   = egads.attributeNum(faces2[1])

  @test attrVal[1] == 1 && nAttr == 1

  egads.attributeAdd!(faces2[1],"faceDouble", [0.1,0.2])

  attrVal = egads.attributeRet(faces2[1], "faceDouble")
  nAttr   = egads.attributeNum(faces2[1])
  @test attrVal == [0.1,0.2]
  @test nAttr == 2

  egads.attributeAddSeq!(box1,"faceDouble", [0.1,0.2])
  egads.attributeAddSeq!(box1,"faceDouble", [0.3,0.4])

  attrVal = egads.attributeRetSeq(box1, "faceDouble", 1)
  nAttr   = egads.attributeNumSeq(box1,"faceDouble")

  @test attrVal == [0.1,0.2] && nAttr == 2

  egads.attributeAdd!(box1, "faceString", "myAttr")

  nAttr   = egads.attributeNum(box1)
  attrVal = egads.attributeRet(box1, "faceString")

  @test attrVal == "myAttr" && nAttr == 3

  info = egads.attributeGet(box1, 3)
  @test info.name == "faceString" && info.attrVal == attrVal

  egads.attributeDup!(box1,box2)
  nAttr2 = egads.attributeNum(box2)
  @test nAttr2 == nAttr


  egads.attributeDel!(box1; name = "faceString")
  @test egads.attributeNum(box1) == 2
  egads.attributeDel!(box1)
  @test  egads.attributeNum(box1) == 0
  # csystem

  attrCsys1 =  egads.CSys([0,0, 0, 1,0, 0,0,1,0])
  egads.attributeAdd!(box1,"faceCSYS", attrCsys1)

  data = egads.attributeRet(box1, "faceCSYS")
  @test data.val == attrCsys1.val && typeof(data) == egads.CSys

  data = egads.attributeRetSeq(box1, "faceCSYS", 1)
  @test data.ortho == attrCsys1.ortho && typeof(data) == egads.CSys

  name, data = egads.attributeGet(box1, 1)
  @test name == "faceCSYS" && typeof(data) == egads.CSys
  return true
end


function test_attribution_for_lite(context)

  egads.setFullAttrs!(context, 1)
  box   = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  faces = egads.getBodyTopos(box, egads.FACE)

  egads.attributeAdd!(faces[1], "faceMap", 2)
  egads.attributeAdd!(faces[1],"faceDouble", [0.1,0.2])
  egads.attributeAddSeq!(box,"faceDouble", [0.1,0.2])
  egads.attributeAdd!(box, "faceString", "myAttr")
  attrCsys1 =  egads.CSys([0,0, 0, 1,0, 0,0,1,0])
  egads.attributeAdd!(box,"faceCSYS", attrCsys1)
  model = egads.makeTopology!(context, egads.MODEL; children = box)
  stream = egads.exportModel(model)
  @test typeof(stream) == Vector{UInt8}
  fp = joinpath(@__DIR__, "../../egadslite/test/facebody.lite") |> normpath
  open(fp, "w") do io
           write(io,stream) ; end
  fp = joinpath(@__DIR__, "../../egadslite/test/facebody.egads") |> normpath
  egads.saveModel!(model,fp,true)
  return true
end

#------------------------------------------------------------------------------
function test_tolerance(context)

  node1 = egads.makeTopology!(context,egads.NODE; reals= [0.0,0.0,0.0])
  node2 = egads.makeTopology!(context,egads.NODE; reals= [1.0,0.0,0.0])
  ldata = [ 0.0,0.0,0.0,1.0,0.0,0.0]
  line  = egads.makeGeometry(context, egads.CURVE,egads.LINE, ldata)
  edge  = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=line, children=[node1,node2], reals = [0.0,1.0])
  # tolerance
  tol1 = egads.tolerance(edge)
  tol2 = egads.getTolerance(edge)
  return true
end

#------------------------------------------------------------------------------
function test_make_export_IO_model(context)
  plane    = egads.makeGeometry(context, egads.SURFACE , egads.PLANE, [0.0,0.0,0.0,1.0,0.0,0.0, 0.0,1.0,0.0])
  pFace    = egads.makeFace(plane, egads.SREVERSE ; rdata = [0.0,1.0,0.0,1.0])
  faceBody = egads.makeTopology!(context, egads.BODY; mtype = egads.FACEBODY , children= pFace)

  egadsfile = joinpath(@__DIR__,"faceBody.egads")

  model    = egads.makeTopology!(context, egads.MODEL; children = faceBody)
  egads.saveModel!(model, egadsfile, true)
  mSB   = egads.loadModel(context, 0, egadsfile)
  mTopo = egads.getTopology(mSB)

  info = egads.getInfo(mTopo.children[1])
  @test info.oclass == egads.BODY && info.mtype == egads.FACEBODY
  face = egads.objectBodyTopo(mTopo.children[1],egads.FACE,1)
  @test egads.isSame(face, pFace)

  bbox   = egads.getBoundingBox(mSB)
  rm(egadsfile)

  finalize(model)
  finalize(faceBody)
  finalize(pFace)
  finalize(plane)
  finalize(mSB)
  return true
end

#------------------------------------------------------------------------------
function test_geoProperties(context)
  #####   CREATE A LINE
  dataL  = [1.0,0.0,0.0,1.0, 1.0, 0.0]
  geoL   = egads.makeGeometry(context, egads.CURVE,egads.LINE, dataL)
  radius = 0.5
  dataC  = [0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,radius]
  dupL   = egads.getGeometry(geoL)
  info   = egads.getInfo(geoL)

  @test info.oclass == dupL.oclass && info.mtype == dupL.mtype && dupL.oclass == egads.CURVE && dupL.mtype == egads.LINE

  #####   CREATE A CIRCLE
  geoC = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, dataC)
  refG = egads.getGeometry(geoC)
  @test refG.oclass == egads.CURVE && refG.mtype == egads.CIRCLE

  @test sum((dataC .- refG.reals).^2) < 1.e-14

  rangeInfo,periodic = egads.getRange(geoC)

  arc = egads.arcLength(geoC,rangeInfo[1], rangeInfo[2])
  xyz = egads.evaluate(geoC,rangeInfo[1])

  @test abs(2.0 * pi * radius - arc) < 1.e-14

  kappa = egads.curvature(geoC, rangeInfo)
  @test abs(kappa.curvature - 1 / radius) < 1.e-15

  dataC = [0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0,radius]
  cylinder = egads.makeGeometry(context, egads.SURFACE, egads.CYLINDRICAL, dataC)

  rInfo = egads.getRange(cylinder)


  mp    = [0.5 * sum(rInfo.range[1,:]), 0.5 * sum(rInfo.range[2,:])]
  kappa = egads.curvature(cylinder, mp)
  @test (abs(kappa.curvature[1] + 1.0 / radius) < 1.e-15 || abs(kappa.curvature[2] + 1.0 / radius) < 1.e-15)

  isoU     = egads.isoCline(cylinder, 0, mp[1]) ;  isoV = egads.isoCline(cylinder, 1, mp[2])
  isoInfoU = egads.getInfo(isoU) ;  isoInfoV = egads.getInfo(isoV)
  @test isoInfoU.oclass == egads.CURVE && isoInfoU.mtype == egads.LINE && isoInfoV.oclass == egads.CURVE && isoInfoV.mtype == egads.CIRCLE

  sphere = egads.makeSolidBody(context, egads.SPHERE,[0.0,0.0,0.0,radius])

  topF   = egads.getBodyTopos(sphere, egads.FACE)

  infoR  = egads.getRange(topF[1])

  area = egads.getArea(topF[1]; limits = infoR.range)
  @test abs(2.0 *pi * radius^2 - area) < 1.e-13

  massP = egads.getMassProperties(sphere)

  @test abs(massP.area - 4 * pi * radius^2) < 1.e-13 && abs(massP.volume - 4 / 3 * pi * radius^3) < 1.e-13 && sum(massP.CG.^2) < 1.e-13

  edges = egads.getBodyTopos(sphere, egads.EDGE)

  rInfo = egads.getRange(edges[2])

  t     = (rInfo.range[1] + rInfo.range[2]) * 0.5
  @test egads.getWindingAngle(edges[2], t) == 180.0
  return true
end

#------------------------------------------------------------------------------
function test_evaluate(context)
  wireBody = make_wire(context, [[0,0,0]  [1,0,0] [1,1,0] [0,1,0]], egads.BODY)
  info     = egads.getInfo(wireBody)
  nodes    = egads.getBodyTopos(wireBody,egads.NODE)
  @test egads.evaluate(nodes[1]) == [0.0,0.0,0.0]
  edges    = egads.getBodyTopos(wireBody,egads.EDGE)

  rE    = egads.getRange(edges[1])
  eval  = egads.evaluate(edges[1],0.5 * sum(rE.range))
  inve  = egads.invEvaluate(edges[1], [0.5, 0., 0.])

  @test eval.X      == [0.5,0.0,0.0]
  @test eval.dX     == [1.0,0.0,0.0]
  @test eval.ddX    == [0.0,0.0,0.0]
  @test inve.result == [0.5,0.0,0.0]

  faces = egads.getBodyTopos(wireBody,egads.FACE)
  rF    = egads.getRange(faces[1])
  uv    = 0.5 .* [sum(rF.range[1,:]) ,sum(rF.range[2,:]) ]

  eval  = egads.evaluate(faces[1], uv)

  inve  = egads.invEvaluate(faces[1], [0.5, 0.5, 0.])
  @test isapprox(inve.result,[0.5,0.5,0.0],atol = 1.e-10)
  @test isapprox(eval.X,[0.5,0.5,0.0],atol = 1.e-10)

  dv = cos(pi*0.25)
  @test isapprox(eval.dX[1],[-dv,-dv,0], atol = 1.e-10) && isapprox(eval.dX[2],[dv,-dv,0], atol = 1.e-10) &&
        isapprox(eval.ddX[1],[0,0,0], atol = 1.e-10)    && isapprox(eval.ddX[2],[0,0,0], atol = 1.e-10) &&
        isapprox(eval.ddX[3],[0,0,0], atol = 1.e-10)    && isapprox(inve.result,[0.5,0.5,0.0],atol = 1.e-10) &&
        isapprox(inve.params,uv,atol = 1.e-10)
  return true
end

#------------------------------------------------------------------------------
function test_otherCurve(context)
  box   = egads.makeSolidBody(context, egads.BOX, [1,0,0,-1,1,1,2])
  edges = egads.getBodyTopos(box, egads.EDGE)
  idx   = [egads.indexBodyTopo(box,edges[i])  for i = 1:length(edges)]
  @test all(idx .> 0)
  edge1 = egads.objectBodyTopo(box,egads.EDGE,1)
  @test egads.isSame(edge1, edges[1])

  faces = egads.getBodyTopos(box, egads.FACE ; ref = edges[1])
  pcurv = egads.otherCurve(faces[1], edges[1])
  data = egads.getInfo(pcurv)
  @test data.oclass == egads.PCURVE && data.mtype == egads.LINE
  return true
end

#------------------------------------------------------------------------------
function test_isoCline(context)

  box0 = egads.makeSolidBody(context, egads.BOX, [1,0,0,-1,1,1,2])

  faces = egads.getBodyTopos(box0, egads.FACE)
  uvrange, periodic = egads.getRange(faces[1])

  mTopo = egads.getTopology(faces[1])

  ucurv = egads.isoCline(mTopo.geom, egads.UISO, (uvrange[1]+uvrange[2])/2)
  data = egads.getInfo(ucurv)
  @test data.oclass == egads.CURVE && data.mtype == egads.LINE

  vcurv = egads.isoCline(mTopo.geom, egads.VISO, (uvrange[3]+uvrange[4])/2)
  data = egads.getInfo(vcurv)
  @test data.oclass == egads.CURVE && data.mtype == egads.LINE

end

#------------------------------------------------------------------------------
function test_isIsoPCurve(context)

  box0 = egads.makeSolidBody(context, egads.BOX, [1,0,0,-1,1,1,2])

  faces = egads.getBodyTopos(box0, egads.FACE)
  uvrange, periodic = egads.getRange(faces[1])

  mTopo = egads.getTopology(faces[1])

  ucurv = egads.isoCline(mTopo.geom, egads.UISO, (uvrange[1]+uvrange[2])/2)
  pcurve = egads.otherCurve(faces[1], ucurv)
  isIso, iUV, value, fwd = egads.isIsoPCurve(pcurve)
  @test isIso == true
  @test iUV   == egads.UISO
  @test value == (uvrange[1]+uvrange[2])/2
  @test fwd   == egads.SFORWARD

  vcurv = egads.isoCline(mTopo.geom, egads.VISO, (uvrange[3]+uvrange[4])/2)
  pcurve = egads.otherCurve(faces[1], vcurv)
  isIso, iUV, value, fwd = egads.isIsoPCurve(pcurve)
  @test isIso == true
  @test iUV   == egads.VISO
  @test value == (uvrange[3]+uvrange[4])/2
  @test fwd   == egads.SFORWARD
end

#------------------------------------------------------------------------------
function test_getEdgeUV(context)

  wireBody = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.BODY)
  faces    = egads.getBodyTopos(wireBody,egads.FACE)
  edges    = egads.getBodyTopos(wireBody,egads.EDGE)
  rInfo    = egads.getRange(edges[1])
  t        = 0.5 * sum(rInfo.range)
  uvs      = egads.getEdgeUV(faces[1], edges[1], 1,t)

  xyz1  = egads.evaluate(faces[1], uvs)
  xyz2  = egads.evaluate(edges[1], t)
  @test sum((xyz1.X .- xyz2.X).^2) < 1.e-14
  return true
end

#------------------------------------------------------------------------------
function test_removeNodes(context)

  # construction data
  x0 = [0.0,0.0,0.0]
  x1 = [1.0,0.0,0.0]
  x2 = [1.0,1.0,0.0]
  x3 = [0.0,1.0,0.0]
  x4 = [0.0,0.5,0.0] # redundant node to be removed

  # make Nodes
  nodes = Array{egads.Ego}(undef, 5)
  
  nodes[1] = egads.makeTopology!(context, egads.NODE; reals=x0)
  nodes[2] = egads.makeTopology!(context, egads.NODE; reals=x1)
  nodes[3] = egads.makeTopology!(context, egads.NODE; reals=x2)
  nodes[4] = egads.makeTopology!(context, egads.NODE; reals=x3)
  nodes[5] = egads.makeTopology!(context, egads.NODE; reals=x4)

  # make lines and edges
  lines = Array{egads.Ego}(undef, 4)
  edges = Array{egads.Ego}(undef, 5)
  tdata = [0.0, 1.0]

  # Line data (point and direction)
  rdata = fill(0.0, 6)
  rdata[1] = x0[1]
  rdata[2] = x0[2]
  rdata[3] = x0[3]
  rdata[4] = x1[1] - x0[1]
  rdata[5] = x1[2] - x0[2]
  rdata[6] = x1[3] - x0[3]

  lines[1] = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)
  edges[1] = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=lines[1], children=nodes[1:2], reals=tdata)

  rdata[1] = x1[1]
  rdata[2] = x1[2]
  rdata[3] = x1[3]
  rdata[4] = x2[1] - x1[1]
  rdata[5] = x2[2] - x1[2]
  rdata[6] = x2[3] - x1[3]

  lines[2] = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)
  edges[2] = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=lines[2], children=nodes[2:3], reals=tdata)

  rdata[1] = x2[1]
  rdata[2] = x2[2]
  rdata[3] = x2[3]
  rdata[4] = x3[1] - x2[1]
  rdata[5] = x3[2] - x2[2]
  rdata[6] = x3[3] - x2[3]

  lines[3] = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)
  edges[3] = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=lines[3], children=nodes[3:4], reals=tdata)

  rdata[1] = x3[1]
  rdata[2] = x3[2]
  rdata[3] = x3[3]
  rdata[4] = x0[1] - x3[1]
  rdata[5] = x0[2] - x3[2]
  rdata[6] = x0[3] - x3[3]

  lines[4] = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)

  edges[4] = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=lines[4], children=[nodes[4], nodes[5]], reals=[0.0,0.5])
  edges[5] = egads.makeTopology!(context, egads.EDGE; mtype=egads.TWONODE, geom=lines[4], children=[nodes[5], nodes[1]], reals=[0.5,1.0])

  # make the loop
  sens = fill(egads.SFORWARD, 5)
  loop = egads.makeTopology!(context, egads.LOOP; mtype=egads.CLOSED, children=edges, senses=sens)

  # create Plane
  rdata = fill(0.0, 9)
  rdata[1] = 0 # center
  rdata[2] = 0
  rdata[3] = 0
  rdata[4] = 1 # x-axis
  rdata[5] = 0
  rdata[6] = 0
  rdata[7] = 0 # y-axis
  rdata[8] = 1
  rdata[9] = 0

  plane = egads.makeGeometry(context, egads.SURFACE, egads.PLANE, rdata)
  face = egads.makeTopology!(context, egads.FACE; mtype=egads.SFORWARD, geom=plane, children=[loop])

  body = egads.makeTopology!(context, egads.BODY; mtype=egads.FACEBODY, children=[face])
  
  nodes = egads.getBodyTopos(body, egads.NODE)

  newbody = egads.removeNodes(body, [nodes[5]])
  
  newnodes = egads.getBodyTopos(newbody, egads.NODE)
  newedges = egads.getBodyTopos(newbody, egads.EDGE)
  @test length(newnodes) == 4
  @test length(newedges) == 4
end

#------------------------------------------------------------------------------
function test_match_sew_replace_topology(context)

  box1 = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context,egads.BOX, [1.0,0.0,0.0,1.0,1.0,1.0])

  matchB = egads.matchBodyEdges(box1,box2)
  @test matchB == [[5, 1]  [6, 2]  [7, 3] [8, 4]]

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  matchB = egads.matchBodyFaces(box1,box2)
  @test matchB == [2,1]

  piv1 = deleteat!([1:1:length(faces1);], matchB[1,:])
  piv2 = deleteat!([1:1:length(faces2);], matchB[2,:])
  faces12 = vcat(faces1[piv1],faces2[matchB[2,:]])

  sew     = egads.sewFaces(faces12)

  bodyT = egads.getTopology(sew)

  topo3 = egads.getBodyTopos(bodyT.children[1], egads.FACE)

  box3 = egads.getBody(topo3[1])

  stat = egads.isSame(box3, bodyT.children[1])

  sheet  = egads.replaceFaces(box1,[ [faces1[1], nothing] [faces1[2],nothing] ])
  rFaces = egads.getBodyTopos(sheet, egads.FACE)
  @test length(rFaces) == 4

  radius = 3.0
  node1  = egads.makeTopology!(context, egads.NODE;reals=[ radius,0.0,0.0])
  node2  = egads.makeTopology!(context, egads.NODE;reals=[-radius,0.0,0.0])

  geoL   = egads.makeGeometry(context, egads.CURVE, egads.LINE, [radius,0.0,0.0,-2.0 * radius,0.0,0.0])

  t0,x0  = egads.invEvaluate(geoL,[ radius,0.0,0.0])
  t1,x1  = egads.invEvaluate(geoL,[-radius,0.0,0.0])

  edgeL  = egads.makeTopology!(context, egads.EDGE ;geom=geoL, mtype=egads.TWONODE, reals=[t0,t1], children=[node1,node2])


  stat   =  egads.inTopology(edgeL, x0)

  dataC  = [0.0, 0.0, 0.0, 1.0, 0.0,0.0, 0.0,1.0,0.0, radius]
  geoC   = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, dataC)

  t0,x0  = egads.invEvaluate(geoC,[ radius,0.0,0.0])
  t1,x1  = egads.invEvaluate(geoC,[-radius,0.0,0.0])
  edgeC  = egads.makeTopology!(context, egads.EDGE; geom=geoC, mtype=egads.TWONODE, reals=[t0,t1], children=[node1,node2])
  loop,edges   = egads.makeLoop!([edgeC,edgeL])

  face  = egads.makeFace(loop, egads.SFORWARD)
  uv,x0 = egads.invEvaluate(face,[ radius,0.0,0.0])
  stat  = egads.inFace(face,uv)
  return true
end

#------------------------------------------------------------------------------
function test_inTopology(context)
  box   = egads.makeSolidBody(context, egads.BOX, [0,0,0, 1,1,1])
  faces = egads.getBodyTopos(box,egads.FACE)
  info  = egads.getRange(faces[1])

  bbox  = egads.getBoundingBox(faces[1])
  @test egads.inTopology(faces[1],0.5 .* [sum(bbox[1]),sum(bbox[2]),sum(bbox[3])])
  @test egads.inFace(faces[1],0.5 .* [sum(info.range[1]),sum(info.range[2])])
  return true
end

#------------------------------------------------------------------------------
function test_booleanOperations(context)

  box1    = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2    = egads.makeSolidBody(context,egads.BOX, [-0.5,-0.5,-0.5,1.0,1.0,1.0])
  bt      = [egads.SUBTRACTION, egads.INTERSECTION, egads.FUSION]
  #boolbox = Array{egads.Ego}(undef, 3)
  #boolbox[1] = egads.generalBoolean(box1, box2, bt[1])
  #boolbox[2] = egads.generalBoolean(box1, box2, bt[2])
  #boolbox[3] = egads.generalBoolean(box1, box2, bt[3])
  #boolbox = egads.generalBoolean(box1, box2, bt[2])
  #boolbox = egads.generalBoolean(box1, box2, bt[3])
  boolbox = [egads.generalBoolean(box1, box2, oper) for oper in bt]
  return true
end
#------------------------------------------------------------------------------
function test_intersection(context)

  wBody  = make_wire(context, [[0,0,0] [1,0,0] [1,1,0] [0,1,0]], egads.BODY)
  box    = egads.makeSolidBody(context,egads.BOX, [0.25,0.25,-0.25, 0.5,0.5,0.5])
  inter  = egads.intersection(box, wBody)
  box2   = egads.imprintBody(box, inter.pairs)
  faces1 = egads.getBodyTopos(box2, egads.FACE)
  @test length(faces1) == 10
  return true
end

#------------------------------------------------------------------------------
function test_filletBody(context)

  box1   = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  edges1 = egads.getBodyTopos(box1, egads.EDGE)
  box4   = egads.filletBody(box1, edges1, 0.1)
  info   = egads.getInfo(box4.result)
  @test info.oclass == egads.BODY && info.mtype == egads.SOLIDBODY
  return true
end

#------------------------------------------------------------------------------
function test_fuseSheets(context)
  box1    = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2    = egads.makeSolidBody(context,egads.BOX, [-0.5,-0.5,-0.5,1.0,1.0,1.0])

  shell1 = egads.getBodyTopos(box1, egads.SHELL)
  shell2 = egads.getBodyTopos(box2, egads.SHELL)

  sheet1 = egads.makeTopology!(context,egads.BODY;mtype = egads.SHEETBODY, children = shell1[1])
  sheet2 = egads.makeTopology!(context,egads.BODY;mtype = egads.SHEETBODY, children = shell2[1])

  fuse = egads.fuseSheets(sheet1, sheet2)
  return true
end

#------------------------------------------------------------------------------
function test_chamferHollow(context)

  box1     = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  edges1   = egads.getBodyTopos(box1,egads.EDGE)
  b1_faces = [egads.getBodyTopos(box1, egads.FACE ; ref = k)[1] for k in edges1]
  cham     = egads.chamferBody(box1, edges1, b1_faces, 0.1, 0.15)
  info     = egads.getInfo(cham.result)
  @test info.oclass == egads.BODY && info.mtype == egads.SOLIDBODY


  holl = egads.hollowBody(box1, nothing, 0.1, 0)
  info = egads.getInfo(holl.result)
  @test info.oclass == egads.BODY && info.mtype == egads.SOLIDBODY

  face = egads.getBodyTopos(box1, egads.FACE)[1]
  holl2 = egads.hollowBody(box1, face, 0.1, 0)

  info = egads.getInfo(holl2.result)
  @test info.oclass == egads.BODY && info.mtype == egads.SOLIDBODY
  return true
end

#------------------------------------------------------------------------------
function test_sweepExtrude(context)

  faceBODY = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.SURFACE)

  xyz   = [ [0.0,0.0,0.0] [0.0,0.0,1.0]]
  nodes = [egads.makeTopology!(context, egads.NODE; reals = xyz[:,k]) for k = 1:2]

  rdata = [0.0, 0.0,0.0, 0.0,0.0,1.0]      # Line data (point and direction)
  line  = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)
  edge  = egads.makeTopology!(context, egads.EDGE ; mtype = egads.TWONODE, geom=line, children=nodes, reals=[0.0,1.0])


  sweep  = egads.sweep(faceBODY, edge, 0)
  sFaces = egads.getBodyTopos(sweep,egads.FACE)
  lf = length(sFaces)
  @test lf == 6

  extrude = egads.extrude(faceBODY, 1.0, [0.0,0.0,1.0])
  eFaces  = egads.getBodyTopos(extrude, egads.FACE)
  le      = length(eFaces)
  @test le == 6
  return true
end

#------------------------------------------------------------------------------

function test_effective(context)

  dataC  = [0.0, 0.0, 0.0, 1.0, 0.0,0.0, 0.0,1.0,0.0, 1.0]
  circle = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, dataC)

  points = [[1.0,0.0,0.0] [0.0,1.0,0.0] [-1.0,0.0,0.0] [0.0,-1.0,0.0] ]
  nodes  = [egads.makeTopology!(context, egads.NODE ; reals = points[:,k]) for k = 1:4]

  ts     = [0.0:0.5 * pi:2 * pi ;]
  cyc    = vcat([1:1:4;],1)
  edges  = [egads.makeTopology!(context, egads.EDGE; geom = circle, mtype = egads.TWONODE, children = nodes[cyc[k:k+1]], reals=ts[k:k+1]) for k = 1:4]

  loop   = egads.makeTopology!(context, egads.LOOP; mtype = egads.CLOSED, children = edges, senses=[1,1,1,1])

  wire   = egads.makeTopology!(context,egads.BODY;mtype = egads.WIREBODY, children = loop)

  mat1   = [ [1.0,0.0,0.0,0.0] [0.0,1.0,0.0,0.0] [0.0,0.0,1.0,1.0]]
  xform  = egads.makeTransform(context,mat1)

  secs = [wire]
  [push!(secs, egads.copyObject(secs[i-1] ,other=xform)) for i = 2:4]

  ruled = egads.ruled(secs)
  faces = egads.getBodyTopos(ruled, egads.FACE)
  lf    = length(faces)
  @test 12 == lf
  [egads.attributeAdd!(faces[k], "faceMap",1) for k = 1:2]

  tess  = egads.makeTessBody(ruled, [0.1, 0.1, 15])

  ebody = egads.initEBody(tess, 10.0)

  eface = egads.makeEFace(ebody, faces[10:11])
  info = egads.getInfo(eface)
  @test info.oclass == egads.EFACE && info.mtype == egads.SFORWARD

  efaces = egads.makeAttrEFaces(ebody,"faceMap")
  le     = length(efaces)
  @test 1 == le

  info = egads.getInfo(efaces[1])
  @test info.oclass == egads.EFACE && info.mtype == egads.SFORWARD
  egads.finishEBody!(ebody)
  bTopo = egads.getTopology(ebody)

  bl    = length(bTopo.children)
  bc    = bTopo.oclass
  bg    = bTopo.geom


  @test bc == egads.EBODY && bl == 1 && bg._ego  == ruled._ego

  for shell in bTopo.children
    sTopo = egads.getTopology(shell)
    ls =  length(sTopo.children)
    @test sTopo.oclass == egads.ESHELL &&  ls == 10
    for face in sTopo.children
      fTopo = egads.getTopology(face)
      lf = length(fTopo.children)
      @test fTopo.oclass == egads.EFACE && lf == 1
      for loop in fTopo.children
        lTopo = egads.getTopology(loop)
        @test lTopo.oclass == egads.ELOOPX
      end
    end
  end

  # check getting the original FACES
  faces  = egads.getBodyTopos(ebody, egads.FACE)
  lf     = length(faces)
  @test 12 == lf

  efaces = egads.getBodyTopos(ebody, egads.EFACE)
  le     = length(efaces)
  @test 10 == length(efaces)

  eUV    = egads.getRange(efaces[1])
  fMap   = egads.effectiveMap(efaces[1], 0.5 .* [sum(eUV.range[1,:]), sum(eUV.range[2,:])])

  eedges = egads.getBodyTopos(ebody, egads.EEDGE)
  eList  = egads.effectiveEdgeList(eedges[1])

  eT   = egads.getRange(eedges[1])
  eMap = egads.effectiveMap(eedges[1], 0.5 .* sum(eT.range))
  return true
end

#------------------------------------------------------------------------------
function test_blend_and_ruled(context)

  # first section is a NODE
  node1 = egads.makeTopology!(context, egads.NODE ; reals=[0.5,0.5,-1])

  # FaceBody for 2nd section
  face1 = make_wire(context, [[0,0,0] [1,0,0] [1,1,0] [0,1,0]], egads.SURFACE)
  mat1     = [ [1.0,0.0,0.0,0.0] [0.0,1.0,0.0,0.0] [0.0,0.0,1.0,1.0]]
  xform    = egads.makeTransform(context, mat1)

  # Translate 2nd section to make 3rd
  face2 = egads.copyObject(face1,other=xform)

  # Finish with a NODE
  node2 = egads.makeTopology!(context, egads.NODE ; reals=[0.5,0.5,2])

  # sections for blend
  secs = [node1,face1,face2,node2]
  blend1 = egads.blend(secs)
  data   = egads.getInfo(blend1)
  @test data.oclass == egads.BODY && data.mtype == egads.SOLIDBODY

  rc = [0.2, 1, 0, 0,  0.4, 0, 1, 0]
  blend2 = egads.blend(secs; rc1=rc)
  data = egads.getInfo(blend2)
  @test data.oclass == egads.BODY && data.mtype == egads.SOLIDBODY

  blend3 = egads.blend(secs; rc1=rc, rc2=rc)
  data = egads.getInfo(blend3)
  @test data.oclass == egads.BODY && data.mtype == egads.SOLIDBODY

  blend4 = egads.blend(secs; rc2=rc)
  data = egads.getInfo(blend4)
  @test data.oclass == egads.BODY && data.mtype == egads.SOLIDBODY

  ruled = egads.ruled(secs)
  data = egads.getInfo(ruled)
  @test data.oclass == egads.BODY && data.mtype == egads.SOLIDBODY
  return true
end

#------------------------------------------------------------------------------
function test_rotate(context)
  faceBODY = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.SURFACE)

  body     = egads.rotate(faceBODY, 180, [[0.0,-1.0,0.0]  [1.0,0.0,0.0]])

  faces = egads.getBodyTopos(body,egads.FACE)
  @test 6 == length(faces)
  return true
end

#------------------------------------------------------------------------------
function test_mapBody(context)

  box1 = egads.makeSolidBody(context, egads.BOX,[0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context, egads.BOX,[0.1,0.1,0.1,1.0,1.0,1.0])

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  for i=1:length(faces1)
    egads.attributeAdd!(faces1[i], "faceMap", [i])
    egads.attributeAdd!(faces2[i], "faceMap", [i])
  end

  mapB  = egads.mapBody(box1,box2, "faceMap")
  tess1 = egads.makeTessBody(box1,[0.1, 0.1, 15.])

  tess2   = egads.mapTessBody(tess1, box2)
  tessF   = egads.getTessFace(tess1, 1)
  ifaces  = fill(1, Int(floor(length(tessF.uv[1,:])/2)))
  locTess = egads.locateTessBody(tess2, ifaces, tessF.uv)
  return true
end

#------------------------------------------------------------------------------
function test_fitTriangles(context)

  box  = egads.makeSolidBody(context, egads.BOX,[0.0,0.0,0.0,1.0,1.0,1.0])
  tess = egads.makeTessBody(box,[0.05, 0.1, 30.])
  tessData = egads.getTessFace(tess, 1)

  triBox = egads.fitTriangles(context, tessData.xyz, tessData.tris; tol = 1.e-8)
  info   = egads.getInfo(triBox)

  @test info.mtype == egads.BSPLINE && info.oclass == egads.SURFACE

  triBox2 = egads.fitTriangles(context, tessData.xyz, tessData.tris; tessData.tric, tol = 1.e-8)
  info    = egads.getInfo(triBox2)
  @test info.mtype == egads.BSPLINE && info.oclass == egads.SURFACE
  return true
end

#------------------------------------------------------------------------------
function test_addKnots(context)

  #---------------
  # Curve example
  #---------------
  dataE = [1.0,0.0,0.0,1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4.0, 0.5]
  geoE  = egads.makeGeometry(context, egads.CURVE, egads.ELLIPSE, dataE)

  rInfo = egads.getRange(geoE)
  bspline = egads.convertToBSplineRange(geoE, rInfo.range)
  oclass, mtype, geom, ints, reals = egads.getGeometry(bspline)

  degU    = ints[2]
  nKnotsU = ints[4]
  Us = zeros(nKnotsU-2*degU+1)
  for i = degU:nKnotsU-degU
      Us[i-degU+1] = (reals[i+1]+reals[i])/2
  end

  # Only add knots
  minDegU = degU
  newCurve = egads.addKnots(bspline, minDegU, Us, 0, nothing)
  
  new_oclass, new_mtype, geom, new_ints, new_reals = egads.getGeometry(newCurve)

  @test oclass == new_oclass
  @test mtype == new_mtype
  @test minDegU == new_ints[2]
  @test ints[4]+length(Us) == new_ints[4]

  # Increase degree
  minDegU = degU+1
  newCurve = egads.addKnots(bspline, minDegU, nothing, 0, nothing)
  
  new_oclass, new_mtype, geom, new_ints, new_reals = egads.getGeometry(newCurve)

  @test oclass == new_oclass
  @test mtype == new_mtype
  @test minDegU == new_ints[2]

  #---------------
  # Surface example
  #---------------
  nCPu = 4;
  nCPv = 4;
  pts = [[0.00  0.00  0.00] 
         [1.00  0.00  0.10] 
         [1.50  1.00  0.70] 
         [0.25  0.75  0.60]

         [0.00  0.00  1.00] 
         [1.00  0.00  1.10] 
         [1.50  1.00  1.70] 
         [0.25  0.75  1.60] 

         [0.00  0.00  2.00] 
         [1.00  0.00  2.10] 
         [1.50  1.00  2.70] 
         [0.25  0.75  2.60] 

         [0.00  0.00  3.00] 
         [1.00  0.00  3.10] 
         [1.50  1.00  3.70] 
         [0.25  0.75  3.60]]

  bspline = egads.approximate(context, [nCPu, nCPv], pts)
  oclass, mtype, geom, ints, reals = egads.getGeometry(bspline)

  degU    = ints[2]
  nKnotsU = ints[4]
  Us = zeros(nKnotsU-2*degU-1)
  for i = degU+1:nKnotsU-degU-1
    Us[i-degU] = (reals[i+1]+reals[i])/2
  end

  degV    = ints[5]
  nKnotsV = ints[7]
  Vs = zeros(nKnotsV-2*degV-1)
  for i = nKnotsU+degV+1:nKnotsU+nKnotsV-degV-1
      Vs[i-degV-nKnotsU] = (reals[i+1]+reals[i])/2
  end

  # # Only add knots
  minDegU = degU
  minDegV = degV
  newSurface = egads.addKnots(bspline, minDegU, Us, minDegV, Vs)
  
  new_oclass, new_mtype, geom, new_ints, new_reals = egads.getGeometry(newSurface)

  @test oclass == new_oclass
  @test mtype == new_mtype
  @test minDegU == new_ints[2]
  @test minDegV == new_ints[5]
  @test ints[4]+length(Us) == new_ints[4]
  @test ints[7]+length(Us) == new_ints[7]

  # Increase degree
  minDegU = degU+1
  minDegV = degV+2
  newSurface = egads.addKnots(bspline, minDegU, nothing, minDegV, nothing)
  
  new_oclass, new_mtype, geom, new_ints, new_reals = egads.getGeometry(newSurface)

  @test oclass == new_oclass
  @test mtype == new_mtype
  @test minDegU == new_ints[2]
  @test minDegV == new_ints[5]
end

#------------------------------------------------------------------------------
function test_quadTess(context)
  box   = egads.makeSolidBody(context, egads.BOX,[0.0,0.0,0.0,1.0,1.0,1.0])
  tess1 = egads.makeTessBody(box,[0.1, 0.1, 30.])

  edges  = egads.getBodyTopos(box, egads.EDGE)
  facedg = egads.getBodyTopos(box, egads.FACE; ref = edges[1])
  egads.remakeTess!(tess1, vcat(facedg, edges[1]) ,[0.2,0.1,10.0])

  tessData = egads.getTessFace(tess1, 5)
  quadBox  = egads.quadTess(tess1)

  tess2    = egads.makeTessBody(box,[0.1, 0.1, 30.])
  egads.makeQuads!(tess2,[0.1, 0.1, 30.], 1)

  qTess = egads.getQuads(tess2, 1)

  patch = egads.getPatch(tess2,1,qTess.npatch)

  fList = egads.getTessQuads(tess2)
  @test length(fList) == 1 && fList[1] == 1

  plane = egads.makeGeometry(context,egads.SURFACE,egads.PLANE,[0.0,0.0,0.0, 1.0,0.0,0.0 ,0.0,1.0,0.0])
  pTess = egads.makeTessGeom( plane, [0.0,1.0, 0.0,1.0], [5,5])

  planeT = egads.getTessGeom(pTess)
  @test planeT.sizes == [5,5]
  return true
end

#------------------------------------------------------------------------------
function test_tessBody(context)
  sphere  = egads.makeSolidBody(context,egads.SPHERE,[0.0,0.0,0.0,1.0])

  tParams = [0.2, 0.2, 30.]
  tess    = egads.makeTessBody(sphere, tParams)

  tEdge  = egads.getTessEdge(tess, 2)
  infoL  = egads.getTessLoops(tess, 1)

  infoTess = egads.statusTessBody(tess)

  globalT  = [egads.getGlobal(tess,j) for j = 1:infoTess.npts]
  faces    = egads.getBodyTopos(sphere, egads.FACE)

  # Build up a tessellation (by copying)
  tess2 = egads.initTessBody(sphere)

  edges = egads.getBodyTopos(sphere, egads.EDGE)

  for i=1:length(edges)
      info = egads.getInfo(edges[i])
      info.mtype == egads.DEGENERATE && continue
      vals = egads.getTessEdge(tess,i)
      egads.setTessEdge!(tess2, i, vals.xyzs, vals.ts)
      vals2 = egads.getTessEdge(tess2, i)
      tdiff   = sum(abs.(   vals.ts .- vals2.ts))
      xyzdiff = sum((vals.xyzs .- vals2.xyzs).^2)
      @test tdiff + xyzdiff < 1.e-14
  end
  for j = 1:length(faces)
    faceT = egads.getTessFace(tess, j)
    glob  = [egads.localToGlobal(tess,j, faceT.tris[i]) for i=1:length(faceT.tris)]
    egads.setTessFace!(tess2,j, faceT.xyz, faceT.uv, faceT.tris)
  end

  infoTess2 = egads.statusTessBody(tess2)

  # Open up a tessellation and finish it
  egads.openTessBody(tess2)
  egads.finishTess!(tess2,[0.1,0.1,20.0])
  return true
end

#------------------------------------------------------------------------------
function test_modifyTess(context)
  radius = 1.25
  body   = egads.makeSolidBody(context, egads.SPHERE, [0.0,0.0,0.0,radius])
  tess   = egads.makeTessBody(body,[0.1, 0.1, 20.])
  eIndex = 2
  tEdge  = egads.getTessEdge(tess,eIndex)
  vIndex = Int(floor(length(tEdge.ts)/2))
  egads.moveEdgeVert!(tess, eIndex, vIndex,0.5*(tEdge.ts[vIndex-1] + tEdge.ts[vIndex]))

  egads.deleteEdgeVert!(tess, eIndex, vIndex, 1)

  tEdge = egads.getTessEdge(tess,eIndex)

  egads.insertEdgeVerts!(tess,eIndex, vIndex, 0.5 * (tEdge.ts[vIndex] + tEdge.ts[vIndex + 1]))

  massP = egads.tessMassProps(tess)
  eA = abs(massP.area - 4 * pi * radius^2)
  eV = abs(massP.volume - 4 / 3 * pi * radius^3)
  @test eA < 0.1 && eV < 0.1 && sum((massP.CG).^2) < 0.1
  return true
end

#------------------------------------------------------------------------------
function test_memo(context)
  ctxt  = egads.Context()
  level = egads.setOutLevel(ctxt,1)
  face = make_wire(ctxt, [[0.0,0,0] [1.0,0,0] [1.0,1.0,0] [0,1.0,0]], egads.SURFACE)
  face = egads.makeSolidBody(context, egads.SPHERE, [0.0,0.0,0.0,1.0])
  face = nothing
  GC.gc()
  info = egads.getInfo(ctxt)
  @test info.next == nothing
  return true
end

#------------------------------------------------------------------------------
function test_NMwireBody(context)
	# construction data
	context  = egads.Context()
        xy  = [[0.0,0.0,0.0] [1.0,0.0,0.0] [1.0,1.0,0.0] [0.0,1.0,0.0]]
	c14 = vcat([1:4;],1)

        # make Nodes
        nodes = [egads.makeTopology!(context,egads.NODE; reals=xy[:,j]) for j = 1 : 4]
  	dir   = hcat([ xy[:,c14[j+1]] .- xy[:,c14[j]] for j = 1:4]...)
        lines = [egads.makeGeometry(context,egads.CURVE, egads.LINE, vcat(xy[:,j],dir[:,j])) for j = 1:4]
	range = zeros(2,4)
	for j in 1:4
    		aux = egads.invEvaluate(lines[j],xy[:,c14[j]])
		range[1,j] = aux.params
		aux = egads.invEvaluate(lines[j],xy[:,c14[j+1]])
		range[2,j] = aux.params
	end
	edges = [egads.makeTopology!(context,egads.EDGE; mtype = egads.TWONODE, geom=lines[j],
                 children = [nodes[c14[j]], nodes[c14[j+1]]], reals = range[:,j]) for j = 1:4]

	wbody = egads.makeNmWireBody(edges)
        data  = egads.getInfo(wbody)
        @test data.mtype == egads.WIREBODY && data.oclass == egads.BODY
	return true
end


@time @testset " Julia egads API " begin
  # RUN ALL FUNCTION TESTS
  for name in names(@__MODULE__; all = true)
    if !startswith("$(name)", "test_")
        continue
    end
    @info " RUNNING TEST $(name) "
    context = egads.Context()
    level = egads.setOutLevel(context,0)
    getfield(@__MODULE__, name)(context)
    #GC.gc()
    #info = egads.getInfo(context)
    #@test info.next === nothing
  end
end
