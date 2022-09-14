###########################################################################
#                                                                         #
#                         unit testing for jlEGADS                        #
#                                                                         #
#                                                                         #
#      Copyright 2011-2021, Massachusetts Institute of Technology         #
#      Licensed under The GNU Lesser General Public License, version 2.1  #
#      See http://www.opensource.org/licenses/lgpl-2.1.php                #
#                                                                         #
#                                                                         #
###########################################################################

using Test

push!(LOAD_PATH, ".")
using egads

loc_folder = ENV["jlEGADS"]*"/test"

function make_wire(context, xyz, class = egads.LOOP)
  d,n = size(xyz)
  d != 3 && @error " make_wire data needst to be a [3][n] matrix "
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

function test_revision()
  maxR, minR, occ = egads.revision()
  @info "********* EGADS version $maxR.$minR OCC $occ *********"
  return true
end

function test_makeLoopFace()
  xyz   = Cdouble.([[0,0,0] [0,1,0] [1,1,0] [1,0,0]])
  d, n  = size(xyz)
  c1n   = vcat([1:1:n;],1)
  nodes = [egads.makeTopology!(context, egads.NODE ; reals=xyz[:,j], name = "NODE_"*string(j)) for j = 1:n]

  dir   = hcat([ xyz[:,c1n[j+1]] .- xyz[:,c1n[j]] for j = 1:n]...)
  lines = [egads.makeGeometry(context,egads.CURVE, egads.LINE, vcat(xyz[:,j],dir[:,j]) ; name = "line_"*string(j) ) for j = 1:n]
  range = zeros(2,n)

  for j in 1:n
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j]])
    range[1,j] = aux.params
    aux = egads.invEvaluate(lines[j],xyz[:,c1n[j+1]])
    range[2,j] = aux.params
  end

  edges = [egads.makeTopology!(context,egads.EDGE; mtype = egads.TWONODE, geom=lines[j],
           children = [nodes[c1n[j]], nodes[c1n[j+1]]], reals = range[:,j], name = "EDGE_"*string(j)) for j = 1:n]

  loop1 = egads.makeTopology!(context,egads.LOOP; mtype=egads.CLOSED,
          children = edges, senses=fill(egads.SFORWARD,n), name = "LOOP1")

  loop2 = egads.makeLoop!(edges; name = "LOOP2")


  a = egads.getArea(loop1) ; b = egads.getArea(loop2)
  if (a - b > 1.e-10)
     @warn "Loops don't give same area !! $a != $b "
  end

  faceL1 = egads.makeFace(loop1, egads.SFORWARD)
  faceL2 = egads.makeFace(loop2, egads.SFORWARD)

  data   = egads.getInfo(faceL1)
  (data.oclass != egads.FACE || data.mtype != egads.SFORWARD) &&
    @warn "makeGeometry oclass (FACE, SFORWARD) != $(data.oclass),$(data.mtype) !"

  a1 = egads.getArea(faceL1) ; a2 = egads.getArea(faceL2)
  (a1 - a2 > 1.e-10) && @warn "Faces don't give the same area: $a1 != $a2"

  plane = egads.makeGeometry(context, egads.SURFACE, egads.PLANE,[0,0,0,1,0,0,0,1,0]; name = "FACE PLANE")
  facePlane = egads.makeFace(plane, egads.SFORWARD ; rdata = [0,1, 0,1])

  data = egads.getInfo(facePlane)
  (data.oclass != egads.FACE || data.mtype != egads.SFORWARD) &&
    @warn "makeGeometry oclass (FACE, SFORWARD) != $(data.oclass),$(data.mtype) !"
  return true
end

function test_makeTransform()
  mat1  = [ [1.0,0.0,0.0,4.0] [0.0,1.0,0.0,8.0] [0.0,0.0,1.0,12.0]]
  trans = egads.makeTransform(context,mat1)
  mat2  = egads.getTransform(trans)
  mat1 != mat2 && @warn "getTransform matrices are different! "
  egads.updateThread!(context)
return true
end

function test_makeGeometry_curves()

  xyz    = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  cent   =  [0.0,0.0,0.0]
  xyzAX  = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]

  radius = 0.5
  rM     = 1.5
  rm     = 1.0

  #CIRCLE
  circle = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, vcat(xyz[1], xyzAX[1], radius))
  geoC   = egads.getGeometry(circle)
  geoC.oclass != egads.CURVE  && @warn "getGeometry oclass $(geoC.oclass)!= CURVE "
  geoC.mtype  != egads.CIRCLE && @warn "getGeometry mtype $(geoC.mtype) != CIRCLE "
  #ELLIPSE
  ellipse = egads.makeGeometry(context, egads.CURVE, egads.ELLIPSE,  vcat(cent, xyzAX[1], xyzAX[2], rM,rm))
  geoE    = egads.getGeometry(ellipse)
  geoE.oclass != egads.CURVE   && @warn "getGeometry oclass $(geoE.oclass)!= CURVE "
  geoE.mtype  != egads.ELLIPSE && @warn "getGeometry mtype $(geoE.mtype) != CIRCLE "
  geoE.reals[1:3] != xyz[1]    && @warn "getGeometry reals $(geoE.reals[1:3])!= [0.0,0.0,0.0]"

  #PARABOLA
  parabola = egads.makeGeometry(context, egads.CURVE, egads.PARABOLA,  vcat(cent, xyzAX[1], xyzAX[2],rM))
  geoP     = egads.getGeometry(parabola)
  geoP.oclass != egads.CURVE    && @warn "getGeometry oclass $(geoP.oclass)!=CURVE "
  geoP.mtype  != egads.PARABOLA && @warn "getGeometry mtype $(geoP.mtype) !=PARABOLA "
  geoP.reals[end] != rM         && @warn "getGeometry reals[end] $(geoP.reals[end])!=$rM"

  #HIPERBOLA
  hiperbola = egads.makeGeometry(context, egads.CURVE, egads.HYPERBOLA,  vcat(cent, xyzAX[1], xyzAX[2],rm, rM))
  geoH      = egads.getGeometry(hiperbola)
  geoH.oclass != egads.CURVE     && @warn "getGeometry oclass $(geoH.oclass)!= CURVE "
  geoH.mtype  != egads.HYPERBOLA && @warn "getGeometry mtype $(geoH.mtype) != HIPERBOLA "
  geoH.reals[1:3] != xyz[1]      && @warn "getGeometry reals[1:3] $(geoH.reals[1:3]) != [0.0,0.0,0.0]"

  #OFFSET CURVE

  offset  = egads.makeGeometry(context, egads.CURVE, egads.OFFSET,  vcat(xyzAX[1], 2.0), rGeom = circle)

  geoO    = egads.getGeometry(offset)
  geoO.oclass != egads.CURVE        && @warn "getGeometry oclass $(geoO.oclass) != CURVE "
  geoO.mtype  != egads.OFFSET       && @warn "getGeometry mtype $(geoO.mtype) != OFFSET "
  geoO.reals  != vcat(xyzAX[1],2.0) && @warn "getGeometry reals[1:3] $(geoO.reals[1:3]) != [1.0,0.0,0.0,2.0]"
  return true
end

function test_makeGeometry_surfaces()

  xyz    = [[0.0,0.0,0.0], [1.0,0.0,0.0]]
  cent   =  [0.0,0.0,0.0]
  xyzAX  = [[1.0,0.0,0.0], [0.0,1.0,0.0], [0.0,0.0,1.0]]
  radius = 0.5
  rM     = 1.5
  rm     = 1.0
  angle  = 30.0 * 3.14/180.
  #SPHERE
  sphere = egads.makeGeometry(context, egads.SURFACE, egads.SPHERICAL, vcat(xyz[1], xyzAX[1], radius))
  geoS   = egads.getGeometry(sphere)
  geoS.oclass != egads.SURFACE   && @warn "getGeometry oclass $(geoS.oclass) != SURFACE "
  geoS.mtype  != egads.SPHERICAL && @warn "getGeometry mtype $(geoS.mtype) != SPHERICAL "

  #CONICAL
  conical = egads.makeGeometry(context, egads.SURFACE, egads.CONICAL,  vcat(cent, xyzAX[1], xyzAX[2],xyzAX[3], angle, rM))
  geoC    = egads.getGeometry(conical)
  geoC.oclass != egads.SURFACE && @warn "getGeometry oclass $(geoC.oclass) != SURFACE "
  geoC.mtype  != egads.CONICAL && @warn "getGeometry mtype $(geoC.mtype) != CONICAL "
  geoC.reals[4:6] != xyzAX[1]  && @warn "getGeometry reals[4:6] $(geoC.reals[4:6]) != [1.0,0.0,0.0]"

  #CYLINDER
  cylinder = egads.makeGeometry(context, egads.SURFACE, egads.CYLINDRICAL,  vcat(cent, xyzAX[1], xyzAX[2],xyzAX[3], rM))
  geoC     = egads.getGeometry(cylinder)
  geoC.oclass != egads.SURFACE     && @warn "getGeometry oclass $(geoC.oclass) != SURFACE "
  geoC.mtype  != egads.CYLINDRICAL && @warn "getGeometry mtype $(geoC.mtype) != CYLINDRICAL"

  #TOROIDAL
  toroidal = egads.makeGeometry(context, egads.SURFACE, egads.TOROIDAL,  vcat(cent, xyzAX[1], xyzAX[2],xyzAX[3], rm, rM))
  geoT     = egads.getGeometry(toroidal)
  geoT.oclass != egads.SURFACE       && @warn "getGeometry oclass $(geoT.oclass) != SURFACE "
  geoT.mtype  != egads.TOROIDAL      && @warn "getGeometry mtype $(geoT.mtype) != TOROIDAL "
  geoT.reals[4:12] != vcat(xyzAX...) && @warn "getGeometry reals[4:12] $(geoT.reals[4:12]) != $xyzAX"

  #SURFACE REVOLUTION
  circle = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, vcat(xyz[1], xyzAX[1], radius))
  revolution = egads.makeGeometry(context, egads.SURFACE, egads.REVOLUTION,  vcat(cent, xyzAX[3]); rGeom = circle)
  geoR       = egads.getGeometry(revolution)
  geoR.oclass != egads.SURFACE    && @warn "getGeometry oclass $(geoR.oclass) != SURFACE "
  geoR.mtype  != egads.REVOLUTION && @warn "getGeometry mtype $(geoR.mtype) != REVOLUTION "
  geoR.reals[1:3] != cent         && @warn "getGeometry reals[1:3] $(geoR.reals[1:3]) != [0.0,0.0,0.0]"

  #EXTRUSION SURFACE
  extrusion  = egads.makeGeometry(context, egads.SURFACE, egads.EXTRUSION, xyzAX[3]; rGeom = circle)
  geoE       = egads.getGeometry(extrusion)
  geoE.oclass != egads.SURFACE   && @warn "getGeometry oclass $(geoE.oclass) != SURFACE "
  geoE.mtype  != egads.EXTRUSION && @warn "getGeometry mtype $(geoE.mtype) != REVOLUTION "
  geoE.reals  != xyzAX[3]        && @warn "getGeometry reals $(geoE.reals) != [0.0,0.0,1.0]"


  #OFFSET CURVE
  offset  = egads.makeGeometry(context, egads.SURFACE, egads.OFFSET, 2.0; rGeom = cylinder)
  geoO    = egads.getGeometry(offset)
  geoO.oclass != egads.SURFACE  && @warn "getGeometry oclass $(geoO.oclass) != SURFACE "
  geoO.mtype  != egads.OFFSET   && @warn "getGeometry mtype $(geoO.mtype) != OFFSET "
  geoO.reals  != 2.0            && @warn "getGeometry reals $(geoO.reals) != 2.0"
  return true
end

function test_makeBezier()
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
  geoB.oclass != egads.SURFACE && @warn "getGeometry oclass $(geoB.oclass)!= SURFACE "
  geoB.mtype  != egads.BEZIER  && @warn "getGeometry mtype  $(geoB.mtype) != BEZIER "
  geoB.ints   != idata         && @warn "getGeometry ints  $(geoB.ints) != $idata"

  # make Nodes
  rn    = hcat(rdata[:,1], rdata[:,4], rdata[:,16] , rdata[:,12] )
  nodes = [egads.makeTopology!(context, egads.NODE ; reals=rn[:,j]) for j = 1:4]
  #BEZIER CURVE
  idata = [0, nCPu-1,nCPu]
  idx   = [[1,2,3,4], [4,8,12,16], [16,15,14,13],[13,9,5,1]]
  beC   = [egads.makeGeometry(context, egads.CURVE, egads.BEZIER, rdata[:,idx[j]] ; ints = idata) for j =1:4]
  geoB   = egads.getGeometry(beC[1])
  geoB.oclass != egads.CURVE  && @warn "getGeometry oclass $(geoB.oclass)!=CURVE "
  geoB.mtype  != egads.BEZIER && @warn "getGeometry mtype $(geoB.mtype) !=BEZIER "
  geoB.ints   != idata        && @warn "getGeometry ints $(geoB.ints) != $idata"
  c1n = vcat([1:4;],1)
  edges   = [egads.makeTopology!(context, egads.EDGE; mtype = egads.TWONODE, geom = beC[j], children=vcat(nodes[c1n[j]],nodes[c1n[j+1]]), reals=[0.0,1.0]) for j = 1:4]
  pcurves = [egads.otherCurve(beS, e) for e in edges]
  edges   = vcat(edges, pcurves)

  loopF   = egads.makeTopology!(context, egads.LOOP ; mtype = egads.CLOSED, geom=beS, children= edges, senses = fill(egads.SFORWARD,4) )
  (egads.getArea(loopF) < 0.0) &&  @warn "getArea of forward loop is < 0 !!  "

  faceF = egads.makeTopology!(context, egads.FACE; mtype = egads.SFORWARD, geom=beS, children= loopF)
  # making a face with the loop reversed
  piv = [3,2,1,4,7,6,5,8]


  loopR = egads.makeTopology!(context, egads.LOOP ; mtype = egads.CLOSED, geom=beS, children=edges[piv], senses=fill(egads.SREVERSE,4))

  (egads.getArea(loopR) > 0.0) &&  @warn "getArea of reverse loop is > 0 !!  "
  faceR = egads.makeTopology!(context, egads.FACE; mtype = egads.SREVERSE, geom=beS, children= loopR)
  (egads.getArea(faceF) - egads.getArea(faceR) > 1.e-10) &&  @warn "getArea of forward and reverse face are different ! "

  plane = egads.makeGeometry(context, egads.SURFACE , egads.PLANE, [0.0,0.0,0.0,1.0,0.0,0.0, 0.0,1.0,0.0])

  planeFace = egads.makeFace(plane, egads.SREVERSE ; rdata = [0 1 0 1])
  info = egads.getInfo(planeFace)

  (info.oclass != egads.FACE || info.mtype != egads.SREVERSE) &&
    @warn "makeFace wrong $(info.oclass),$(info.mtype) != (FACE, REVERSE) ! "

  a1 = egads.getArea(planeFace)
  a2 = egads.getArea(faceR)
  (a1 - a2 > 1.e-10) && @warn "getArea of reverse faces are different $a1 = $a2 ! "


  flip  = egads.flipObject(planeFace)

  a1   = egads.getArea(planeFace)

  (a1 - a2 > 1.e-10) && @warn "getArea of forward faces are different $a1 = $a2 ! "

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

function test_approximate()

  dataE = [1.0,0.0,0.0,1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 4.0, 0.5]
  geoE  = egads.makeGeometry(context, egads.CURVE, egads.ELLIPSE, dataE)

  nCP   = 20
  rInfo = egads.getRange(geoE)

  geoSpline = egads.convertToBSplineRange(geoE, rInfo.range)

  geoInfo = egads.getGeometry(geoSpline)
  if geoInfo.oclass != egads.CURVE || geoInfo.mtype != egads.BSPLINE
    @warn "convertToBsplineEdge (CURVE,BSPLINE) != $(geoInfo.oclass),$(geoInfo.mtype)"
  end
  h     = (rInfo.range[2] - rInfo.range[1]) / (nCP + 1)
  jeval = [egads.evaluate(geoE, rInfo.range[1] + h * j) for j in 1:nCP]
  data3D = hcat(getfield.(jeval,:X)...)

  # Fit spline through data points
  bspline = egads.approximate(context, [nCP,0], data3D)

  # LOOK BACKWARDS STUFF
  geoInfo = egads.getGeometry(bspline)

  if geoInfo.oclass != egads.CURVE || geoInfo.mtype != egads.BSPLINE
    @warn "convertToBsplineEdge (CURVE,BSPLINE) != $(geoInfo.oclass),$(geoInfo.mtype)"
  end

  for j in 1:nCP
    invJ  = egads.invEvaluate(bspline, data3D[:,j])
    e_xyz = sum((invJ.result .- data3D[:,j]).^2)
    e_xyz > 1.e-5 && @warn "invEvaluate: fit point $j error interpol $e_xyz"
  end
  dupSpline = egads.convertToBSpline(bspline)
  stat      = egads.isSame(bspline, dupSpline)
  return true
end

function test_skinning()

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

  egadsfile = loc_folder*"/skinning.egads"

  range, per = egads.getRange(bspline)
  face       = egads.makeFace(bspline, egads.SFORWARD ; rdata =  range)
  splineBody = egads.makeTopology!(context, egads.BODY;mtype=egads.FACEBODY,children=face)

  faces      = egads.getBodyTopos(splineBody,egads.FACE)

  body2 = egads.getBody(faces[1])

  stat = egads.isEquivalent(body2, splineBody)
  egads.saveModel!(body2,egadsfile,true)
  return true
end

function test_attribution()

  egads.setFullAttrs!(context, 1)

  box1 = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context,egads.BOX, [1.0,0.0,0.0,1.0,1.0,1.0])

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  egads.attributeAdd!(faces1[2], "faceMap", 2)
  egads.attributeAdd!(faces2[1], "faceMap", 1)

  mapB    = egads.mapBody(box1,box2, "faceMap")

  attrVal = egads.attributeRet(faces1[2], "faceMap")
  nAttr   = egads.attributeNum(faces1[2])

  if attrVal[1] != 2 || nAttr != 1
    @warn "Box 1 face 2: attributeAdd faceMap should be 2 != $attrVal TOTAL ATTRIBUTES $nAttr should be = 1"
  end
  attrVal = egads.attributeRet(faces2[1], "faceMap")
  nAttr   = egads.attributeNum(faces2[1])

  if attrVal[1] != 1 || nAttr != 1
    @warn "Box 1 face 2: attributeAdd faceMap should be 2 != $attrVal TOTAL ATTRIBUTES $nAttr should be = 1"
  end

  egads.attributeAddSeq!(box1,"faceDouble", [0.1,0.2])
  egads.attributeAddSeq!(box1,"faceDouble", [0.3,0.4])


  attrVal = egads.attributeRetSeq(box1, "faceDouble", 1)

  nAttr   = egads.attributeNumSeq(box1,"faceDouble")

  if attrVal != [0.1,0.2] || nAttr != 2
    @warn "Box 1 face 2: attributeAddSeq faceDouble should be [0.1,0.2] != $attrVal  TOTAL ATTRIBUTES $nAttr should be = 2"
  end

  egads.attributeAdd!(box1, "faceString", "myAttr")

  nAttr   = egads.attributeNum(box1)
  attrVal = egads.attributeRet(box1, "faceString")

  if attrVal != "myAttr" || nAttr != 3
    @warn "Box 1 face 2: attributeAdd faceString should be myAttr != $attrVal TOTAL ATTRIBUTES $nAttr  should be = 3"
  end

  info = egads.attributeGet(box1, 3)
  if info.name != "faceString" || info.attrVal != attrVal
    @warn "egads.attributeGet & Ret are different : NAMES $info.name != faceString ; VALS $(info.attrVal) != $attrVal"
  end

  egads.attributeDup!(box1,box2)
  nAttr2 = egads.attributeNum(box2)
  nAttr2 != nAttr && @warn "attributeDup didn't give the same number $nAttr2 != $nAttr"


  egads.attributeDel!(box1; name = "faceString")

  nAttr = egads.attributeNum(box1)

  nAttr != 2 && @warn "attributeDel left ATTRIBUTES $nAttr should be = 2"


  egads.attributeDel!(box1)
  nAttr = egads.attributeNum(box1)

  nAttr != 0 && @warn "attributeDel left ATTRIBUTES $nAttr should be = 0"
  # csystem

  attrCsys1 =  egads.CSys([0,0, 0, 1,0, 0,0,1,0])
  egads.attributeAdd!(box1,"faceCSYS", attrCsys1)

  data = egads.attributeRet(box1, "faceCSYS")
  (data.val!= attrCsys1.val || typeof(data) != egads.CSys)      && @warn "attributeRet for faceCYS should be $attrCsys1 != $data"

  data = egads.attributeRetSeq(box1, "faceCSYS", 1)
  (data.ortho != attrCsys1.ortho || typeof(data) != egads.CSys) && @warn " egads.attributeRetSeq for faceCYS should be $attrCsys1 != $data"

  name, data = egads.attributeGet(box1, 1)
  (name != "faceCSYS" || typeof(data) != egads.CSys) && @warn "attributeGet should be name faceCSYS type CSYS. Got $data"
  return true
end

function test_tolerance()

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

function test_make_export_IO_model()
  plane    = egads.makeGeometry(context, egads.SURFACE , egads.PLANE, [0.0,0.0,0.0,1.0,0.0,0.0, 0.0,1.0,0.0])
  pFace    = egads.makeFace(plane, egads.SREVERSE ; rdata = [0,1,0,1])
  faceBody = egads.makeTopology!(context, egads.BODY; mtype = egads.FACEBODY , children= pFace)

  model    = egads.makeTopology!(context, egads.MODEL; children = faceBody)
  egads.saveModel!(model, loc_folder*"/faceBody.egads", true)
  mSB   = egads.loadModel(context, 0, loc_folder*"/faceBody.egads")
  mTopo = egads.getTopology(mSB)

  info = egads.getInfo(mTopo.children[1])
  (info.oclass !=egads.BODY || info.mtype != egads.FACEBODY) && @warn "loadModel $(info.oclass), $(info.mtype) != (BODY, FACEBODY) !"
  face = egads.objectBodyTopo(mTopo.children[1],egads.FACE,1)
  !egads.isSame(face, pFace) && @warn "loadModel different body faces"

  bbox   = egads.getBoundingBox(mSB)
  stream = egads.exportModel(mSB)
  !(typeof(stream) <: String) && @warn "exportModel return type is $stream !=a string"

  finalize(model)
  finalize(faceBody)
  finalize(pFace)
  finalize(plane)
  finalize(mSB)
  return true
end

function test_geoProperties()
  #####   CREATE A LINE
  dataL  = [1.0,0.0,0.0,1.0, 1.0, 0.0]
  geoL   = egads.makeGeometry(context, egads.CURVE,egads.LINE, dataL)
  radius = 0.5
  dataC  = [0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0,0.0,radius]
  dupL   = egads.getGeometry(geoL)
  info   = egads.getInfo(geoL)

  if (info.oclass != dupL.oclass || info.mtype != dupL.mtype ||
      dupL.oclass != egads.CURVE || dupL.mtype != egads.LINE)
     @warn " getInfo + getGeometry: CURVE != $(info.oclass) || $(dupL.oclass) LINE != $(info.oclass)  $(dupL.oclass)"
  end

  #####   CREATE A CIRCLE
  geoC = egads.makeGeometry(context, egads.CURVE, egads.CIRCLE, dataC)
  refG = egads.getGeometry(geoC)
  if refG.oclass != egads.CURVE || refG.mtype != egads.CIRCLE
    @warn "getGeometry: CURVE != $(refG.oclass)  CIRCLE != $(refG.oclass)"
  end

  sum((dataC .- refG.reals).^2) > 1.e-14 && @warn "makeGeometry: reals diff $(sum((dataC .- refG.reals).^2)) "

  rangeInfo,periodic = egads.getRange(geoC)

  arc = egads.arcLength(geoC,rangeInfo[1], rangeInfo[2])
  xyz = egads.evaluate(geoC,rangeInfo[1])

  abs(2.0 * pi * radius - arc) > 1.e-14 && @warn "arcLength Circle perimeter  $(2 * pi * radius) != $arc"

  kappa = egads.curvature(geoC, rangeInfo)
  abs(kappa.curvature - 1 / radius) > 1.e-15 && @warn "curvature Circle $(1.0 / radius) != $(kappa.curvature)"

  dataC = [0.0,0.0,0.0, 1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0,radius]
  cylinder = egads.makeGeometry(context, egads.SURFACE, egads.CYLINDRICAL, dataC)

  rInfo = egads.getRange(cylinder)


  mp    = [0.5 * sum(rInfo.range[1,:]), 0.5 * sum(rInfo.range[2,:])]
  kappa = egads.curvature(cylinder, mp)
  (abs(kappa.curvature[1] + 1.0 / radius) > 1.e-15 && abs(kappa.curvature[2] + 1.0 / radius) > 1.e-15) &&
  @warn "curvature Cylinder $(-1.0 / radius), 0 != $(kappa.curvature[1]),$(kappa.curvature[2])"

  isoU     = egads.isoCline(cylinder, 0, mp[1]) ;  isoV = egads.isoCline(cylinder, 1, mp[2])
  isoInfoU = egads.getInfo(isoU) ;  isoInfoV = egads.getInfo(isoV)
  if (isoInfoU.oclass != egads.CURVE || isoInfoU.mtype != egads.LINE ||
      isoInfoV.oclass != egads.CURVE || isoInfoV.mtype != egads.CIRCLE)
     @warn "isoClines U: (CURVE, LINE) != $(isoInfoU.oclass), $(isoInfoU.mtype) V: (CURVE, CIRCLE) != $(isoInfoV.oclass) $(isoInfoV.mtype)"
  end
  sphere = egads.makeSolidBody(context, egads.SPHERE,[0.0,0.0,0.0,radius])

  topF   = egads.getBodyTopos(sphere, egads.FACE)

  infoR  = egads.getRange(topF[1])

  area = egads.getArea(topF[1]; limits = infoR.range)
  abs(2.0 *pi * radius^2 - area) > 1.e-13 && @warn "getArea sphere loop area $(2.0 *pi * radius^2) != $area"

  massP = egads.getMassProperties(sphere)

  if (abs(massP.area - 4 * pi * radius^2) > 1.e-13 || 
      abs(massP.volume - 4 / 3 * pi * radius^3) > 1.e-13||
      sum((massP.CG).^2) > 1.e-13)
     @warn "AREA   $(massP.area)  should be  $(4 * pi * radius^2)  ERR $(abs(massP.area - 4 * pi * radius^2))"
     @warn "VOLUME $(massP.volume) should be $(4 /3 * pi * radius^3) ERR $(abs(massP.volume-4 / 3 * pi * radius^3))"
     @warn "CG     $(massP.CG) should be [0,0,0] ERR $(sum((massP.CG).^2))"
   end

  edges = egads.getBodyTopos(sphere, egads.EDGE)

  rInfo = egads.getRange(edges[2])

  t     = (rInfo.range[1] + rInfo.range[2]) * 0.5
  angle = egads.getWindingAngle(edges[2], t)
  angle != 180.0 && @warn "getWindingAngle 180 != $(angle.angle)"
  return true
end

function test_evaluate()
  wireBody = make_wire(context, [[0,0,0]  [1,0,0] [1,1,0] [0,1,0]], egads.BODY)
  info     = egads.getInfo(wireBody)
  nodes    = egads.getBodyTopos(wireBody,egads.NODE)
  eval     = egads.evaluate(nodes[1])
  eval    != [0.0,0.0,0.0] && @warn "evaluate node returns wrong value: $eval != [0.0,0.0,0.0]"
  edges    = egads.getBodyTopos(wireBody,egads.EDGE)

  rE    = egads.getRange(edges[1])
  eval  = egads.evaluate(edges[1],0.5 * sum(rE.range))
  inve  = egads.invEvaluate(edges[1], [0.5, 0., 0.])

  eval.X      != [0.5,0.0,0.0] && @warn "evaluate edge returns wrong value: $(eval.X) != [0.5,0.0,0.0]"
  eval.dX     != [1.0,0.0,0.0] && @warn "evaluate edge returns wrong 1st derivative: $(eval.dX) != [1.0,0.0,0.0]"
  eval.ddX    != [0.0,0.0,0.0] && @warn "evaluate edge returns wrong 2nd derivative: $(eval.ddX)!= [0.0,0.0,0.0]"
  inve.result != [0.5,0.0,0.0] && @warn "invEvaluate edge returns wrong value: $(inve.result) != [0.5,0.0,0.0]"

  faces = egads.getBodyTopos(wireBody,egads.FACE)
  rF    = egads.getRange(faces[1])
  uv    = 0.5 .* [sum(rF.range[1,:]) ,sum(rF.range[2,:]) ]

  eval  = egads.evaluate(faces[1], uv)

  inve  = egads.invEvaluate(faces[1], [0.5, 0.5, 0.])
  !(isapprox(inve.result,[0.5,0.5,0.0],atol = 1.e-10)) &&
    @warn "invEvaluate face returns wrong result: $(inve.result) != [0.5,0.5,0.0]"

    !(isapprox(eval.X,[0.5,0.5,0.0],atol = 1.e-10)) &&
      @warn "evaluate face returns wrong value: $(eval.X) != [0.5,0.5,0.0]"

  dv = cos(pi*0.25)
  !(isapprox(eval.dX[1],[-dv,-dv,0], atol = 1.e-10)) &&
    @warn "evaluate face returns wrong 1st derivative in u: $(eval.dX) != [1.0,0.0,0.0]"
  !(isapprox(eval.dX[2],[dv,-dv,0], atol = 1.e-10)) &&
    @warn "evaluate face returns wrong 1st derivative in v: $(eval.dX) != [0.0,1.0,0.0]"
 !(isapprox(eval.ddX[1],[0,0,0], atol = 1.e-10)) &&
    @warn "evaluate face returns wrong 2nd derivative (duu): $(eval.dX) != [0.0,0.0,0.0]"
  !(isapprox(eval.ddX[2],[0,0,0], atol = 1.e-10)) &&
    @warn "evaluate face returns wrong cross derivatives (duv): $(eval.dX) != [0.0,0.0,0.0]"
  !(isapprox(eval.ddX[3],[0,0,0], atol = 1.e-10)) &&
    @warn "evaluate face returns wrong 2nd derivative (dvv): $(eval.dX) != [0.0,0.0,0.0]"

  !(isapprox(inve.result,[0.5,0.5,0.0],atol = 1.e-10)) &&
    @warn " egads.invEvaluate: face returns wrong value: $(inve.result) != [0.5,0.5,0.0]"
  !(isapprox(inve.params,uv,atol = 1.e-10)) &&
    @warn " egads.invEvaluate: face returns wrong param: $(inve.params) != $uv"
  return true
end

function test_otherCurve()
  box   = egads.makeSolidBody(context, egads.BOX, [1,0,0,-1,1,1,2])
  edges = egads.getBodyTopos(box, egads.EDGE)
  idx   = [egads.indexBodyTopo(box,edges[i])  for i = 1:length(edges)]
  any(idx .== false) && @warn "indexBodyTopo: edges wrong id"
  edge1 = egads.objectBodyTopo(box,egads.EDGE,1)
  !egads.isSame(edge1, edges[1]) && @warn "isSame != same edge"

  faces = egads.getBodyTopos(box, egads.FACE ; ref = edges[1])
  pcurv = egads.otherCurve(faces[1], edges[1])
  data = egads.getInfo(pcurv)
  (data.oclass != egads.PCURVE || data.mtype != egads.LINE) &&
    @warn "otherCurve (PCURVE, LINE) != $(data.oclass), $(data.mtype)"
  return true
end

function test_getEdgeUV()

  wireBody = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.BODY)
  faces    = egads.getBodyTopos(wireBody,egads.FACE)
  edges    = egads.getBodyTopos(wireBody,egads.EDGE)
  rInfo    = egads.getRange(edges[1])
  t        = 0.5 * sum(rInfo.range)
  uvs      = egads.getEdgeUV(faces[1], edges[1], 1,t)

  xyz1  = egads.evaluate(faces[1], uvs)
  xyz2  = egads.evaluate(edges[1], t)
    sum((xyz1.X .- xyz2.X).^2) > 1.e-14 && @warn "getEdgeUV $(uvs.uv) in face != $(t) IN EDGE"
  return true
end

function test_match_sew_replace_topology()

  box1 = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context,egads.BOX, [1.0,0.0,0.0,1.0,1.0,1.0])

  matchB = egads.matchBodyEdges(box1,box2)
  if matchB != [[5, 1]  [6, 2]  [7, 3] [8, 4]]
    @warn "matchBodyEdges [[5, 1]  [6, 2]  [7, 3] [8, 4]] != $matchB"
  end

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  matchB = egads.matchBodyFaces(box1,box2)
  (matchB != [2,1]) && @warn "matchBodyEdges [2,1] != $matchB"

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
  nF = length(rFaces)
  nF != 4 && @warn "replaceFaces new body has $nF != 4"

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
  loop   = egads.makeLoop!([edgeC,edgeL])

  face  = egads.makeFace(loop, egads.SFORWARD)
  uv,x0 = egads.invEvaluate(face,[ radius,0.0,0.0])
  stat  = egads.inFace(face,uv)
  return true
end

function test_inTopology()
  box   = egads.makeSolidBody(context, egads.BOX, [0,0,0, 1,1,1])
  faces = egads.getBodyTopos(box,egads.FACE)
  info  = egads.getRange(faces[1])
  bbox  = egads.getBoundingBox(faces[1])

  !egads.inTopology(faces[1],0.5 .* [sum(bbox[1]),sum(bbox[2]),sum(bbox[3])]) &&
    @warn "inTopology center bbox not in face !!! "
  !egads.inFace(faces[1],0.5 .* [sum(info.range[1]),sum(info.range[2])]) &&
    @warn "inFace center range not in face !!! "
  return true
end

function test_booleanOperations()

  box1    = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2    = egads.makeSolidBody(context,egads.BOX, [-0.5,-0.5,-0.5,1.0,1.0,1.0])

  bt      = [egads.SUBTRACTION, egads.INTERSECTION, egads.FUSION]
  boolbox = [egads.generalBoolean(box1, box2, bt[k]) for k = 1:3]
  return true
end

function test_intersection()

  wBody  = make_wire(context, [[0,0,0] [1,0,0] [1,1,0] [0,1,0]], egads.BODY)
  box    = egads.makeSolidBody(context,egads.BOX, [0.25,0.25,-0.25, 0.5,0.5,0.5])
  inter  = egads.intersection(box,wBody)
  box2   = egads.imprintBody( box, inter.pairs)
  faces1 = egads.getBodyTopos(box2, egads.FACE)
  length(faces1) != 10 && @warn "imprintBody topos 10 != $(length(faces1))"
  return true
end

function test_filletBody()

  box1   = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  edges1 = egads.getBodyTopos(box1, egads.EDGE)
  box4   = egads.filletBody(box1, edges1, 0.1)
  info   = egads.getInfo(box4.result)
  (info.oclass != egads.BODY || info.mtype != egads.SOLIDBODY) &&
    @warn "filletBody: (BODY,SOLIDBODY) != $(info.oclass),$(info.mtype)"
  return true
end

function test_fuseSheets()
  box1    = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  box2    = egads.makeSolidBody(context,egads.BOX, [-0.5,-0.5,-0.5,1.0,1.0,1.0])

  shell1 = egads.getBodyTopos(box1, egads.SHELL)
  shell2 = egads.getBodyTopos(box2, egads.SHELL)

  sheet1 = egads.makeTopology!(context,egads.BODY;mtype = egads.SHEETBODY, children = shell1[1])
  sheet2 = egads.makeTopology!(context,egads.BODY;mtype = egads.SHEETBODY, children = shell2[1])

  fuse = egads.fuseSheets(sheet1, sheet2)
  return true
end

function test_chamferHollow()

  box1     = egads.makeSolidBody(context,egads.BOX, [0.0,0.0,0.0,1.0,1.0,1.0])
  edges1   = egads.getBodyTopos(box1,egads.EDGE)
  b1_faces = [egads.getBodyTopos(box1, egads.FACE ; ref = k)[1] for k in edges1]
  cham     = egads.chamferBody(box1, edges1, b1_faces, 0.1, 0.15)
  info     = egads.getInfo(cham.result)
  (info.oclass != egads.BODY || info.mtype != egads.SOLIDBODY) &&
    @warn "chamferBody (BODY, SOLIDBODY) != $(info.oclass),$(info.mtype) !"


  holl = egads.hollowBody(box1, nothing, 0.1, 0)
  info = egads.getInfo(holl.result)
  (info.oclass != egads.BODY || info.mtype != egads.SOLIDBODY) &&
    @warn "hollowBody (BODY, SOLIDBODY) != $(info.oclass),$(info.mtype) !"

  face = egads.getBodyTopos(box1, egads.FACE)[1]
  holl2 = egads.hollowBody(box1, face, 0.1, 0)

  info = egads.getInfo(holl2.result)
  (info.oclass != egads.BODY || info.mtype != egads.SOLIDBODY) &&
    @warn "hollowBody (BODY, SOLIDBODY) != $(info.oclass),$(info.mtype) !"
  return true
end

function test_sweepExtrude()

  faceBODY = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.SURFACE)

  xyz   = [ [0.0,0.0,0.0] [0.0,0.0,1.0]]
  nodes = [egads.makeTopology!(context, egads.NODE; reals = xyz[:,k]) for k = 1:2]

  rdata = [0.0, 0.0,0.0, 0.0,0.0,1.0]      # Line data (point and direction)
  line  = egads.makeGeometry(context, egads.CURVE, egads.LINE, rdata)
  edge  = egads.makeTopology!(context, egads.EDGE ; mtype = egads.TWONODE, geom=line, children=nodes, reals=[0.0,1.0])


  sweep  = egads.sweep(faceBODY, edge, 0)
  sFaces = egads.getBodyTopos(sweep,egads.FACE)
  lf = length(sFaces)
  lf != 6 && @warn "sweep body topos 6 != $lf"

  extrude = egads.extrude(faceBODY, 1.0, [0.0,0.0,1.0])
  eFaces  = egads.getBodyTopos(extrude, egads.FACE)
  le      = length(eFaces)
  le != 6 && @warn "extrude body topos 6 != $le"
  return true
end

function test_effective()

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
  [push!(secs, egads.copyObject(secs[i-1] ,xform)) for i = 2:4]

  ruled = egads.ruled(secs)
  faces = egads.getBodyTopos(ruled, egads.FACE)
  lf    = length(faces)
  12 != lf && @warn "ruled surface has $lf!= 12 faces !"
  [egads.attributeAdd!(faces[k], "faceMap",1) for k = 1:2]

  tess  = egads.makeTessBody(ruled, [0.1, 0.1, 15])

  ebody = egads.initEBody(tess, 10.0)

  eface = egads.makeEFace(ebody, faces[10:11])
  info = egads.getInfo(eface)
  (info.oclass != egads.EFACE || info.mtype != egads.SFORWARD) &&
    @warn " egads.makeEFace (EFACE, SFORWARD) != $(info.oclass), $(info.mtype) !"

  efaces = egads.makeAttrEFaces(ebody,"faceMap")
  le     = length(efaces)
  1 != le && @warn "makeAttrEFaces : # faces 1 != $le !"

  info = egads.getInfo(efaces[1])
  (info.oclass != egads.EFACE || info.mtype != egads.SFORWARD) &&
    @warn " egads.makeAttrEFaces (EFACE, SFORWARD) != $(info.oclass), $(info.mtype) !"
  egads.finishEBody!(ebody)
  bTopo = egads.getTopology(ebody)

  bl    = length(bTopo.children)
  bc    = bTopo.oclass
  bg    = bTopo.geom

  (bc != egads.EBODY || bl != 1 || bg.obj != ruled.obj) &&
    @warn " egads.getTopology EBODY != $bc, 1 != $bl , $ruled) != $bg"

  for shell in bTopo.children
    sTopo = egads.getTopology(shell)
    ls =  length(sTopo.children)
    if sTopo.oclass != egads.ESHELL ||  ls != 10
      @warn "getTopology (ESHELL, 10) != $(sTopo.oclass),$ls"
      continue
    end
    for face in sTopo.children
      fTopo = egads.getTopology(face)
      lf = length(fTopo.children)
      if fTopo.oclass != egads.EFACE || lf != 1
        @warn "getTopology (EFACE, 1) != $(fTopo.oclass),$lf"
        continue
      end
      for loop in fTopo.children
        lTopo = egads.getTopology(loop)
        lTopo.oclass != egads.ELOOPX &&
          @warn "getTopology ELOOPX != $(lTopo.oclass)"
      end
    end
  end

  # check getting the original FACES
  faces  = egads.getBodyTopos(ebody, egads.FACE)
  lf     = length(faces)
  12    != lf && @warn "ruled surface has $lf != 12 faces "

  efaces = egads.getBodyTopos(ebody, egads.EFACE)
  le     = length(efaces)
  10    != length(efaces) && @warn "ruled surface has $le != 10 faces "

  eUV    = egads.getRange(efaces[1])
  fMap   = egads.effectiveMap(efaces[1], 0.5 .* [sum(eUV.range[1,:]), sum(eUV.range[2,:])])

  eedges = egads.getBodyTopos(ebody, egads.EEDGE)
  eList  = egads.effectiveEdgeList(eedges[1])

  eT   = egads.getRange(eedges[1])
  eMap = egads.effectiveMap(eedges[1], 0.5 .* sum(eT.range))
  return true
end

function test_blend_and_ruled()

  # first section is a NODE
  node1 = egads.makeTopology!(context, egads.NODE ; reals=[0.5,0.5,-1])

  # FaceBody for 2nd section
  face1 = make_wire(context, [[0,0,0] [1,0,0] [1,1,0] [0,1,0]], egads.SURFACE)
  mat1     = [ [1.0,0.0,0.0,0.0] [0.0,1.0,0.0,0.0] [0.0,0.0,1.0,1.0]]
  xform    = egads.makeTransform(context, mat1)

  # Translate 2nd section to make 3rd
  face2 = egads.copyObject(face1,xform)

  # Finish with a NODE
  node2 = egads.makeTopology!(context, egads.NODE ; reals=[0.5,0.5,2])

  # sections for blend
  secs = [node1,face1,face2,node2]
  blend1 = egads.blend(secs)
  data   = egads.getInfo(blend1)
  (data.oclass != egads.BODY || data.mtype != egads.SOLIDBODY) &&
    @warn "blend (BODY,SOLIDBODY) != $(data.oclass),$(data.mtype) !"

  rc = [0.2, 1, 0, 0,  0.4, 0, 1, 0]
  blend2 = egads.blend(secs; rc1=rc)
  data = egads.getInfo(blend2)
  (data.oclass != egads.BODY || data.mtype != egads.SOLIDBODY) &&
    @warn "blend (BODY,SOLIDBODY) != $(data.oclass),$(data.mtype) !"

  blend3 = egads.blend(secs; rc1=rc, rc2=rc)
  data = egads.getInfo(blend3)
  (data.oclass != egads.BODY || data.mtype != egads.SOLIDBODY) &&
    @warn "blend (BODY,SOLIDBODY) != $(data.oclass),$(data.mtype) !"

  blend4 = egads.blend(secs; rc2=rc)
  data = egads.getInfo(blend4)
  (data.oclass != egads.BODY || data.mtype != egads.SOLIDBODY) &&
    @warn "blend (BODY,SOLIDBODY) != $(data.oclass),$(data.mtype) !"

  ruled = egads.ruled(secs)
  data = egads.getInfo(ruled)
  (data.oclass != egads.BODY || data.mtype != egads.SOLIDBODY) &&
    @warn "ruled (BODY,SOLIDBODY) != $(data.oclass),$(data.mtype) !"
  return true
end

function test_rotate()
  faceBODY = make_wire(context, [[0.0,0.0,0.5] [1.0,0.0,0.5] [1.0,1.0,0.5] [0.0,1.0,0.5] ],egads.SURFACE)

  body    = egads.rotate(faceBODY, 180, [[0,-1,0]  [1,0,0]])

  faces = egads.getBodyTopos(body,egads.FACE)
  lf = length(faces)
  lf != 6 && @warn "rotate body faces 6 != $le"
  return true
end

function test_mapBody()

  box1 = egads.makeSolidBody(context, egads.BOX,[0.0,0.0,0.0,1.0,1.0,1.0])
  box2 = egads.makeSolidBody(context, egads.BOX,[0.1,0.1,0.1,1.0,1.0,1.0])

  faces1 = egads.getBodyTopos(box1, egads.FACE)
  faces2 = egads.getBodyTopos(box2, egads.FACE)

  egads.attributeAdd!(faces1[2], "faceMap", [1])
  egads.attributeAdd!(faces2[1], "faceMap", [1])

  mapB  = egads.mapBody(box1,box2, "faceMap")
  tess1 = egads.makeTessBody(box1,[0.1, 0.1, 15.])

  tess2   = egads.mapTessBody(tess1, box2)
  tessF   = egads.getTessFace(tess1, 1)
  ifaces  = fill(1, Int(floor(length(tessF.uv[1,:])/2)))
  locTess = egads.locateTessBody(tess2, ifaces, tessF.uv)
  return true
end

function test_fitTriangles()

  box  = egads.makeSolidBody(context, egads.BOX,[0.0,0.0,0.0,1.0,1.0,1.0])
  tess = egads.makeTessBody(box,[0.05, 0.1, 30.])
  tessData = egads.getTessFace(tess, 1)

  triBox = egads.fitTriangles(context, tessData.xyz, tessData.tris; tol = 1.e-8)
  info   = egads.getInfo(triBox)

  (info.mtype != egads.BSPLINE || info.oclass != egads.SURFACE) &&
    @warn "fitTriangles (BSPLINE, SURFACE) != $(info.mtype),$(info.oclass)"

  triBox2 = egads.fitTriangles(context, tessData.xyz, tessData.tris; tessData.tric, tol = 1.e-8)
  info    = egads.getInfo(triBox2)
  (info.mtype != egads.BSPLINE || info.oclass != egads.SURFACE) &&
    @warn "fitTriangles (BSPLINE, SURFACE) != $(info.mtype),$(info.oclass)"
  return true
end

function test_quadTess()
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
  if length(fList) != 1 || fList[1] != 1
    @warn "getTessQuad should be 1 face with index 1 but got $fList"
  end

  plane = egads.makeGeometry(context,egads.SURFACE,egads.PLANE,[0.0,0.0,0.0, 1.0,0.0,0.0 ,0.0,1.0,0.0])
  pTess = egads.makeTessGeom( plane, [0.0,1.0, 0.0,1.0], [5,5])

  planeT = egads.getTessGeom(pTess)
  if planeT.sizes != [5,5]
    @warn "getTessGeom should have given [5,5] sizes but got $(planeT.sizes) !"
  end
  return true
end

function test_tessBody()
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
      tdiff + xyzdiff > 1.e-14 && @warn "Set tess2  $i gave different $t, $xyz in edges"
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

function test_modifyTess()
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
  if (eA > 0.1 || eV > 0.1||
      sum((massP.CG).^2) >0.1)
       @warn " LINE    $(@__LINE__) egads.tessMassProps"
       @warn " AREA    $massP.area   != $(4 * pi * radius^2) err $eA"
       @warn " VOLUME  $massP.volume != $(4 /3 * pi * radius^3) err $eV"
       @warn " CG  $massP.CG [0,0,0] $(sum((massP.CG).^2))"
   end
  return true
end


function test_memo()
  ctxt  = egads.Ccontext()
  level = egads.setOutLevel(ctxt,2)
  face = make_wire(ctxt, [[0.0,0,0] [1.0,0,0] [1.0,1.0,0] [0,1.0,0]], egads.SURFACE)
  GC.gc()
  info = egads.getInfo(context)
  if info.next.obj != C_NULL
    @warn " FINALIZERS NOT WORKING. CONTEXT SHOULD BE CLEAN !! "
    return false
  else return true
  end
end


@time @testset " Julia egads API " begin
  # GET ALL FUNCTION TESTS
  matches = []
  open(ENV["jlEGADS"]*"/test/runtests.jl","r") do f
    ext = "test_"
    while ! eof(f)
        # read a new / next line for every iteration
        a   = readline(f)
        idJ = findfirst(ext, a)
        (idJ == nothing || length(idJ) == 0 ) && continue
        push!(matches, a[idJ[1]:end])
    end
  end
  fnames = matches[1:end-1]
  ### NOW RUN
  global context = egads.Ccontext(name  = "ctxt_glob")
  level = egads.setOutLevel(context,1)
  for j = 1:length(fnames)
    x =  eval(Meta.parse(fnames[j]))
    @test x
  end
  rm(loc_folder*"/skinning.egads")
  rm(loc_folder*"/faceBody.egads")
end
