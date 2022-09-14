#
#     EGADS: Engineering Geometry Aircraft Design System
#
#            Python Construction Example
#
#     Copyright 2011-2022, Massachusetts Institute of Technology
#     Licensed under The GNU Lesser General Public License, version 2.1
#     See http://www.opensource.org/licenses/lgpl-2.1.php
#
#
#

from pyEGADS import egads
import math
import sys

TWOPI = 2*math.pi

def periodicSeam(eedge,    # Edge associated with seam
                 sense):   # sense of PCurve

    context = eedge.context

    range, per = eedge.getRange()

    # set up u and v at ends
    data = [TWOPI, range[0],
                0,    sense]
    if sense == egads.SREVERSE: data[1] = range[1]

    # return (linear) PCurve
    return context.makeGeometry(egads.PCURVE, egads.LINE, data)

# set the parameter values
width  =  5.0
minrad =  8.0
maxrad = 12.0
fillet =  2.0
thick  =  0.5
bolts  =  5
crad   =  5.0
brad   =  1.0

# define a context
context = egads.Context()

# Node locations
node1 = [-minrad, 0.0, -width / 2.0]
node2 = [-minrad, 0.0,  width / 2.0]
node3 = [-maxrad, 0.0,  width / 2.0]
node4 = [-maxrad, 0.0, -width / 2.0]
node5 = [ minrad, 0.0, -width / 2.0]
node6 = [ maxrad, 0.0, -width / 2.0]
node7 = [ maxrad, 0.0,  width / 2.0]
node8 = [ minrad, 0.0,  width / 2.0]

# make Nodes */
enodes = [None]*8
enodes[0] = context.makeTopology(egads.NODE, reals=node1)
enodes[1] = context.makeTopology(egads.NODE, reals=node2)
enodes[2] = context.makeTopology(egads.NODE, reals=node3)
enodes[3] = context.makeTopology(egads.NODE, reals=node4)
enodes[4] = context.makeTopology(egads.NODE, reals=node5)
enodes[5] = context.makeTopology(egads.NODE, reals=node6)
enodes[6] = context.makeTopology(egads.NODE, reals=node7)
enodes[7] = context.makeTopology(egads.NODE, reals=node8)

ecurve = [None]*16
eedges = [None]*16
trange = [None]*2

# make (linear) Edge 1
data = [None]*6
data[0] = node1[0]
data[1] = node1[1]
data[2] = node1[2]
data[3] = node2[0] - node1[0]
data[4] = node2[1] - node1[1]
data[5] = node2[2] - node1[2]
ecurve[0] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[0].invEvaluate(node1)
trange[1], xyz = ecurve[0].invEvaluate(node2)

elist = [enodes[0], enodes[1]]
eedges[0] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=ecurve[0], reals=trange, children=elist)

# make (linear) Edge 2
data = [None]*6
data[0] = node2[0]
data[1] = node2[1]
data[2] = node2[2]
data[3] = node3[0] - node2[0]
data[4] = node3[1] - node2[1]
data[5] = node3[2] - node2[2]
ecurve[1] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[1].invEvaluate(node2)
trange[1], xyz = ecurve[1].invEvaluate(node3)

elist = [enodes[1], enodes[2]]
eedges[1] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[1], trange, elist)

# make (linear) Edge 3
data = [None]*6
data[0] = node3[0]
data[1] = node3[1]
data[2] = node3[2]
data[3] = node4[0] - node3[0]
data[4] = node4[1] - node3[1]
data[5] = node4[2] - node3[2]
ecurve[2] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[2].invEvaluate(node3)
trange[1], xyz = ecurve[2].invEvaluate(node4)

elist = [enodes[2], enodes[3]]
eedges[2] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[2], trange, elist)

# make (linear) Edge 4
data = [None]*6
data[0] = node4[0]
data[1] = node4[1]
data[2] = node4[2]
data[3] = node1[0] - node4[0]
data[4] = node1[1] - node4[1]
data[5] = node1[2] - node4[2]
ecurve[3] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[3].invEvaluate(node4)
trange[1], xyz = ecurve[3].invEvaluate(node1)

elist = [enodes[3], enodes[0]]
eedges[3] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[3], trange, elist)

# make (linear) Edge 5
data = [None]*6
data[0] = node5[0]
data[1] = node5[1]
data[2] = node5[2]
data[3] = node6[0] - node5[0]
data[4] = node6[1] - node5[1]
data[5] = node6[2] - node5[2]
ecurve[4] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[4].invEvaluate(node5)
trange[1], xyz = ecurve[4].invEvaluate(node6)

elist = [enodes[4], enodes[5]]
eedges[4] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[4], trange, elist)

# make (linear) Edge 6
data = [None]*6
data[0] = node6[0]
data[1] = node6[1]
data[2] = node6[2]
data[3] = node7[0] - node6[0]
data[4] = node7[1] - node6[1]
data[5] = node7[2] - node6[2]
ecurve[5] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[5].invEvaluate(node6)
trange[1], xyz = ecurve[5].invEvaluate(node7)

elist = [enodes[5], enodes[6]]
eedges[5] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[5], trange, elist)

# make (linear) Edge 7
data = [None]*6
data[0] = node7[0]
data[1] = node7[1]
data[2] = node7[2]
data[3] = node8[0] - node7[0]
data[4] = node8[1] - node7[1]
data[5] = node8[2] - node7[2]
ecurve[6] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[6].invEvaluate(node7)
trange[1], xyz = ecurve[6].invEvaluate(node8)

elist = [enodes[6], enodes[7]]
eedges[6] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[6], trange, elist)

# make (linear) Edge 8
data = [None]*6
data[0] = node8[0]
data[1] = node8[1]
data[2] = node8[2]
data[3] = node5[0] - node8[0]
data[4] = node5[1] - node8[1]
data[5] = node5[2] - node8[2]
ecurve[7] = context.makeGeometry(egads.CURVE, egads.LINE, data)

trange[0], xyz = ecurve[7].invEvaluate(node8)
trange[1], xyz = ecurve[7].invEvaluate(node5)

elist = [enodes[7], enodes[4]]
eedges[7] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[7], trange, elist)

# data used in creating the arcs
axis1 = [1, 0, 0]
axis2 = [0, 1, 0]
axis3 = [0, 0, 1]

cent1 = [0, 0, -width / 2.0]
cent2 = [0, 0,  width / 2.0]

# make (circular) Edge 9
data = [cent1[0], cent1[1], cent1[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2], minrad]
ecurve[8] = context.makeGeometry(egads.CURVE, egads.CIRCLE, data)

trange[0], xyz = ecurve[8].invEvaluate(node5)
trange[1], xyz = ecurve[8].invEvaluate(node1)

if trange[0] > trange[1]: trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[4], enodes[0]]
eedges[8] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[8], trange, elist)

# make (circular) Edge 10
data = [cent2[0], cent2[1], cent2[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2], minrad]
ecurve[9] = context.makeGeometry(egads.CURVE, egads.CIRCLE, data)

trange[0], xyz = ecurve[9].invEvaluate(node8)
trange[1], xyz = ecurve[9].invEvaluate(node2)

if trange[0] > trange[1]: trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[7], enodes[1]]
eedges[9] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[9], trange, elist)

# make (circular) Edge 11
data = [cent1[0], cent1[1], cent1[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2], maxrad]
ecurve[10] = context.makeGeometry(egads.CURVE, egads.CIRCLE, data)

trange[0], xyz = ecurve[10].invEvaluate(node6)
trange[1], xyz = ecurve[10].invEvaluate(node4)

if trange[0] > trange[1]: trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[5], enodes[3]]
eedges[10] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[10], trange, elist)

# make (circular) Edge 12
data = [cent2[0], cent2[1], cent2[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2], maxrad]
ecurve[11] = context.makeGeometry(egads.CURVE, egads.CIRCLE, data)

trange[0], xyz = ecurve[11].invEvaluate(node7)
trange[1], xyz = ecurve[11].invEvaluate(node3)

if trange[0] > trange[1]: trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[6], enodes[2]]
eedges[11] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[11], trange, elist)

# make the outer cylindrical surface
esurface = [None]*4
data = [cent1[0], cent1[1], cent1[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2],
        axis3[0], axis3[1], axis3[2], maxrad]
esurface[0] = context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, data)
esurface[2] = context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, data)

# make the inner cylindrical surface
data = [cent1[0], cent1[1], cent1[2],
        axis1[0], axis1[1], axis1[2],
        axis2[0], axis2[1], axis2[2],
        axis3[0], axis3[1], axis3[2], minrad]
esurface[1] = context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, data)
esurface[3] = context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, data)

efaces = [None]*8

# make (planar) Face 1
sense = [egads.SFORWARD, egads.SREVERSE, egads.SFORWARD, egads.SFORWARD]
elist = [     eedges[3],      eedges[8],      eedges[4],     eedges[10]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, children=elist, senses=sense)

efaces[0] = eloop.makeFace(egads.SFORWARD)

# make (planar) Face 2
sense = [egads.SFORWARD, egads.SREVERSE, egads.SFORWARD, egads.SFORWARD]
elist = [     eedges[1],     eedges[11],      eedges[6],      eedges[9]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, children=elist, senses=sense)

efaces[1] = eloop.makeFace(egads.SFORWARD)

# make (cylindrical) Face 3
epcurve = [esurface[0].otherCurve(ecurve[ 2]),
           esurface[0].otherCurve(ecurve[10]),
           esurface[0].otherCurve(ecurve[ 5]),
           esurface[0].otherCurve(ecurve[11])]

sense = [egads.SFORWARD, egads.SREVERSE, egads.SFORWARD, egads.SFORWARD]
elist = [     eedges[2],     eedges[10],      eedges[5],     eedges[11],
             epcurve[0],    epcurve[ 1],     epcurve[2],    epcurve[ 3]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, esurface[0], children=elist, senses=sense)

efaces[2] = context.makeTopology(egads.FACE, egads.SREVERSE, esurface[0], children=[eloop], senses=[egads.SFORWARD])

# make (cylindrical) Face 4
epcurve = [esurface[1].otherCurve(ecurve[0]),
           esurface[1].otherCurve(ecurve[9]),
           esurface[1].otherCurve(ecurve[7]),
           esurface[1].otherCurve(ecurve[8])]

sense = [egads.SFORWARD, egads.SREVERSE, egads.SFORWARD, egads.SFORWARD]
elist = [     eedges[0],      eedges[9],      eedges[7],      eedges[8],
             epcurve[0],     epcurve[1],     epcurve[2],     epcurve[3]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, esurface[1], children=elist, senses=sense)

efaces[3] = context.makeTopology(egads.FACE, egads.SFORWARD, esurface[1], children=[eloop], senses=[egads.SFORWARD])

# make (circular) Edge 13
ecurve[12] = ecurve[8]   # reuse

trange[0], xyz = ecurve[12].invEvaluate(node1)
trange[1], xyz = ecurve[12].invEvaluate(node5)

if (trange[0] > trange[1]): trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[0], enodes[4]]
eedges[12] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[12], trange, elist)

# make (circular) Edge 14
ecurve[13] = ecurve[9]   # reuse

trange[0], xyz = ecurve[13].invEvaluate(node2)
trange[1], xyz = ecurve[13].invEvaluate(node8)

if (trange[0] > trange[1]): trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[1], enodes[7]]
eedges[13] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[13], trange, elist)

# make (circular) Edge 15
ecurve[14] = ecurve[10]

trange[0], xyz = ecurve[13].invEvaluate(node4)
trange[1], xyz = ecurve[13].invEvaluate(node6)

if (trange[0] > trange[1]): trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[3], enodes[5]]
eedges[14] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[14], trange, elist)

# make (circular) Edge 16
ecurve[15] = ecurve[11]

trange[0], xyz = ecurve[15].invEvaluate(node4)
trange[1], xyz = ecurve[15].invEvaluate(node6)

if (trange[0] > trange[1]): trange[1] += TWOPI   # ensure trange[1] > trange[0]

elist = [enodes[2], enodes[6]]
eedges[15] = context.makeTopology(egads.EDGE, egads.TWONODE, ecurve[15], trange, elist)

# make (planar) Face 5
sense = [egads.SFORWARD, egads.SREVERSE, egads.SREVERSE, egads.SREVERSE]
elist = [    eedges[14],      eedges[4],     eedges[12],      eedges[3]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, children=elist, senses=sense)

efaces[4] = eloop.makeFace(egads.SFORWARD)

# make (planar) Face 6
sense = [egads.SFORWARD, egads.SREVERSE, egads.SFORWARD, egads.SFORWARD]
elist = [     eedges[6],     eedges[13],      eedges[1],     eedges[15]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, children=elist, senses=sense)

efaces[5] = eloop.makeFace(egads.SFORWARD)

# make (cylindrical) Face 7
epcurve = [esurface[2].otherCurve(ecurve[ 2]),
           esurface[2].otherCurve(ecurve[14]),
# otherCurve cannot be used along periodic seam
#          esurface[2].otherCurve(ecurve[5]),
           periodicSeam(eedges[5], egads.SFORWARD),
           esurface[2].otherCurve(ecurve[15])]

sense = [egads.SFORWARD, egads.SFORWARD, egads.SFORWARD, egads.SREVERSE]
elist = [     eedges[2],     eedges[14],      eedges[5],     eedges[15],
             epcurve[0],    epcurve[ 1],     epcurve[2],    epcurve[ 3]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, esurface[2], children=elist, senses=sense)

efaces[6] = context.makeTopology(egads.FACE, egads.SFORWARD, esurface[2], children=[eloop], senses=[egads.SFORWARD])

# make (cylindrical) Face 8
epcurve = [esurface[3].otherCurve(ecurve[ 0]),
           esurface[3].otherCurve(ecurve[13]),
# otherCurve cannot be used along periodic seam
#          esurface[3].otherCurve(ecurve[ 7]),
           periodicSeam(eedges[7], egads.SREVERSE),
           esurface[3].otherCurve(ecurve[12])]

sense = [egads.SFORWARD, egads.SFORWARD, egads.SFORWARD, egads.SREVERSE]
elist = [     eedges[0],     eedges[13],      eedges[7],     eedges[12],
             epcurve[0],    epcurve[ 1],     epcurve[2],    epcurve[ 3]]
eloop = context.makeTopology(egads.LOOP, egads.CLOSED, esurface[3], children=elist, senses=sense)

efaces[7] = context.makeTopology(egads.FACE, egads.SREVERSE, esurface[3], children=[eloop], senses=[egads.SFORWARD])

# make the shell and initial Body
eshell = context.makeTopology(egads.SHELL, egads.CLOSED, children=efaces)
ebody  = context.makeTopology(egads.BODY, egads.SOLIDBODY, children=[eshell])

# add fillets if desired (result is ebody2)
if (fillet > 0):
    elist = [eedges[10], eedges[11], eedges[14], eedges[15]]
    ebody, maps = ebody.filletBody(elist, fillet)

# add wheel if desired
if (thick > 0):
    data = [0, 0,  thick / 2.0,
            0, 0, -thick / 2.0,
          (minrad + maxrad) / 2.0]
    ecylbody = context.makeSolidBody(egads.CYLINDER, data)

    emodel = ebody.generalBoolean(ecylbody, egads.FUSION)

    oclass, mtype, eref, data, echilds, senses = emodel.getTopology()

    if (oclass != egads.MODEL or len(echilds) != 1):
        print("You didn't input a model or are returning more than one body ochild = " + str(oclass) + " nchild = " + str(len(echilds)))
        sys.exit(-999)

    ebody = echilds[0].copyObject()

    # add bolt holes
    for i  in range(bolts):
        data = [None]*7
        data[0] = crad * math.cos(i * (TWOPI / bolts))
        data[1] = crad * math.sin(i * (TWOPI / bolts))
        data[2] =  thick / 2.0
        data[3] = crad * math.cos(i * (TWOPI / bolts))
        data[4] = crad * math.sin(i * (TWOPI / bolts))
        data[5] = -thick / 2.0
        data[6] = brad

        ecylbody = context.makeSolidBody(egads.CYLINDER, data)

        emodel = ebody.generalBoolean(ecylbody, egads.SUBTRACTION)

        oclass, mtype, eref, data, echilds, senses = emodel.getTopology()

        if (oclass != egads.MODEL or len(echilds) != 1):
            print("You didn't input a model or are returning more than one body ochild = " + str(oclass) + " nchild = " + str(len(echilds)))
            sys.exit(-999)

        ebody = echilds[0].copyObject()

# make and dump the model
emodel = context.makeTopology(egads.MODEL, children=[ebody])
emodel.saveModel("tire.egads", overwrite=True)
