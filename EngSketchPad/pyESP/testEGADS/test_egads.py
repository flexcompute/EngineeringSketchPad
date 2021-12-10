import unittest
from pyEGADS import egads
import os
import math
import sys

class TestEGADS(unittest.TestCase):

    # setUp a new context for each test
    def setUp(self):
        self.context = egads.Context()
        self.context.setOutLevel(0)

    # make sure the context is clean before it is closed
    def tearDown(self):
        oclass, mtype, topObj, prev, next = self.context.getInfo()
        self.assertFalse( next , "Context is not properly clean!")
        del self.context

#==============================================================================
    def test_revision(self):
        imajor, iminor, OCCrev = egads.revision()
        self.assertIsInstance(imajor,int)
        self.assertIsInstance(iminor,int)
        self.assertIsInstance(OCCrev,str)

#==============================================================================
    def test_OpenClose(self):
        # Test open/close a context

        # Test creating a reference to a context c_ego
        # This reference will not close the context
        ref_context = egads.c_to_py(self.context.py_to_c())

#==============================================================================
    def test_updateThread(self):
        self.context.updateThread()

#==============================================================================
    def test_makeTransform(self):

        mat1 = ((1,  0,  0,  4),
                (0,  1,  0,  8),
                (0,  0,  1, 12))
        xform = self.context.makeTransform(mat1)
        mat2 = xform.getTransform()

        self.assertEqual(mat1, mat2)

        # Test creating a refernce to an ego
        ref_ego = egads.c_to_py(xform.py_to_c())

#==============================================================================
    def test_makeGeometry(self):

        # construction data
        # vectors are orthogonalized by makeGeometry
        # using orthagonal input vectors means getGeometry will match makeGeometry
        x0     = [0,0,0]
        x1     = [1,0,0]
        xcent  = [0,0,0]
        xax    = [1,0,0]
        yax    = [0,1,0]
        zax    = [0,0,1]
        r      = 1
        majr   = 1.5
        minr   = 1.0
        focus  = 1.0
        vec    = [0, 1, 0]
        offset = 1.0
        angle  = 30.*3.14/180.

        # create the Line (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        self.assertEqual(self.context, line.context)
        data = line.getGeometry()
        self.assertEqual( (egads.CURVE, egads.LINE), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Circle
        rdata = [None]*10
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]
        rdata[2] = xcent[2]
        rdata[3] = xax[0]   # x-axis
        rdata[4] = xax[1]
        rdata[5] = xax[2]
        rdata[6] = yax[0]   # y-axis
        rdata[7] = yax[1]
        rdata[8] = yax[2]
        rdata[9] = r        # radius

        circle = self.context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
        data = circle.getGeometry()
        self.assertEqual( (egads.CURVE, egads.CIRCLE), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Ellipse
        rdata = [None]*11
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]
        rdata[2] = xcent[2]
        rdata[3] = xax[0]   # x-axis
        rdata[4] = xax[1]
        rdata[5] = xax[2]
        rdata[6] = yax[0]   # y-axis
        rdata[7] = yax[1]
        rdata[8] = yax[2]
        rdata[ 9] = majr    # major radius
        rdata[10] = minr    # minor radius

        ellipse = self.context.makeGeometry(egads.CURVE, egads.ELLIPSE, rdata)
        data = ellipse.getGeometry()
        self.assertEqual( (egads.CURVE, egads.ELLIPSE), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Parabola #
        rdata = [None]*10
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]
        rdata[2] = xcent[2]
        rdata[3] = xax[0]   # x-axis
        rdata[4] = xax[1]
        rdata[5] = xax[2]
        rdata[6] = yax[0]   # y-axis
        rdata[7] = yax[1]
        rdata[8] = yax[2]
        rdata[9] = focus    # focus

        parabola = self.context.makeGeometry(egads.CURVE, egads.PARABOLA, rdata)
        data = parabola.getGeometry()
        self.assertEqual( (egads.CURVE, egads.PARABOLA), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Hyperbola
        rdata = [None]*11
        rdata[ 0] = xcent[0] # center
        rdata[ 1] = xcent[1]
        rdata[ 2] = xcent[2]
        rdata[ 3] = xax[0]   # x-axis
        rdata[ 4] = xax[1]
        rdata[ 5] = xax[2]
        rdata[ 6] = yax[0]   # y-axis
        rdata[ 7] = yax[1]
        rdata[ 8] = yax[2]
        rdata[ 9] = majr     # major radius
        rdata[10] = minr     # minor radius

        hyperbola = self.context.makeGeometry(egads.CURVE, egads.HYPERBOLA, rdata)
        data = hyperbola.getGeometry()
        self.assertEqual( (egads.CURVE, egads.HYPERBOLA), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Offset curve
        rdata = [None]*4
        rdata[0] = vec[0] # offset direction
        rdata[1] = vec[1]
        rdata[2] = vec[2]
        rdata[3] = offset # offset magnitude

        curveOffset = self.context.makeGeometry(egads.CURVE, egads.OFFSET, rdata, geom=line)
        data = curveOffset.getGeometry()
        self.assertEqual( (egads.CURVE, egads.OFFSET), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create bezier curve
        npts = 4

        idata = [None]*3
        idata[0] = 0
        idata[1] = npts-1
        idata[2] = npts

        rdata = [(0.00, 0.00, 0.00),
                 (1.00, 0.00, 0.10),
                 (1.50, 1.00, 0.70),
                 (0.25, 0.75, 0.60)]

        curveBezier = self.context.makeGeometry(egads.CURVE, egads.BEZIER, rdata, idata)
        data = curveBezier.getGeometry()
        self.assertEqual( (egads.CURVE, egads.BEZIER), data[0:2])
        self.assertEqual( rdata, data[2] )
        self.assertEqual( idata, data[3] )

        # create Plane
        rdata = [None]*9
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]
        rdata[2] = xcent[2]
        rdata[3] = xax[0]   # x-axis
        rdata[4] = xax[1]
        rdata[5] = xax[2]
        rdata[6] = yax[0]   # y-axis
        rdata[7] = yax[1]
        rdata[8] = yax[2]

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)
        data = plane.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.PLANE), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Spherical
        rdata = [None]*10
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]
        rdata[2] = xcent[2]
        rdata[3] = xax[0]   # x-axis
        rdata[4] = xax[1]
        rdata[5] = xax[2]
        rdata[6] = yax[0]   # y-axis
        rdata[7] = yax[1]
        rdata[8] = yax[2]
        rdata[9] = r        # radius

        spherical = self.context.makeGeometry(egads.SURFACE, egads.SPHERICAL, rdata)
        data = spherical.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.SPHERICAL), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Conical
        rdata = [None]*14
        rdata[ 0] = xcent[0] # center
        rdata[ 1] = xcent[1]
        rdata[ 2] = xcent[2]
        rdata[ 3] = xax[0]   # x-axis
        rdata[ 4] = xax[1]
        rdata[ 5] = xax[2]
        rdata[ 6] = yax[0]   # y-axis
        rdata[ 7] = yax[1]
        rdata[ 8] = yax[2]
        rdata[ 9] = zax[0]   # z-axis
        rdata[10] = zax[1]
        rdata[11] = zax[2]
        rdata[12] = angle    # angle
        rdata[13] = r        # radius

        conical = self.context.makeGeometry(egads.SURFACE, egads.CONICAL, rdata)
        data = conical.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.CONICAL), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create Cylinder
        rdata = [None]*13
        rdata[ 0] = xcent[0] # center
        rdata[ 1] = xcent[1]
        rdata[ 2] = xcent[2]
        rdata[ 3] = xax[0]   # x-axis
        rdata[ 4] = xax[1]
        rdata[ 5] = xax[2]
        rdata[ 6] = yax[0]   # y-axis
        rdata[ 7] = yax[1]
        rdata[ 8] = yax[2]
        rdata[ 9] = zax[0]   # z-axis
        rdata[10] = zax[1]
        rdata[11] = zax[2]
        rdata[12] = r        # radius

        cylindrical = self.context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, rdata)
        data = cylindrical.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.CYLINDRICAL), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create Toroidal
        rdata = [None]*14
        rdata[ 0] = xcent[0] # center
        rdata[ 1] = xcent[1]
        rdata[ 2] = xcent[2]
        rdata[ 3] = xax[0]   # x-axis
        rdata[ 4] = xax[1]
        rdata[ 5] = xax[2]
        rdata[ 6] = yax[0]   # y-axis
        rdata[ 7] = yax[1]
        rdata[ 8] = yax[2]
        rdata[ 9] = zax[0]   # z-axis
        rdata[10] = zax[1]
        rdata[11] = zax[2]
        rdata[12] = majr     # major radius
        rdata[13] = minr     # minor radius

        toroidal = self.context.makeGeometry(egads.SURFACE, egads.TOROIDAL, rdata)
        data = toroidal.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.TOROIDAL), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Revolution surface
        rdata = [None]*6
        rdata[0] = xcent[0] # center
        rdata[1] = xcent[1]-1
        rdata[2] = xcent[2]
        rdata[3] = xax[0]    # direction
        rdata[4] = xax[1]
        rdata[5] = xax[2]

        revolution = self.context.makeGeometry(egads.SURFACE, egads.REVOLUTION, rdata, geom=line)
        data = revolution.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.REVOLUTION), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create the Extrusion surface
        rdata = [None]*3
        rdata[0] = zax[0] # direction
        rdata[1] = zax[1]
        rdata[2] = zax[2]

        extrusion = self.context.makeGeometry(egads.SURFACE, egads.EXTRUSION, rdata, geom=circle)
        data = extrusion.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.EXTRUSION), data[0:2])
        self.assertEqual( rdata, data[2] )

        # create Bezier surface
        nCPu = 4
        nCPv = 4

        idata = [None]*5
        idata[0] = 0
        idata[1] = nCPu-1
        idata[2] = nCPu
        idata[3] = nCPv-1
        idata[4] = nCPv

        rdata = [(0.00, 0.00, 0.00),
                 (1.00, 0.00, 0.10),
                 (1.50, 1.00, 0.70),
                 (0.25, 0.75, 0.60),

                 (0.00, 0.00, 1.00),
                 (1.00, 0.00, 1.10),
                 (1.50, 1.00, 1.70),
                 (0.25, 0.75, 1.60),

                 (0.00, 0.00, 2.00),
                 (1.00, 0.00, 2.10),
                 (1.50, 1.00, 2.70),
                 (0.25, 0.75, 2.60),

                 (0.00, 0.00, 3.00),
                 (1.00, 0.00, 3.10),
                 (1.50, 1.00, 3.70),
                 (0.25, 0.75, 3.60)]

        surfBezier = self.context.makeGeometry(egads.SURFACE, egads.BEZIER, rdata, idata)
        data = surfBezier.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.BEZIER), data[0:2])
        self.assertEqual( rdata, data[2] )
        self.assertEqual( idata, data[3] )

        # create the Offset surface
        rdata = [offset]

        surfBezier = self.context.makeGeometry(egads.SURFACE, egads.OFFSET, rdata, geom=plane)
        data = surfBezier.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.OFFSET), data[0:2])
        self.assertEqual( rdata, data[2] )

#==============================================================================
    def _makePlanarLoop(self, x0, x1, x2, x3):

        nodes = [self.context.makeTopology(egads.NODE, reals=x0),
                 self.context.makeTopology(egads.NODE, reals=x1),
                 self.context.makeTopology(egads.NODE, reals=x2),
                 self.context.makeTopology(egads.NODE, reals=x3)]

        lines = [None]*4
        edges = [None]*4

        # create the Line (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2)**0.5]

        lines[0] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[0], children=nodes[0:2], reals=tdata)

        rdata[0] = x1[0]
        rdata[1] = x1[1]
        rdata[2] = x1[2]
        rdata[3] = x2[0] - x1[0]
        rdata[4] = x2[1] - x1[1]
        rdata[5] = x2[2] - x1[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2)**0.5]

        lines[1] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[1] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[1], children=nodes[1:3], reals=tdata)

        rdata[0] = x2[0]
        rdata[1] = x2[1]
        rdata[2] = x2[2]
        rdata[3] = x3[0] - x2[0]
        rdata[4] = x3[1] - x2[1]
        rdata[5] = x3[2] - x2[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2)**0.5]

        lines[2] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[2] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[2], children=nodes[2:4], reals=tdata)

        rdata[0] = x3[0]
        rdata[1] = x3[1]
        rdata[2] = x3[2]
        rdata[3] = x0[0] - x3[0]
        rdata[4] = x0[1] - x3[1]
        rdata[5] = x0[2] - x3[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2)**0.5]

        lines[3] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[3] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[3], children=[nodes[3], nodes[0]], reals=tdata)

        return self.context.makeTopology(egads.LOOP, egads.CLOSED, children=edges, senses=[egads.SFORWARD]*4)

#==============================================================================
    def test_makeTopology(self):

        # create Bezier surface
        nCPu = 4
        nCPv = 4

        idata = [None]*5
        idata[0] = 0
        idata[1] = nCPu-1
        idata[2] = nCPu
        idata[3] = nCPv-1
        idata[4] = nCPv

        rdata = [(0.00, 0.00, 0.00),
                 (1.00, 0.00, 0.10),
                 (1.50, 1.00, 0.70),
                 (0.25, 0.75, 0.60),

                 (0.00, 0.00, 1.00),
                 (1.00, 0.00, 1.10),
                 (1.50, 1.00, 1.70),
                 (0.25, 0.75, 1.60),

                 (0.00, 0.00, 2.00),
                 (1.00, 0.00, 2.10),
                 (1.50, 1.00, 2.70),
                 (0.25, 0.75, 2.60),

                 (0.00, 0.00, 3.00),
                 (1.00, 0.00, 3.10),
                 (1.50, 1.00, 3.70),
                 (0.25, 0.75, 3.60)]

        surface = self.context.makeGeometry(egads.SURFACE, egads.BEZIER, rdata, idata)

        # make Nodes
        nodes = [None]*4
        nodes[0] = self.context.makeTopology(egads.NODE, reals=[0.00, 0.00, 0.00])
        nodes[1] = self.context.makeTopology(egads.NODE, reals=[0.25, 0.75, 0.60])
        nodes[2] = self.context.makeTopology(egads.NODE, reals=[0.25, 0.75, 3.60])
        nodes[3] = self.context.makeTopology(egads.NODE, reals=[0.00, 0.00, 3.00])

        self.context.makeTopology(*nodes[0].getTopology())
        self.context.makeTopology(*nodes[1].getTopology())
        self.context.makeTopology(*nodes[2].getTopology())
        self.context.makeTopology(*nodes[3].getTopology())


        npts = 4
        idata = [None]*3
        idata[0] = 0
        idata[1] = npts-1
        idata[2] = npts

        curves = [None]*4

        rdata = [(0.00, 0.00, 0.00),
                 (1.00, 0.00, 0.10),
                 (1.50, 1.00, 0.70),
                 (0.25, 0.75, 0.60)]

        curves[0] = self.context.makeGeometry(egads.CURVE, egads.BEZIER, rdata, idata)

        rdata = [(0.25, 0.75, 0.60),
                 (0.25, 0.75, 1.60),
                 (0.25, 0.75, 2.60),
                 (0.25, 0.75, 3.60)]

        curves[1] = self.context.makeGeometry(egads.CURVE, egads.BEZIER, rdata, idata)

        rdata = [(0.25, 0.75, 3.60),
                 (1.50, 1.00, 3.70),
                 (1.00, 0.00, 3.10),
                 (0.00, 0.00, 3.00)]

        curves[2] = self.context.makeGeometry(egads.CURVE, egads.BEZIER, rdata, idata)


        rdata = [(0.00, 0.00, 3.00),
                 (0.00, 0.00, 2.00),
                 (0.00, 0.00, 1.00),
                 (0.00, 0.00, 0.00)]

        curves[3] = self.context.makeGeometry(egads.CURVE, egads.BEZIER, rdata, idata)


        tdata = [0, 1]
        edges = [None]*4
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=curves[0], children=nodes[0:2], reals=tdata)
        edges[1] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=curves[1], children=nodes[1:3], reals=tdata)
        edges[2] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=curves[2], children=nodes[2:4], reals=tdata)
        edges[3] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=curves[3], children=[nodes[3],nodes[0]], reals=tdata)

        self.context.makeTopology(*edges[0].getTopology())
        self.context.makeTopology(*edges[1].getTopology())
        self.context.makeTopology(*edges[2].getTopology())
        self.context.makeTopology(*edges[3].getTopology())

        pcurves = [surface.otherCurve(edge) for edge in edges]

        edges.extend(pcurves)

        # making a face with the loop forward
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, geom=surface, children=edges, senses=[egads.SFORWARD]*4)
        self.context.makeTopology(*loop.getTopology())
        
        self.assertGreater(loop.getArea(), 0)
        self.assertLess(loop.flipObject().getArea(), 0)  # Area is negative as loop is reversed relative to the surface
        
        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=surface, children=[loop])
        self.context.makeTopology(*face.getTopology())
        self.assertGreater(face.getArea(), 0)
        
        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=surface, children=[loop])
        self.context.makeTopology(*face.getTopology())
        self.assertGreater(face.getArea(), 0)

        # making a face with the loop reversed
        edges[2], edges[0] = edges[0], edges[2]
        edges[6], edges[4] = edges[4], edges[6]
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, geom=surface, children=edges, senses=[egads.SREVERSE]*4)

        self.assertLess(loop.getArea(), 0) # Area is negative as loop is reversed relative to the surface
        self.assertGreater(loop.flipObject().getArea(), 0)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=surface, children=[loop])
        self.context.makeTopology(*face.getTopology())
        self.assertGreater(face.getArea(), 0)

        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=surface, children=[loop])
        self.context.makeTopology(*face.getTopology())
        self.assertGreater(face.getArea(), 0)

        # The orientation of the Loops does not matter when making a Face
        # A single SOUTER loop must be marked with the senses input, all other loops are SINNER

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)

        # construct forward faces
        # right handed loop
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3)
        
        # left handed loop
        x0 = [0.25,0.25,0]
        x1 = [0.25,0.75,0]
        x2 = [0.75,0.75,0]
        x3 = [0.75,0.25,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3)
        
        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face0.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face0.egads", True)

        # right handed loop
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3)

        # right handed loop flipped
        x0 = [0.25,0.25,0]
        x1 = [0.75,0.25,0]
        x2 = [0.75,0.75,0]
        x3 = [0.25,0.75,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3).flipObject()

        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face1.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face1.egads", True)

        # left handed loop flipped
        x0 = [0,0,0]
        x1 = [0,1,0]
        x2 = [1,1,0]
        x3 = [1,0,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3).flipObject()

        # left handed loop
        x0 = [0.25,0.25,0]
        x1 = [0.25,0.75,0]
        x2 = [0.75,0.75,0]
        x3 = [0.75,0.25,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3)

        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face2.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face2.egads", True)


        # construct reversed faces
        # left handed loop
        x0 = [0,0,0]
        x1 = [0,1,0]
        x2 = [1,1,0]
        x3 = [1,0,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3)
        
        # right handed loop
        x0 = [0.25,0.25,0]
        x1 = [0.75,0.25,0]
        x2 = [0.75,0.75,0]
        x3 = [0.25,0.75,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face3.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face3.egads", True)

        # left handed loop
        x0 = [0,0,0]
        x1 = [0,1,0]
        x2 = [1,1,0]
        x3 = [1,0,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3)

        # left handed loop flipped
        x0 = [0.25,0.25,0]
        x1 = [0.25,0.75,0]
        x2 = [0.75,0.75,0]
        x3 = [0.75,0.25,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3).flipObject()

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face4.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face4.egads", True)

        # right handed loop flipped
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]
        loop_outer = self._makePlanarLoop(x0, x1, x2, x3).flipObject()

        # right handed loop
        x0 = [0.25,0.25,0]
        x1 = [0.75,0.25,0]
        x2 = [0.75,0.75,0]
        x3 = [0.25,0.75,0]
        loop_inner = self._makePlanarLoop(x0, x1, x2, x3)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face5.egads", True)

        face = self.context.makeTopology(egads.FACE, egads.SREVERSE, geom=plane, children=[loop_outer, loop_inner], senses=[egads.SOUTER, egads.SINNER])
        self.context.makeTopology(*face.getTopology())
        #body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
        #body.saveModel("face5.egads", True)


#==============================================================================
    def test_makeTopology_periodic(self):

        # radius
        r = 1
        xcent = [0, 0, 0]

        # create the Cylindrical
        rdata = [None]*13
        rdata[ 0] = xcent[0] # center
        rdata[ 1] = xcent[1]
        rdata[ 2] = xcent[2]
        rdata[ 3] = 1   # x-axis
        rdata[ 4] = 0
        rdata[ 5] = 0
        rdata[ 6] = 0   # y-axis
        rdata[ 7] = 1
        rdata[ 8] = 0
        rdata[ 9] = 0   # z-axis
        rdata[10] = 0
        rdata[11] = 1
        rdata[12] = r   # radius

        cylindrical = self.context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, rdata)

        (oclass, mtype, rvec, ivec, eref) = cylindrical.getGeometry()
        
        # get the axes for the cylinder
        dx = rvec[3:6]
        dy = rvec[6:9]
        dz = rvec[9:12]
        
        circle = [None]*2

        # create the Circle curve for the base
        rdata = [None]*10
        rdata[0] = rvec[0]; # center
        rdata[1] = rvec[1];
        rdata[2] = rvec[2];
        rdata[3] = dx[0];   # x-axis
        rdata[4] = dx[1];
        rdata[5] = dx[2];
        rdata[6] = dy[0];   # y-axis
        rdata[7] = dy[1];
        rdata[8] = dy[2];
        rdata[9] = r;       # radius
        circle[0] = self.context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)

        # create the Circle curve for the top 
        rdata[0] = rvec[0] + dz[0]*r; # center
        rdata[1] = rvec[1] + dz[1]*r;
        rdata[2] = rvec[2] + dz[2]*r;
        rdata[3] = dx[0];   # x-axis
        rdata[4] = dx[1];
        rdata[5] = dx[2];
        rdata[6] = dy[0];   # y-axis
        rdata[7] = dy[1];
        rdata[8] = dy[2];
        rdata[9] = r;       # radius
        circle[1] = self.context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)

        nodes = [None]*2
 
        # create the Node on the Base
        x0 = [None]*3
        x0[0] = xcent[0] + dx[0]*r;
        x0[1] = xcent[1] + dx[1]*r;
        x0[2] = xcent[2] + dx[2]*r;
        nodes[0] = self.context.makeTopology(egads.NODE, reals=x0)
        
        # create the Node on the Top
        x1 = [None]*3
        x1[0] = xcent[0] + dx[0]*r + dz[0]*r;
        x1[1] = xcent[1] + dx[1]*r + dz[1]*r;
        x1[2] = xcent[2] + dx[2]*r + dz[2]*r;
        nodes[1] = self.context.makeTopology(egads.NODE, reals=x1)

        # create Line and Edge connecting the nodes
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2 + rdata[5]**2)**0.5]

        edges = [None]*8

        # create Line connecting Base and Top
        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=line, children=nodes, reals=tdata)

        TWOPI = 2*math.pi

        # create Edge on the Base
        tdata = [0, TWOPI]
        edges[1] = self.context.makeTopology(egads.EDGE, egads.ONENODE, geom=circle[0], children=[nodes[0]], reals=tdata)
        
        edges[2] = edges[0]; # repeat the line edge */

        # make Edge on the Top
        edges[3] = self.context.makeTopology(egads.EDGE, egads.ONENODE, geom=circle[1], children=[nodes[1]], reals=tdata)

        # create P-curves 
        rdata = [None]*4
        rdata[0] = 0.;    rdata[1] = 0.;    # u == 0 UMIN
        rdata[2] = 0.;    rdata[3] = 1.;
        edges[4+0] = self.context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
        
        rdata[0] = 0.;    rdata[1] = 0;     # v == 0 VMIN
        rdata[2] = 1.;    rdata[3] = 0.;
        edges[4+1] = self.context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
        
        rdata[0] = TWOPI; rdata[1] = 0.;    # u == TWOPI UMAX
        rdata[2] = 0.;    rdata[3] = 1.;
        edges[4+2] = self.context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
        
        rdata[0] = 0.;    rdata[1] = r;      # v == r VMAX
        rdata[2] = 1.;    rdata[3] = 0.;
        edges[4+3] = self.context.makeGeometry(egads.PCURVE, egads.LINE, rdata)

        # create the loop
        senses = [egads.SREVERSE, egads.SFORWARD, egads.SFORWARD, egads.SREVERSE]
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, geom=cylindrical, children=edges, senses=senses)

        # create the face/FaceBody with alternating mtype
        for mtype in [egads.SFORWARD, egads.SREVERSE]:

            face = self.context.makeTopology(egads.FACE, mtype, geom=cylindrical, children=[loop])

            body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])

        # reverse the Line Edge and make again

        # create Line connecting Top and Base
        # This cannot be done with flipObject
        rdata = [None]*6
        rdata[0] = x1[0]
        rdata[1] = x1[1]
        rdata[2] = x1[2]
        rdata[3] = x0[0] - x1[0]
        rdata[4] = x0[1] - x1[1]
        rdata[5] = x0[2] - x1[2]
        tdata = [0, (rdata[3]**2 + rdata[4]**2 + rdata[5]**2)**0.5]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=line, children=[nodes[1], nodes[0]], reals=tdata)

        edges[2]   = edges[0]

        # the PCurves can be reversed with flipObject
        edges[4+0] = edges[4+0].flipObject()
        edges[4+2] = edges[4+2].flipObject()

        # create the loop with the Line Edge reversed (senses[0]
        senses[0] = -senses[0]
        senses[2] = -senses[2]
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, geom=cylindrical, children=edges, senses=senses)

        # create the face/FaceBody with alternating mtype
        for mtype in [egads.SFORWARD, egads.SREVERSE]:

            face = self.context.makeTopology(egads.FACE, mtype, geom=cylindrical, children=[loop])

            body = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])


#==============================================================================
    def test_approximate(self):

        npts = 10
        pts = [None]*npts
        for i in range(npts):
            pts[i] = (math.cos(math.pi*i/npts), math.sin(math.pi*i/npts), 0.0)

        bspline = self.context.approximate([npts, 0], pts)
        data = bspline.getGeometry()
        self.assertEqual( (egads.CURVE, egads.BSPLINE), data[0:2])

        bsplineCopy = self.context.makeGeometry(*data)
        self.assertTrue( bspline.isSame(bsplineCopy) )


        nCPu = 4;
        nCPv = 4;
        pts = [(0.00, 0.00, 0.00),
               (1.00, 0.00, 0.10),
               (1.50, 1.00, 0.70),
               (0.25, 0.75, 0.60),

               (0.00, 0.00, 1.00),
               (1.00, 0.00, 1.10),
               (1.50, 1.00, 1.70),
               (0.25, 0.75, 1.60),

               (0.00, 0.00, 2.00),
               (1.00, 0.00, 2.10),
               (1.50, 1.00, 2.70),
               (0.25, 0.75, 2.60),

               (0.00, 0.00, 3.00),
               (1.00, 0.00, 3.10),
               (1.50, 1.00, 3.70),
               (0.25, 0.75, 3.60)]

        bspline = self.context.approximate([nCPu, nCPv], pts)
        data = bspline.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.BSPLINE), data[0:2])

        bsplineCopy = self.context.makeGeometry(*data)
        self.assertTrue( bspline.isSame(bsplineCopy) )

#==============================================================================
    def test_skinning(self):

        curves = [None]*2

        npts = 10
        pts = [None]*npts
        for i in range(npts):
            pts[i] = (math.cos(math.pi*i/npts), math.sin(math.pi*i/npts), 0.0)

        curves[0] = self.context.approximate([npts, 0], pts)

        pts = [None]*npts
        for i in range(npts):
            pts[i] = (math.cos(math.pi*i/npts), math.sin(math.pi*i/npts), 1.0)

        curves[1] = self.context.approximate([npts, 0], pts)

        bspline = egads.skinning(curves, degree=1)

        data = bspline.getGeometry()
        self.assertEqual( (egads.SURFACE, egads.BSPLINE), data[0:2])

#==============================================================================
    def test_attribute(self):

        # construction data
        x0     = [0,0,0]
        x1     = [1,0,0]

        # create the Line (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)

        # ints
        attrInt = [1, 2]
        line.attributeAdd("name", attrInt)
        data = line.attributeRet("name")
        self.assertEqual( attrInt, data )

        self.assertEqual(1, line.attributeNum())

        name, data = line.attributeGet(1)
        self.assertEqual( "name", name )
        self.assertEqual( attrInt, data )

        # float
        attrFloat= [1.0, 2.0]
        line.attributeAdd("name", attrFloat)
        data = line.attributeRet("name")
        self.assertEqual( attrFloat, data )

        name, data = line.attributeGet(1)
        self.assertEqual( "name", name )
        self.assertEqual( attrFloat, data )

        # string
        attrString = "attrString"
        line.attributeAdd("name", attrString)
        data = line.attributeRet("name")
        self.assertEqual( attrString, data )

        name, data = line.attributeGet(1)
        self.assertEqual( "name", name )
        self.assertEqual( attrString, data )

        # csystem
        attrCsys = egads.csystem([0,1, 0,1,0])
        line.attributeAdd("name", attrCsys)

        data = line.attributeRet("name")
        self.assertEqual( attrCsys, data )
        self.assertIsInstance( data, egads.csystem )

        name, data = line.attributeGet(1)
        self.assertEqual( "name", name )
        self.assertEqual( attrCsys, data )
        self.assertIsInstance( data, egads.csystem )

        # duplicate attributes
        line2 = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        line.attributeDup(line2)
        self.assertEqual(1, line2.attributeNum())

        # Add a second attribute
        line.attributeAdd("name2", attrInt)
        self.assertEqual(2, line.attributeNum())

        # Delete it
        line.attributeDel("name2")
        self.assertEqual(1, line.attributeNum())

        # Delete all
        line.attributeDel()
        self.assertEqual(0, line.attributeNum())


#==============================================================================
    def test_fullAttribute(self):

        self.context.setFullAttrs(1)

        # construction data
        x0     = [0,0,0]
        x1     = [1,0,0]

        # create the Line (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)


        # ints
        attrInt1 = [1, 2]
        line.attributeAddSeq("name", attrInt1)
        attrInt2 = [3, 4]
        line.attributeAddSeq("name", attrInt2)

        data = line.attributeRetSeq("name", 1)
        self.assertEqual( attrInt1, data )
        data = line.attributeRetSeq("name", 2)
        self.assertEqual( attrInt2, data )

        self.assertEqual(2, line.attributeNum())
        self.assertEqual(2, line.attributeNumSeq("name"))

        name, data = line.attributeGet(1)
        self.assertEqual( "name 1", name )
        self.assertEqual( attrInt1, data )

        name, data = line.attributeGet(2)
        self.assertEqual( "name 2", name )
        self.assertEqual( attrInt2, data )

        line.attributeDel("name")
        self.assertEqual(0, line.attributeNum())


        # float
        attrFloat1 = [1.0, 2.0]
        line.attributeAddSeq("name", attrFloat1)
        attrFloat2 = [3.0, 4.0]
        line.attributeAddSeq("name", attrFloat2)

        data = line.attributeRetSeq("name", 1)
        self.assertEqual( attrFloat1, data )

        data = line.attributeRetSeq("name", 2)
        self.assertEqual( attrFloat2, data )

        name, data = line.attributeGet(1)
        self.assertEqual( "name 1", name )
        self.assertEqual( attrFloat1, data )

        name, data = line.attributeGet(2)
        self.assertEqual( "name 2", name )
        self.assertEqual( attrFloat2, data )

        line.attributeDel("name")
        self.assertEqual(0, line.attributeNum())


        # string
        attrString1 = "attrString1"
        line.attributeAddSeq("name", attrString1)

        attrString2 = "attrString2"
        line.attributeAddSeq("name", attrString2)

        data = line.attributeRetSeq("name", 1)
        self.assertEqual( attrString1, data )

        data = line.attributeRetSeq("name", 2)
        self.assertEqual( attrString2, data )

        name, data = line.attributeGet(1)
        self.assertEqual( "name 1", name )
        self.assertEqual( attrString1, data )

        name, data = line.attributeGet(2)
        self.assertEqual( "name 2", name )
        self.assertEqual( attrString2, data )

        line.attributeDel("name")
        self.assertEqual(0, line.attributeNum())


        # csystem
        attrCsys1 = egads.csystem([0,1, 0,1,0])
        line.attributeAddSeq("name", attrCsys1)

        attrCsys2 = egads.csystem([0,1, 0,1,1])
        line.attributeAddSeq("name", attrCsys2)

        data = line.attributeRetSeq("name", 1)
        self.assertEqual( attrCsys1, data )
        self.assertIsInstance( data, egads.csystem )

        data = line.attributeRetSeq("name", 2)
        self.assertEqual( attrCsys2, data )
        self.assertIsInstance( data, egads.csystem )

        name, data = line.attributeGet(1)
        self.assertEqual( "name 1", name )
        self.assertEqual( attrCsys1, data )
        self.assertIsInstance( data, egads.csystem )

        name, data = line.attributeGet(2)
        self.assertEqual( "name 2", name )
        self.assertEqual( attrCsys2, data )
        self.assertIsInstance( data, egads.csystem )


        # duplicate attributes
        line2 = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        line.attributeDup(line2)
        self.assertEqual(2, line2.attributeNum())

        # Add a second attribute
        line.attributeAdd("name2", attrInt1)
        self.assertEqual(3, line.attributeNum())

        # Delete it
        line.attributeDel("name2")
        self.assertEqual(2, line.attributeNum())

        # Delete all
        line.attributeDel()
        self.assertEqual(0, line.attributeNum())

#==============================================================================
    def test_topology(self):

        # construction data
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]

        nodes = [None]*4
        nodes[0] = self.context.makeTopology(egads.NODE, reals=x0)
        nodes[1] = self.context.makeTopology(egads.NODE, reals=x1)
        nodes[2] = self.context.makeTopology(egads.NODE, reals=x2)
        nodes[3] = self.context.makeTopology(egads.NODE, reals=x3)
        for i in range(4): nodes[i].name = "node" + str(i)

        (oclass, mtype, geom, lim, children, sens) = nodes[0].getTopology()

        self.assertEqual(geom, None)
        self.assertEqual((egads.NODE, 0), (oclass, mtype))
        self.assertEqual(lim, x0)
        self.assertEqual(children, [])
        self.assertEqual(sens, None)

        lines = [None]*4
        edges = [None]*4
        tdata = [0, 1]

        # create the Line (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        lines[0] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[0], children=nodes[0:2], reals=tdata)

        rdata[0] = x1[0]
        rdata[1] = x1[1]
        rdata[2] = x1[2]
        rdata[3] = x2[0] - x1[0]
        rdata[4] = x2[1] - x1[1]
        rdata[5] = x2[2] - x1[2]

        lines[1] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[1] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[1], children=nodes[1:3], reals=tdata)

        rdata[0] = x2[0]
        rdata[1] = x2[1]
        rdata[2] = x2[2]
        rdata[3] = x3[0] - x2[0]
        rdata[4] = x3[1] - x2[1]
        rdata[5] = x3[2] - x2[2]

        lines[2] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[2] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[2], children=nodes[2:4], reals=tdata)

        rdata[0] = x3[0]
        rdata[1] = x3[1]
        rdata[2] = x3[2]
        rdata[3] = x0[0] - x3[0]
        rdata[4] = x0[1] - x3[1]
        rdata[5] = x0[2] - x3[2]

        lines[3] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[3] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[3], children=[nodes[3], nodes[0]], reals=tdata)

        for i in range(4): lines[i].name = "line" + str(i)
        for i in range(4): edges[i].name = "edge" + str(i)

        (oclass, mtype, geom, lim, children, sens) = edges[0].getTopology()

        self.assertTrue(lines[0].isSame(geom))
        self.assertTrue(lines[0].isEquivalent(geom))
        self.assertEqual((egads.EDGE, egads.TWONODE), (oclass, mtype))
        self.assertEqual(lim, tdata)
        self.assertTrue(nodes[0].isSame(children[0]))
        self.assertTrue(nodes[1].isSame(children[1]))
        self.assertEqual(sens, None)

        sens = [1,1,1,1]
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, children=edges, senses=sens)
        loop.name = "loop"

        (oclass, mtype, geom, lim, children, sens) = loop.getTopology()

        self.assertEqual(geom, None)
        self.assertEqual((egads.LOOP, egads.CLOSED), (oclass, mtype))
        self.assertEqual(lim, [])
        self.assertTrue(children[0].isEquivalent(edges[0]))
        self.assertEqual(sens, [1,1,1,1])

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)
        plane.name = "plane"

        face = self.context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop])
        face.name = "face"

        (oclass, mtype, geom, lim, children, sens) = face.getTopology()

        self.assertTrue(plane.isSame(geom))
        self.assertEqual((egads.FACE, egads.SFORWARD), (oclass, mtype))
        self.assertEqual(lim, [0.,1.,0.,1.])
        self.assertTrue(children[0].isSame(loop))
        self.assertEqual(sens, [1])

        faceBody = self.context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])

        bodyNodes = faceBody.getBodyTopos(egads.NODE)
        bodyEdges = faceBody.getBodyTopos(egads.EDGE)
        bodyFaces = faceBody.getBodyTopos(egads.FACE)

        self.assertTrue( all([nodes[i].isEquivalent(bodyNodes[i]) for i in range(4)]) )
        self.assertTrue( all([edges[i].isEquivalent(bodyEdges[i]) for i in range(4)]) )
        self.assertTrue( face.isEquivalent(bodyFaces[0]) )

        body = bodyFaces[0].getBody()
        self.assertTrue( body.isEquivalent(faceBody) )
        self.assertEqual( body, faceBody )

        for i in range(4):
            self.assertEqual( faceBody.indexBodyTopo(bodyNodes[i]), i+1 )
            self.assertEqual( faceBody.indexBodyTopo(bodyEdges[i]), i+1 )

        for i in range(4):
            self.assertEqual( faceBody.objectBodyTopo(egads.NODE, i+1), bodyNodes[i] )
            self.assertEqual( faceBody.objectBodyTopo(egads.EDGE, i+1), bodyEdges[i] )

        # Make a model
        # Modules take ownership if bodies, so a body placed in a model can no longer be deleted
        model = self.context.makeTopology(egads.MODEL, children=[faceBody])

        # exlicitly test deleting in reverse order to test referencing
        faceBody.name = "faceBody"
        model.name = "model"
        del model
        del faceBody
        del nodes
        del lines
        del edges
        del loop
        del plane
        del face

#==============================================================================
    def _makeFaceBody(self, context):

        # construction data
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]

        # make Nodes
        nodes = [None]*4
        nodes[0] = context.makeTopology(egads.NODE, reals=x0)
        nodes[1] = context.makeTopology(egads.NODE, reals=x1)
        nodes[2] = context.makeTopology(egads.NODE, reals=x2)
        nodes[3] = context.makeTopology(egads.NODE, reals=x3)

        # make lines and edges
        lines = [None]*4
        edges = [None]*4
        tdata = [0, 1]

        # Line data (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        lines[0] = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[0], children=nodes[0:2], reals=tdata)

        rdata[0] = x1[0]
        rdata[1] = x1[1]
        rdata[2] = x1[2]
        rdata[3] = x2[0] - x1[0]
        rdata[4] = x2[1] - x1[1]
        rdata[5] = x2[2] - x1[2]

        lines[1] = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[1] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[1], children=nodes[1:3], reals=tdata)

        rdata[0] = x2[0]
        rdata[1] = x2[1]
        rdata[2] = x2[2]
        rdata[3] = x3[0] - x2[0]
        rdata[4] = x3[1] - x2[1]
        rdata[5] = x3[2] - x2[2]

        lines[2] = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[2] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[2], children=nodes[2:4], reals=tdata)

        rdata[0] = x3[0]
        rdata[1] = x3[1]
        rdata[2] = x3[2]
        rdata[3] = x0[0] - x3[0]
        rdata[4] = x0[1] - x3[1]
        rdata[5] = x0[2] - x3[2]

        lines[3] = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[3] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[3], children=[nodes[3], nodes[0]], reals=tdata)

        # make the loop
        sens = [egads.SFORWARD]*4
        loop = context.makeTopology(egads.LOOP, egads.CLOSED, children=edges, senses=sens)

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)
        face = context.makeTopology(egads.FACE, egads.SFORWARD, geom=plane, children=[loop])

        return context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])

#==============================================================================
    def test_tolerance(self):

        # construction data
        x0 = [0,0,0]
        x1 = [1,0,0]

        # make Nodes
        nodes = [None]*2
        nodes[0] = self.context.makeTopology(egads.NODE, reals=x0)
        nodes[1] = self.context.makeTopology(egads.NODE, reals=x1)

        # Line data (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]
        tdata = [0, 1]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edge = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=line, children=nodes, reals=tdata)

        # tolerance
        tol = edge.tolerance()
        tol = edge.getTolerance()
        #edge.setTolerance(2*tol)
        #self.assertAlmostEqual(2*tol, edge.getTolerance(), 8)

#==============================================================================
    def test_getInfo(self):

        faceBody = self._makeFaceBody(self.context)

        # getInfo
        (oclass, mtype, topObj, prev, next) = faceBody.getInfo()
        self.assertEqual((egads.BODY, egads.FACEBODY), (oclass, mtype))

#==============================================================================
    def test_flip(self):

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)

        flipped = plane.flipObject()

        # getInfo
        (oclass, mtype, topObj, prev, next) = flipped.getInfo()
        self.assertEqual((egads.SURFACE, egads.PLANE), (oclass, mtype))

#==============================================================================
    def test_properties(self):

        faceBody = self._makeFaceBody(self.context)

        nodes = faceBody.getBodyTopos(egads.NODE)
        edges = faceBody.getBodyTopos(egads.EDGE)
        faces = faceBody.getBodyTopos(egads.FACE)

        (oclass, mtype, line, lim, children, sens) = edges[0].getTopology()

        (trange, periodic) = edges[0].getRange()

        # length
        self.assertEqual(line.arcLength(trange[0],trange[1]), 1.0)

        (oclass, mtype, plane, lim, children, sens) = faces[0].getTopology()

        (urange, periodic) = faces[0].getRange()

        # area
        self.assertEqual(faces[0].getArea(), 1.0)
        self.assertEqual(plane.getArea(urange), 1.0)

        # bounding box
        bbox = faceBody.getBoundingBox()
        self.assertEqual(bbox, [(0,0,0), (1,1,0)])

        # mass properties
        volume, area, CG, I = faceBody.getMassProperties()
        self.assertEqual(volume, 0)
        self.assertEqual(area, 1.0)
        self.assertEqual(CG, (0.5, 0.5, 0.0))
        Itrue = (1/12.,    0.,   0.,
                    0., 1/12.,   0.,
                    0.,    0., 1/6.)
        for i in range(9):
            self.assertAlmostEqual(I[i], Itrue[i], 9)

        # curvature

        # create the Circle
        rdata = [None]*10
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1   # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0   # y-axis
        rdata[7] = 1
        rdata[8] = 0
        rdata[9] = 1        # radius

        circle = self.context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
        data = circle.curvature(0)
        self.assertEqual((1.0, (0.0, 1.0, 0.0)), data)

        # create Cylinder
        rdata = [None]*13
        rdata[ 0] = 0 # center
        rdata[ 1] = 0
        rdata[ 2] = 0
        rdata[ 3] = 1   # x-axis
        rdata[ 4] = 0
        rdata[ 5] = 0
        rdata[ 6] = 0   # y-axis
        rdata[ 7] = 1
        rdata[ 8] = 0
        rdata[ 9] = 0   # z-axis
        rdata[10] = 0
        rdata[11] = 1
        rdata[12] = 1   # radius

        cylindrical = self.context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, rdata)
        data = cylindrical.curvature([0,0])
        self.assertEqual((0.0, (0.0, 0.0, 1.0), -1.0, (0.0, 1.0, 0.0)), data)

#==============================================================================
    def test_modelIO(self):

        faceBody = self._makeFaceBody(self.context)

        # save/load Model
        faceBody.saveModel("test.egads",True)
        model = self.context.loadModel("test.egads")
        self.assertTrue( os.path.exists("test.egads") )
        os.remove("test.egads")

        (oclass, mtype, geom, lim, children, sens) = model.getTopology()

        self.assertTrue( faceBody.isEquivalent(children[0]) )

        # exlicitly test deleting the body before the model
        del faceBody
        del model

#==============================================================================
    def test_exportModel(self):

        box = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        model = self.context.makeTopology(egads.MODEL, children=[box])

        stream = model.exportModel()
        self.assertIsInstance(stream, bytes)

#==============================================================================
    def test_evaluate(self):

        faceBody = self._makeFaceBody(self.context)

        nodes = faceBody.getBodyTopos(egads.NODE)
        edges = faceBody.getBodyTopos(egads.EDGE)
        faces = faceBody.getBodyTopos(egads.FACE)

        eval = nodes[0].evaluate(None)
        self.assertEqual(eval, (0,0,0))

        eval = edges[0].evaluate(0.5)
        self.assertEqual(eval, ((0.5,0.0,0.),  # value
                                (1.0,0.0,0.),  # derivative
                                (0.0,0.0,0.))) # 2nd-derivative

        eval = faces[0].evaluate([0.5,0.5])
        self.assertEqual(eval, ((0.5,0.5,0),                  # value
                                ((1,0,0), (0,1,0)),           # defivative
                                ((0,0,0), (0,0,0), (0,0,0)))) # 2nd derivative

#==============================================================================
    def test_invEvaluate(self):

        faceBody = self._makeFaceBody(self.context)

        nodes = faceBody.getBodyTopos(egads.NODE)
        edges = faceBody.getBodyTopos(egads.EDGE)
        faces = faceBody.getBodyTopos(egads.FACE)

        inveval = edges[0].invEvaluate([0.5, 0., 0.])
        self.assertEqual(inveval, (0.5,          # parameter
                                  [0.5,0.0,0.])) # coordinate

        inveval = faces[0].invEvaluate([0.5, 0.5, 0.])
        self.assertEqual(inveval, ([0.5, 0.5],    # parameter
                                   [0.5,0.5,0.])) # coordinate

#==============================================================================
    def test_makeSolidBody(self):

        box = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        data = box.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

#==============================================================================
    def test_generalBoolean(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        box1 = self.context.makeSolidBody(egads.BOX, [-0.5,-0.5,-0.5, 1,1,1])

        subtraction = box0.generalBoolean(box1, egads.SUBTRACTION)
        Intersection = box0.generalBoolean(box1, egads.SUBTRACTION)
        fusion = box0.generalBoolean(box1, egads.FUSION)

#==============================================================================
    def test_fuseSheets(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        box1 = self.context.makeSolidBody(egads.BOX, [-0.5,-0.5,-0.5, 1,1,1])

        shell0 = box0.getBodyTopos(egads.SHELL)[0]
        shell1 = box1.getBodyTopos(egads.SHELL)[0]

        sheet0 = self.context.makeTopology(egads.BODY, egads.SHEETBODY, children=[shell0])
        sheet1 = self.context.makeTopology(egads.BODY, egads.SHEETBODY, children=[shell1])

        fuse = sheet0.fuseSheets(sheet1)

#==============================================================================
    def test_intersection(self):

        faceBody = self._makeFaceBody(self.context)

        box0 = self.context.makeSolidBody(egads.BOX, [0.25,0.25,-0.25, 0.5,0.5,0.5])

        pairs, model = box0.intersection(faceBody)

        box1 = box0.imprintBody(pairs)

        faces = box1.getBodyTopos(egads.FACE)
        self.assertEqual(10, len(faces))

#==============================================================================
    def test_makeFace(self):

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1   # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0   # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)

        # make Face from a Plane
        face = plane.makeFace(egads.SFORWARD, [0,1, 0,1])

        data = face.getInfo()
        self.assertEqual((egads.FACE, egads.SFORWARD), data[:2])

        face = plane.makeFace(egads.SREVERSE, [0,1, 0,1])

        data = face.getInfo()
        self.assertEqual((egads.FACE, egads.SREVERSE), data[:2])

        # make Face from a Loop. 
        # This creates a Plane that must be tracked by egads.py

        # make the loop
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]
        loop = self._makePlanarLoop(x0, x1, x2, x3)

        face = loop.makeFace(egads.SFORWARD)
        
        data = face.getInfo()
        self.assertEqual((egads.FACE, egads.SFORWARD), data[:2])


#==============================================================================
    def test_filletBody(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        edges = box0.getBodyTopos(egads.EDGE)
        box1, map = box0.filletBody(edges, 0.1)

        data = box1.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

#==============================================================================
    def test_chamferBody(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        edges = box0.getBodyTopos(egads.EDGE)
        faces = []
        for e in edges:
            faces.append(box0.getBodyTopos(egads.FACE, e)[0])
        box1, map = box0.chamferBody(edges, faces, 0.1, 0.15)

        data = box1.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

#==============================================================================
    def test_hollowBody(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        box1, map = box0.hollowBody(None, 0.1, 0)

        data = box1.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

        faces = [box0.getBodyTopos(egads.FACE)[0]]

        box2, map = box0.hollowBody(faces, 0.1, 0)

        data = box2.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])


#==============================================================================
    def test_extrude(self):

        faceBody = self._makeFaceBody(self.context)

        body = faceBody.extrude(1, [0,0,1])

        faces = body.getBodyTopos(egads.FACE)
        self.assertEqual(6, len(faces))


#==============================================================================
    def test_rotate(self):

        faceBody = self._makeFaceBody(self.context)

        body = faceBody.rotate(180, [(0,-1,0), (1,0,0)])

        faces = body.getBodyTopos(egads.FACE)
        self.assertEqual(6, len(faces))


#==============================================================================
    def test_sweep(self):

        # construction data
        x0 = [0,0,0]
        x1 = [0,0,1]

        # make Nodes
        nodes = [None]*2
        nodes[0] = self.context.makeTopology(egads.NODE, reals=x0)
        nodes[1] = self.context.makeTopology(egads.NODE, reals=x1)

        # Line data (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        tdata = [0,1]

        line = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edge = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=line, children=nodes, reals=tdata)

        faceBody = self._makeFaceBody(self.context)

        body = faceBody.sweep(edge, 0)

        faces = body.getBodyTopos(egads.FACE)
        self.assertEqual(6, len(faces))

#==============================================================================
    def test_getEdgeUV(self):

        faceBody = self._makeFaceBody(self.context)

        faces = faceBody.getBodyTopos(egads.FACE)
        edges = faceBody.getBodyTopos(egads.EDGE)

        uv = faces[0].getEdgeUV(edges[0], 0.5)
        self.assertEqual(uv, (0.5, 0.))

#==============================================================================
    def test_replaceFaces(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        faces = box0.getBodyTopos(egads.FACE)

        # Remove two faces
        replace = [(faces[0], None), (faces[1], None)]

        sheet = box0.replaceFaces(replace)
        faces = sheet.getBodyTopos(egads.FACE)
        self.assertEqual(4, len(faces))

#==============================================================================
    def test_tessBody(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        tess = box0.makeTessBody([0.1, 0.1, 30.])

        xyz, t = tess.getTessEdge(1)
        lIndex = tess.getTessLoops(1)

        body, stat, icode, npts = tess.statusTessBody()

        xyz_global = [None]*npts
        for j in range(npts):
            ptype, pindex, xyz_global[j] = tess.getGlobal(j+1)

        tris_global = []
        nfaces = len(box0.getBodyTopos(egads.FACE))
        for fIndex in range(nfaces):
            xyz, uv, ptype, pindex, tris, tric = tess.getTessFace(fIndex+1)

            for j in range(len(tris)):
                for n in range(3):
                    tris_global.append(tess.localToGlobal(fIndex+1, tris[j][n]))

        # Build up a tesselleation (by copying)
        tess2 = box0.initTessBody()

        nedges = len(box0.getBodyTopos(egads.EDGE))
        for i in range(nedges):
            xyz, t = tess.getTessEdge(i+1)
            tess2.setTessEdge(i+1, xyz, t)

        for i in range(nfaces):
            xyz, uv, ptype, pindex, tris, tric = tess.getTessFace(i+1)
            tess2.setTessFace(i+1, xyz, uv, tris)

        # finish the new tessellation
        tess2.statusTessBody()

        #-------
        # Open up a tessellation and close it again
        tess.openTessBody()
        tess.statusTessBody()

        #-------
        # Open up a tessellation and finish it
        tess.openTessBody()
        tess.finishTess([0.1,0.1,10])

        #-------
        # Quads
        tess.makeQuads([0.1,0.1,10], 1)
        xyz, uv, type, index, npatch = tess.getQuads(1)
        self.assertEqual(npatch, 1)

        fList = tess.getTessQuads()
        self.assertEqual(fList, [1])

        n1, n2, pvindex, pbounds = tess.getPatch(1, 1)

        #-------
        # Modifying Edge vertex
        eIndex = 1
        xyz, t = tess.getTessEdge(eIndex)

        vIndex = int(len(t)/2)

        tess.moveEdgeVert(eIndex, vIndex, t[vIndex-1]*1.01)
        tess.deleteEdgeVert(eIndex, vIndex, -1)
        tess.insertEdgeVerts(eIndex, vIndex, [t[vIndex]+.025])

        #-------
        edges = box0.getBodyTopos(egads.EDGE)
        facedg = box0.getBodyTopos(egads.FACE, edges[0])
        facedg.append( edges[0] )

        # Renamke the tessellation
        tess.remakeTess(facedg, [0.2,0.1,10])

        #-------
        # crate a quad tessellation
        quadTess = tess.quadTess()

        #-------
        # mass properties
        volume, area, CG, I = quadTess.tessMassProps()
        self.assertAlmostEqual(volume, 1.0, 9)
        self.assertAlmostEqual(area, 6.0, 9)
        for i in range(3):
            self.assertAlmostEqual(CG[i], (0.5, 0.5, 0.5)[i], 9)
        Itrue = (1/6.,   0.,   0.,
                   0., 1/6.,   0.,
                   0.,   0., 1/6.)
        for i in range(9):
            self.assertAlmostEqual(I[i], Itrue[i], 9)


#==============================================================================
    def test_tessGeom(self):

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)
        tess = plane.makeTessGeom([0,1, 0,1], [5,5])

        sizes, xyz = tess.getTessGeom()
        self.assertEqual(sizes, (5,5))

#==============================================================================
    def test_fitTriangles(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        tess = box0.makeTessBody([0.05, 0.1, 30.])

        xyz, uv, ptype, pindex, tris, tric = tess.getTessFace(1)

        surf = self.context.fitTriangles(xyz, tris)
        data = surf.getInfo()
        self.assertEqual((egads.SURFACE, egads.BSPLINE), data[:2])

        surf = self.context.fitTriangles(xyz, tris, tric)
        data = surf.getInfo()
        self.assertEqual((egads.SURFACE, egads.BSPLINE), data[:2])

#==============================================================================
    def test_mapBody(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        box1 = self.context.makeSolidBody(egads.BOX, [0.1,0.1,0.1, 1,1,1])

        faces0 = box0.getBodyTopos(egads.FACE)
        faces1 = box1.getBodyTopos(egads.FACE)

        faces0[1].attributeAdd("faceMap", [1])
        faces1[0].attributeAdd("faceMap", [1])

        map = box0.mapBody(box1, "faceMap")

        tess0 = box0.makeTessBody([0.1, 0.1, 15])
        tess1 = tess0.mapTessBody(box1)

        xyz, uv, ptype, pindex, tris, tric = tess0.getTessFace(1)

        ifaces = [1]*int(len(uv)/2)
        result, itris = tess1.locateTessBody(ifaces, uv)

#==============================================================================
    def test_getWindingAngle(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        
        edges = box0.getBodyTopos(egads.EDGE)
        
        trange, per = edges[0].getRange()
        t = (trange[0] + trange[1])/2
        
        angle = edges[0].getWindingAngle(t)
        self.assertAlmostEqual(270.0, angle, 5)

#==============================================================================
    def test_otherCurve(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        edges = box0.getBodyTopos(egads.EDGE)

        faces = box0.getBodyTopos(egads.FACE, edges[0])

        pcurv = faces[0].otherCurve(edges[0])
        data = pcurv.getInfo()
        self.assertEqual((egads.PCURVE, egads.LINE), data[:2])

#==============================================================================
    def test_isoCline(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        faces = box0.getBodyTopos(egads.FACE)
        uvrange, periodic = faces[0].getRange()

        (oclass, mtype, geom, lim, children, sens) = faces[0].getTopology()

        ucurv = geom.isoCline(0, (uvrange[0]+uvrange[1])/2)
        data = ucurv.getInfo()
        self.assertEqual((egads.CURVE, egads.LINE), data[:2])

        vcurv = geom.isoCline(1, (uvrange[2]+uvrange[3])/2)
        data = vcurv.getInfo()
        self.assertEqual((egads.CURVE, egads.LINE), data[:2])

#==============================================================================
    def test_convertToBSpline(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        faces = box0.getBodyTopos(egads.FACE)

        bspline = faces[0].convertToBSpline()
        data = bspline.getInfo()
        self.assertEqual((egads.SURFACE, egads.BSPLINE), data[:2])

        # create Plane
        rdata = [None]*9
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1 # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0 # y-axis
        rdata[7] = 1
        rdata[8] = 0

        plane = self.context.makeGeometry(egads.SURFACE, egads.PLANE, rdata)

        bspline = plane.convertToBSplineRange([0,1, 0,1])
        data = bspline.getInfo()
        self.assertEqual((egads.SURFACE, egads.BSPLINE), data[:2])

#==============================================================================
    def test_inTopology(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        faces = box0.getBodyTopos(egads.FACE)
        uvrange, periodic = faces[0].getRange()
        bbox = faces[0].getBoundingBox()

        self.assertTrue( faces[0].inTopology([(bbox[0][0]+bbox[1][0])/2, (bbox[0][1]+bbox[1][1])/2, (bbox[0][2]+bbox[1][2])/2] ) )

        self.assertTrue( faces[0].inFace([(uvrange[0]+uvrange[1])/2, (uvrange[2]+uvrange[3])/2] ) )

#==============================================================================
    def test_makeLoop(self):

        # construction data
        x0 = [0,0,0]
        x1 = [1,0,0]
        x2 = [1,1,0]
        x3 = [0,1,0]

        # make Nodes
        nodes = [None]*4
        nodes[0] = self.context.makeTopology(egads.NODE, reals=x0)
        nodes[1] = self.context.makeTopology(egads.NODE, reals=x1)
        nodes[2] = self.context.makeTopology(egads.NODE, reals=x2)
        nodes[3] = self.context.makeTopology(egads.NODE, reals=x3)

        # make lines and edges
        lines = [None]*4
        edges = [None]*6 #Add some extra None edges
        tdata = [0, 1]

        # Line data (point and direction)
        rdata = [None]*6
        rdata[0] = x0[0]
        rdata[1] = x0[1]
        rdata[2] = x0[2]
        rdata[3] = x1[0] - x0[0]
        rdata[4] = x1[1] - x0[1]
        rdata[5] = x1[2] - x0[2]

        lines[0] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[0], children=nodes[0:2], reals=tdata)

        rdata[0] = x1[0]
        rdata[1] = x1[1]
        rdata[2] = x1[2]
        rdata[3] = x2[0] - x1[0]
        rdata[4] = x2[1] - x1[1]
        rdata[5] = x2[2] - x1[2]

        lines[1] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[1] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[1], children=nodes[1:3], reals=tdata)

        rdata[0] = x2[0]
        rdata[1] = x2[1]
        rdata[2] = x2[2]
        rdata[3] = x3[0] - x2[0]
        rdata[4] = x3[1] - x2[1]
        rdata[5] = x3[2] - x2[2]

        lines[2] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[2] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[2], children=nodes[2:4], reals=tdata)

        rdata[0] = x3[0]
        rdata[1] = x3[1]
        rdata[2] = x3[2]
        rdata[3] = x0[0] - x3[0]
        rdata[4] = x0[1] - x3[1]
        rdata[5] = x0[2] - x3[2]

        lines[3] = self.context.makeGeometry(egads.CURVE, egads.LINE, rdata)
        edges[3] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=lines[3], children=[nodes[3], nodes[0]], reals=tdata)

        loop, edges = egads.makeLoop(edges)
        self.assertEqual([None]*6, edges)

#==============================================================================
    def test_sewFaces(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])

        faces = box0.getBodyTopos(egads.FACE)
        model = egads.sewFaces(faces)
        data = model.getInfo()
        self.assertEqual(egads.MODEL, data[0])

#==============================================================================
    def test_matchBodyEdges(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        box1 = self.context.makeSolidBody(egads.BOX, [1,0,0, 1,1,1])

        matches = box0.matchBodyEdges(box1)
        self.assertEqual(matches, [(5, 1), (6, 2), (7, 3), (8, 4)])

#==============================================================================
    def test_matchBodyFaces(self):

        box0 = self.context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
        box1 = self.context.makeSolidBody(egads.BOX, [1,0,0, 1,1,1])

        matches = box0.matchBodyFaces(box1)
        self.assertEqual(matches, [(2,1)])

        faces0 = box0.getBodyTopos(egads.FACE)
        faces1 = box1.getBodyTopos(egads.FACE)

        for match in matches:
            faces0[match[0]-1] = None
            faces1[match[1]-1] = None
        faces = faces0 + faces1
        faces = [f for f in faces if f]

        model = egads.sewFaces(faces)
        data = model.getInfo()
        self.assertEqual(egads.MODEL, data[0])

#==============================================================================
    def test_blend(self):

        # sections for blend
        secs = [None]*4

        # first section is a NODE
        secs[0] = self.context.makeTopology(egads.NODE, reals=[0.5,0.5,-1])

        # FaceBody for 2nd section
        secs[1] = self._makeFaceBody(self.context)

        mat1 = ((1,  0,  0,  0), #Translate 1 in z-direction
                (0,  1,  0,  0),
                (0,  0,  1,  1))
        xform = self.context.makeTransform(mat1)

        # Translate 2nd section to make 3rd
        secs[2] = secs[1].copyObject(xform)

        # Finish with a NODE
        secs[3] = self.context.makeTopology(egads.NODE, reals=[0.5,0.5,2])

        blend0 = egads.blend(secs)
        data = blend0.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

        rc = [0.2, 1, 0, 0,  0.4, 0, 1, 0]
        blend2 = egads.blend(secs, rc1=rc)
        data = blend0.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

        blend2 = egads.blend(secs, rc1=rc, rc2=rc)
        data = blend2.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

        blend3 = egads.blend(secs, rc2=rc)
        data = blend3.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

#==============================================================================
    def test_ruled(self):

        # sections for blend
        secs = [None]*4

        # first section is a NODE
        secs[0] = self.context.makeTopology(egads.NODE, reals=[0.5,0.5,-1])

        # FaceBody for 2nd section
        secs[1] = self._makeFaceBody(self.context)

        mat1 = ((1,  0,  0,  0), #Translate 1 in z-direction
                (0,  1,  0,  0),
                (0,  0,  1,  1))
        xform = self.context.makeTransform(mat1)

        # Translate 2nd section to make 3rd
        secs[2] = secs[1].copyObject(xform)

        # Finish with a NODE
        secs[3] = self.context.makeTopology(egads.NODE, reals=[0.5,0.5,2])

        ruled0 = egads.ruled(secs)
        data = ruled0.getInfo()
        self.assertEqual((egads.BODY, egads.SOLIDBODY), data[:2])

#==============================================================================
    def test_effective(self):

        # create the Circle
        rdata = [None]*10
        rdata[0] = 0 # center
        rdata[1] = 0
        rdata[2] = 0
        rdata[3] = 1   # x-axis
        rdata[4] = 0
        rdata[5] = 0
        rdata[6] = 0   # y-axis
        rdata[7] = 1
        rdata[8] = 0
        rdata[9] = 1   # radius

        circle = self.context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)

        # make Nodes
        nodes = [None]*4
        nodes[0] = self.context.makeTopology(egads.NODE, reals=[ 1, 0, 0])
        nodes[1] = self.context.makeTopology(egads.NODE, reals=[ 0, 1, 0])
        nodes[2] = self.context.makeTopology(egads.NODE, reals=[-1, 0, 0])
        nodes[3] = self.context.makeTopology(egads.NODE, reals=[ 0,-1, 0])

        edges = [None]*4
        tdata = [0, math.pi/2]
        edges[0] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[0:2], reals=tdata)

        tdata = [math.pi/2, math.pi]
        edges[1] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[1:3], reals=tdata)

        tdata = [math.pi, 3*math.pi/2]
        edges[2] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[2:4], reals=tdata)

        tdata = [3*math.pi/2, 2*math.pi]
        edges[3] = self.context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=[nodes[3], nodes[0]], reals=tdata)

        # make the loop
        sens = [1,1,1,1]
        loop = self.context.makeTopology(egads.LOOP, egads.CLOSED, children=edges, senses=sens)

        # make the WireBody
        wirebody = self.context.makeTopology(egads.BODY, egads.WIREBODY, children=[loop])

        mat1 = ((1,  0,  0,  0), #Translate 1 in z-direction
                (0,  1,  0,  0),
                (0,  0,  1,  1))
        xform = self.context.makeTransform(mat1)

        # sections for ruled
        secs = [None]*4
        secs[0] = wirebody
        for i in range(1,4):
            secs[i] = secs[i-1].copyObject(xform)

        body = egads.ruled(secs)

        faces = body.getBodyTopos(egads.FACE)
        self.assertEqual(12, len(faces))

        faces[0].attributeAdd("faceMap", [1])
        faces[1].attributeAdd("faceMap", [1])

        faceList = faces[9:11]

        tess = body.makeTessBody([0.1, 0.1, 15])
        tess.name = "tess"

        ebody = tess.initEBody(10)
        ebody.name = "ebody0"

        eface = ebody.makeEFace(faceList)
        eface.name = "eface0"

        data = eface.getInfo()
        self.assertEqual((egads.EFACE, egads.SFORWARD), data[:2])

        efaces = ebody.makeAttrEFaces("faceMap")

        self.assertEqual(1, len(efaces))
        efaces[0].name = "eface1"
        data = efaces[0].getInfo()
        self.assertEqual((egads.EFACE, egads.SFORWARD), data[:2])

        # finish the EBody
        ebody.finishEBody()
        
        # check the hierarchy
        (oclass, mtype, geom, lim, eshells, sens) = ebody.getTopology()
        self.assertEqual(egads.EBODY, oclass)
        self.assertEqual(1, len(eshells))
        self.assertEqual(body, geom)

        for eshell in eshells:
            (oclass, mtype, geom, lim, efaces, sens) = eshell.getTopology()
            self.assertEqual(oclass, egads.ESHELL)
            self.assertEqual(10, len(efaces))
            
            for eface in efaces:
                (oclass, mtype, geom, lim, eloops, sens) = eface.getTopology()
                self.assertEqual(oclass, egads.EFACE)
                self.assertEqual(1, len(eloops))
                
                for eface in efaces:
                    (oclass, mtype, geom, lim, eloops, sens) = eface.getTopology()
                    self.assertEqual(oclass, egads.EFACE)
                    self.assertEqual(1, len(eloops))
                    
                    for eloop in eloops:
                        (oclass, mtype, geom, lim, eedges, sens) = eloop.getTopology()
                        self.assertEqual(oclass, egads.ELOOPX)

        # check getting the original FACES
        faces = ebody.getBodyTopos(egads.FACE)
        self.assertEqual(12, len(faces))

        # get the EFACES and check the mapping
        efaces = ebody.getBodyTopos(egads.EFACE)
        self.assertEqual(10, len(efaces))

        uvrange, per = efaces[0].getRange()

        face, uv = efaces[0].effectiveMap([(uvrange[0]+uvrange[1])/2, (uvrange[2]+uvrange[3])/2])

        # get the EEDGES and check the mapping
        eedges = ebody.getBodyTopos(egads.EEDGE)

        edgeList = eedges[0].effectiveEdgeList()

        trange, per = eedges[0].getRange()

        edge, t = eedges[0].effectiveMap((trange[0]+trange[1])/2)


        faces[6].attributeAdd("faceMap", [2])
        faces[7].attributeAdd("faceMap", [2])

        ebody = tess.initEBody(10)
        ebody.name = "ebody1"

        efaces = ebody.makeAttrEFaces("faceMap")

        self.assertEqual(2, len(efaces))
        efaces[0].name = "eface2"
        efaces[1].name = "eface3"
        data = efaces[0].getInfo()
        self.assertEqual((egads.EFACE, egads.SFORWARD), data[:2])

        ebody.finishEBody()
        
        # create a model with Effective Topology body
        children = [body, ebody]

        model = self.context.makeTopology(egads.MODEL, children=children)

        (oclass, mtype, geom, lim, cldrn, sens) = model.getTopology()

        self.assertEqual(children, cldrn)

if __name__ == '__main__':
    unittest.main()
