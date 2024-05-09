import unittest
import os
import math
import sys
import gc

class TestEGADSlite(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        from pyEGADS import egads, egads_common

        def make_box_model():
            
            context = egads.Context()
            context.setOutLevel(0)
    
            box = context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
    
            model = context.makeTopology(egads.MODEL, children=[box])
    
            cls.box_model = model.exportModel()

        make_box_model()

        def make_box_attr_model():

            context = egads.Context()
            context.setOutLevel(0)
    
            box = context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
    
            # ints
            attrInt = [1, 2]
            box.attributeAdd("ints", attrInt)
    
            # float
            attrFloat= [1.0, 2.0]
            box.attributeAdd("float", attrFloat)
    
            # string
            attrString = "attrString"
            box.attributeAdd("string", attrString)
    
            # csystem
            attrCsys = egads.csystem([1,0,0, 0,1,0, 0,0,1])
            box.attributeAdd("csystem", attrCsys)
    
            model = context.makeTopology(egads.MODEL, children=[box])
            
            cls.box_attr_model = model.exportModel()

        make_box_attr_model()

        def make_box_fullattr_model():
            
            context = egads.Context()
            context.setOutLevel(0)
    
            box = context.makeSolidBody(egads.BOX, [0,0,0, 1,1,1])
    
            # ints
            attrInt1 = [1, 2]
            box.attributeAddSeq("ints", attrInt1)
            attrInt2 = [3, 4]
            box.attributeAddSeq("ints", attrInt2)
    
            # float
            attrFloat1= [1.0, 2.0]
            box.attributeAddSeq("float", attrFloat1)
            attrFloat2= [3.0, 4.0]
            box.attributeAddSeq("float", attrFloat2)
    
            # string
            attrString1 = "attrString1"
            box.attributeAddSeq("string", attrString1)
            attrString2 = "attrString2"
            box.attributeAddSeq("string", attrString2)
    
            # csystem
            attrCsys1 = egads.csystem([1,0,0, 0,1,0, 0,0,1])
            box.attributeAddSeq("csystem", attrCsys1)
            attrCsys2 = egads.csystem([0,1,0, 1,0,0, 0,0,1])
            box.attributeAddSeq("csystem", attrCsys2)
    
            model = context.makeTopology(egads.MODEL, children=[box])
    
            cls.box_fullattr_model = model.exportModel()

        make_box_fullattr_model()

        def make_line_tol_model():

            context = egads.Context()
            context.setOutLevel(0)
    
            # construction data
            x0 = [0,0,0]
            x1 = [1,0.01,0] # this should create a loose tolerance
    
            # make Nodes
            nodes = [None]*2
            nodes[0] = context.makeTopology(egads.NODE, reals=x0)
            nodes[1] = context.makeTopology(egads.NODE, reals=x1)
    
            # Line data (point and direction)
            rdata = [None]*6
            rdata[0] = x0[0]
            rdata[1] = x0[1]
            rdata[2] = x0[2]
            rdata[3] = x1[0] - x0[0]
            rdata[4] = x1[1] - x0[1]
            rdata[5] = x1[2] - x0[2]
            tdata = [0, 1]
    
            line     = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
            edge     = context.makeTopology(egads.EDGE , egads.TWONODE, geom=line, children=nodes, reals=tdata)
            loop     = context.makeTopology(egads.LOOP , egads.OPEN    , children=[edge], senses=[egads.SFORWARD])
            wirebody = context.makeTopology(egads.BODY , egads.WIREBODY, children=[loop])
            model    = context.makeTopology(egads.MODEL, children=[wirebody])

            cls.line_tol_model = model.exportModel()
            
            cls.line_tol_model_edge_tol = edge.tolerance()
            cls.line_tol_model_edge_getTol = edge.tolerance()

        make_line_tol_model()

        def make_faceBody_model():
    
            context = egads.Context()
            context.setOutLevel(0)
    
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
    
            faceBody = context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
    
            model = context.makeTopology(egads.MODEL, children=[faceBody])

            cls.faceBody_model = model.exportModel()

        make_faceBody_model()

        def make_curvature_curve():
            
            context = egads.Context()
            context.setOutLevel(0)
    
            # create the Circle
            rdata = [None]*10
            rdata[0] = 0   # center
            rdata[1] = 0
            rdata[2] = 0
            rdata[3] = 1   # x-axis
            rdata[4] = 0
            rdata[5] = 0
            rdata[6] = 0   # y-axis
            rdata[7] = 1
            rdata[8] = 0
            rdata[9] = 1   # radius
    
            circle = context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
            node   = context.makeTopology(egads.NODE, reals=[1,0,0])
            edge   = context.makeTopology(egads.EDGE, egads.ONENODE, geom=circle, children=[node], reals=[0, 2*math.pi])
            loop   = context.makeTopology(egads.LOOP, egads.CLOSED, children=[edge], senses=[egads.SFORWARD])
            body   = context.makeTopology(egads.BODY, egads.WIREBODY, children=[loop])
            model  = context.makeTopology(egads.MODEL, children=[body])

            cls.curvature_curve_model = model.exportModel()

        make_curvature_curve()
        
        def make_curvature_surface():
            
            context = egads.Context()
            context.setOutLevel(0)
    
            r = 1 # radius
            xcent = [0,0,0]
    
            # create Cylinder
            rdata = [None]*13
            rdata[ 0] = xcent[0]   # center
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
    
            cylindrical = context.makeGeometry(egads.SURFACE, egads.CYLINDRICAL, rdata)
            
            (oclass, mtype, rvec, ivec, eref) = cylindrical.getGeometry()
    
            # get the axes for the cylinder
            dx = rvec[3:6]
            dy = rvec[6:9]
            dz = rvec[9:12]
    
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
            circle0 = context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
    
            # create the Circle curve for the top 
            rdata[0] = rvec[0] + dz[0]*r; # center 
            rdata[1] = rvec[1] + dz[1]*r;
            rdata[2] = rvec[2] + dz[2]*r;
            rdata[3] = dx[0];             # x-axis 
            rdata[4] = dx[1];
            rdata[5] = dx[2];
            rdata[6] = dy[0];             # y-axis 
            rdata[7] = dy[1];
            rdata[8] = dy[2];
            rdata[9] = r;                 # radius 
            circle1 = context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
            
            nodes = [None]*2
            edges = [None]*8
    
            # create the Node on the Base 
            x1 = [None]*3
            x1[0] = xcent[0] + dx[0]*r;
            x1[1] = xcent[1] + dx[1]*r;
            x1[2] = xcent[2] + dx[2]*r;
            nodes[0] = context.makeTopology(egads.NODE, reals=x1)
            
            # create the Node on the Top
            x2 = [None]*3
            x2[0] = xcent[0] + dx[0]*r + dz[0]*r;
            x2[1] = xcent[1] + dx[1]*r + dz[1]*r;
            x2[2] = xcent[2] + dx[2]*r + dz[2]*r;
            nodes[1] = context.makeTopology(egads.NODE, reals=x2)
    
            # create the Line (point and direction)
            rdata = [None]*6
            rdata[0] = x1[0];
            rdata[1] = x1[1];
            rdata[2] = x1[2];
            rdata[3] = x2[0] - x1[0];
            rdata[4] = x2[1] - x1[1];
            rdata[5] = x2[2] - x1[2];
            line = context.makeGeometry(egads.CURVE, egads.LINE, rdata)
            
            tdata = [None]*2
            
            # make the Edge on the Line
            tdata[0] = 0
            tdata[1] = (rdata[3]*rdata[3] + rdata[4]*rdata[4] + rdata[5]*rdata[5])**0.5
            
            TWOPI = 2*math.pi
            
            edges[0] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=line, children=nodes, reals=tdata)
            edges[1] = context.makeTopology(egads.EDGE, egads.ONENODE, geom=circle0, children=[nodes[0]], reals=[0,TWOPI])
            edges[2] = edges[0]; # repeat the line edge 
            edges[3] = context.makeTopology(egads.EDGE, egads.ONENODE, geom=circle1, children=[nodes[1]], reals=[0,TWOPI])
    
            rdata = [None]*4
    
            # create P-curves
            rdata[0] = 0.;    rdata[1] = 0.;    # u == 0 UMIN
            rdata[2] = 0.;    rdata[3] = 1.;
            edges[4+0] = context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
            
            rdata[0] = 0.;    rdata[1] = 0;     # v == 0 VMIN 
            rdata[2] = 1.;    rdata[3] = 0.;
            edges[4+1] = context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
            
            rdata[0] = TWOPI; rdata[1] = 0.;    # u == TWOPI UMA
            rdata[2] = 0.;    rdata[3] = 1.;
            edges[4+2] = context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
            
            rdata[0] = 0.;    rdata[1] = r;     # v == r VMAX */
            rdata[2] = 1.;    rdata[3] = 0.;
            edges[4+3] = context.makeGeometry(egads.PCURVE, egads.LINE, rdata)
    
            # make the loop and Face body 
            
            senses = [egads.SREVERSE, egads.SFORWARD, egads.SFORWARD, egads.SREVERSE]
            loop = context.makeTopology(egads.LOOP, egads.CLOSED, geom=cylindrical, children=edges, senses=senses)
            face = context.makeTopology(egads.FACE, egads.SFORWARD, geom=cylindrical, children=[loop])
            body = context.makeTopology(egads.BODY, egads.FACEBODY, children=[face])
            model  = context.makeTopology(egads.MODEL, children=[body])

            cls.curvature_surface_model = model.exportModel()

        make_curvature_surface()
        
        def make_effective_model():
            
            context = egads.Context()
            context.setOutLevel(1)
    
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
    
            circle = context.makeGeometry(egads.CURVE, egads.CIRCLE, rdata)
    
            # make Nodes
            nodes = [None]*4
            nodes[0] = context.makeTopology(egads.NODE, reals=[ 1, 0, 0])
            nodes[1] = context.makeTopology(egads.NODE, reals=[ 0, 1, 0])
            nodes[2] = context.makeTopology(egads.NODE, reals=[-1, 0, 0])
            nodes[3] = context.makeTopology(egads.NODE, reals=[ 0,-1, 0])
    
            edges = [None]*4
            tdata = [0, math.pi/2]
            edges[0] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[0:2], reals=tdata)
    
            tdata = [math.pi/2, math.pi]
            edges[1] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[1:3], reals=tdata)
    
            tdata = [math.pi, 3*math.pi/2]
            edges[2] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=nodes[2:4], reals=tdata)
    
            tdata = [3*math.pi/2, 2*math.pi]
            edges[3] = context.makeTopology(egads.EDGE, egads.TWONODE, geom=circle, children=[nodes[3], nodes[0]], reals=tdata)
    
            # make the loop
            sens = [1,1,1,1]
            loop = context.makeTopology(egads.LOOP, egads.CLOSED, children=edges, senses=sens)
    
            # make the WireBody
            wirebody = context.makeTopology(egads.BODY, egads.WIREBODY, children=[loop])
    
            mat1 = ((1,  0,  0,  0), #Translate 1 in z-direction
                    (0,  1,  0,  0),
                    (0,  0,  1,  1))
            xform = context.makeTransform(mat1)
    
            # sections for ruled
            secs = [None]*4
            secs[0] = wirebody
            for i in range(1,4):
                secs[i] = secs[i-1].copyObject(xform)
    
            body = egads.ruled(secs)
    
            faces = body.getBodyTopos(egads.FACE)
    
            faces[0].attributeAdd("faceMap", [1])
            faces[1].attributeAdd("faceMap", [1])
    
            faceList = faces[9:11]
    
            tess = body.makeTessBody([0.1, 0.1, 15])

            ebody = tess.initEBody(10)
    
            eface = ebody.makeEFace(faceList)
    
            efaces = ebody.makeAttrEFaces("faceMap")
    
            # finish the EBody
            ebody.finishEBody()
    
            model = context.makeTopology(egads.MODEL, children=[body, ebody])
            
            cls.effective_model = model.exportModel()
 
        make_effective_model()

        # DANGER!! This trick allows for egadslite to be 
        # imported after importing egads. If ANY egads Context/ego
        # instances have not been garbage collected you will get a
        # segfault.
        # Don't do it!
        del egads_common._egads
        
        # Reload egads_common to reset Context/ego classes
        # Don't do this!
        import importlib
        importlib.reload(egads_common)
        
        from pyEGADS import egadslite
        globals()["egadslite"] = egadslite
        

    # setUp a new context for each test
    def setUp(self):
        self.context = egadslite.Context()
        self.context.setOutLevel(0)

    # make sure the context is clean before it is closed
    def tearDown(self):
        oclass, mtype, topObj, prev, next = self.context.getInfo()
        while next:
            oclass, mtype, topObj, prev, next = next.getInfo()
            self.assertNotEqual(oclass, egadslite.TESSELLATION, "Context is not properly clean!")
        del self.context

#==============================================================================
    def test_revision(self):
        imajor, iminor, OCCrev = egadslite.revision()
        self.assertIsInstance(imajor,int)
        self.assertIsInstance(iminor,int)
        self.assertIsInstance(OCCrev,str)

#==============================================================================
    def test_OpenClose(self):
        # Test open/close a context

        # Test creating a reference to a context c_ego
        # This reference will not close the context
        ref_context = egadslite.c_to_py(self.context.py_to_c())

#==============================================================================
    def test_importModel(self):

        modellite = self.context.importModel(self.box_model)

        (oclass, mtype, geom, lim, children, sens) = modellite.getTopology()

#==============================================================================
    def test_getInfo(self):

        model = self.context.importModel(self.box_model)

        # getInfo
        (oclass, mtype, topObj, prev, next) = model.getInfo()
        self.assertEqual((egadslite.MODEL, 0), (oclass, mtype))

#==============================================================================
    def test_attribute(self):

        attrInt = [1, 2]
        attrFloat= [1.0, 2.0]
        attrString = "attrString"
        attrCsys = egadslite.csystem([1,0,0, 0,1,0, 0,0,1])

        model = self.context.importModel(self.box_attr_model)

        (oclass, mtype, geom, lim, children, sens) = model.getTopology()
        box = children[0]

        self.assertEqual(4, box.attributeNum())

        # ints
        data = box.attributeRet("ints")
        self.assertEqual( attrInt, data )
        
        name, data = box.attributeGet(1)
        self.assertEqual( "ints", name )
        self.assertEqual( attrInt, data )

        # float
        data = box.attributeRet("float")
        self.assertEqual( attrFloat, data )

        name, data = box.attributeGet(2)
        self.assertEqual( "float", name )
        self.assertEqual( attrFloat, data )

        # string
        data = box.attributeRet("string")
        self.assertEqual( attrString, data )

        name, data = box.attributeGet(3)
        self.assertEqual( "string", name )
        self.assertEqual( attrString, data )

        # csystem
        data = box.attributeRet("csystem")
        self.assertEqual( attrCsys, data )
        self.assertIsInstance( data, egadslite.csystem )

        name, data = box.attributeGet(4)
        self.assertEqual( "csystem", name )
        self.assertEqual( attrCsys, data )
        self.assertIsInstance( data, egadslite.csystem )

#==============================================================================
    def test_fullAttribute(self):

        # ints
        attrInt1 = [1, 2]
        attrInt2 = [3, 4]

        # float
        attrFloat1= [1.0, 2.0]
        attrFloat2= [3.0, 4.0]

        # string
        attrString1 = "attrString1"
        attrString2 = "attrString2"

        # csystem
        attrCsys1 = egadslite.csystem([1,0,0, 0,1,0, 0,0,1])
        attrCsys2 = egadslite.csystem([0,1,0, 1,0,0, 0,0,1])

        model = self.context.importModel(self.box_fullattr_model)

        (oclass, mtype, geom, lim, children, sens) = model.getTopology()
        box = children[0]

        self.assertEqual(8, box.attributeNum())
        self.assertEqual(2, box.attributeNumSeq("ints"))
        self.assertEqual(2, box.attributeNumSeq("float"))
        self.assertEqual(2, box.attributeNumSeq("string"))
        self.assertEqual(2, box.attributeNumSeq("csystem"))

        # ints
        data = box.attributeRetSeq("ints", 1)
        self.assertEqual( attrInt1, data )
        data = box.attributeRetSeq("ints", 2)
        self.assertEqual( attrInt2, data )
        
        name, data = box.attributeGet(1)
        self.assertEqual( "ints 1", name )
        self.assertEqual( attrInt1, data )

        name, data = box.attributeGet(2)
        self.assertEqual( "ints 2", name )
        self.assertEqual( attrInt2, data )

        # float
        data = box.attributeRetSeq("float", 1)
        self.assertEqual( attrFloat1, data )
        data = box.attributeRetSeq("float", 2)
        self.assertEqual( attrFloat2, data )

        name, data = box.attributeGet(3)
        self.assertEqual( "float 1", name )
        self.assertEqual( attrFloat1, data )
        
        name, data = box.attributeGet(4)
        self.assertEqual( "float 2", name )
        self.assertEqual( attrFloat2, data )

        # string
        data = box.attributeRetSeq("string", 1)
        self.assertEqual( attrString1, data )

        data = box.attributeRetSeq("string", 2)
        self.assertEqual( attrString2, data )

        name, data = box.attributeGet(5)
        self.assertEqual( "string 1", name )
        self.assertEqual( attrString1, data )

        name, data = box.attributeGet(6)
        self.assertEqual( "string 2", name )
        self.assertEqual( attrString2, data )

        # csystem
        data = box.attributeRetSeq("csystem", 1)
        self.assertEqual( attrCsys1, data )
        self.assertIsInstance( data, egadslite.csystem )

        data = box.attributeRetSeq("csystem", 2)
        self.assertEqual( attrCsys2, data )
        self.assertIsInstance( data, egadslite.csystem )

        name, data = box.attributeGet(7)
        self.assertEqual( "csystem 1", name )
        self.assertEqual( attrCsys1, data )
        self.assertIsInstance( data, egadslite.csystem )

        name, data = box.attributeGet(8)
        self.assertEqual( "csystem 2", name )
        self.assertEqual( attrCsys2, data )
        self.assertIsInstance( data, egadslite.csystem )

#==============================================================================
    def test_topology(self):

        model = self.context.importModel(self.box_model)

        # check the hierarchy
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        self.assertEqual(egadslite.MODEL, oclass)
        self.assertEqual(1, len(bodies))

        for body in bodies:
            (oclass, mtype, geom, lim, shells, sens) = body.getTopology()
            self.assertEqual(egadslite.BODY, oclass)
            self.assertEqual(1, len(shells))
    
            for shell in shells:
                (oclass, mtype, geom, lim, faces, sens) = shell.getTopology()
                self.assertEqual(oclass, egadslite.SHELL)
                self.assertEqual(6, len(faces))
                
                for face in faces:
                    (oclass, mtype, geom, lim, loops, sens) = face.getTopology()
                    self.assertEqual(oclass, egadslite.FACE)
                    self.assertEqual(1, len(loops))
                    
                    for loop in loops:
                        (oclass, mtype, geom, lim, edges, sens) = loop.getTopology()
                        self.assertEqual(oclass, egadslite.LOOP)
                        self.assertEqual(4, len(edges))

                        for edge in edges:
                            (oclass, mtype, geom, lim, nodes, sens) = edge.getTopology()
                            self.assertEqual(oclass, egadslite.EDGE)
                            self.assertEqual(2, len(nodes))

        box = bodies[0]

        bodyNodes = box.getBodyTopos(egadslite.NODE)
        bodyEdges = box.getBodyTopos(egadslite.EDGE)
        bodyFaces = box.getBodyTopos(egadslite.FACE)

        for face in bodyFaces:
            self.assertEqual( face.getBody(), box )

        for i in range(len(bodyNodes)):
            self.assertEqual( face.getBody(), box )
            self.assertEqual( box.indexBodyTopo(bodyNodes[i]), i+1 )
            self.assertEqual( box.objectBodyTopo(egadslite.NODE, i+1), bodyNodes[i] )

        for i in range(len(bodyEdges)):
            self.assertEqual( box.indexBodyTopo(bodyEdges[i]), i+1 )
            self.assertEqual( box.objectBodyTopo(egadslite.EDGE, i+1), bodyEdges[i] )

#==============================================================================
    def test_tolerance(self):

        model = self.context.importModel(self.line_tol_model)
        
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        edges = bodies[0].getBodyTopos(egadslite.EDGE)

        # tolerance
        self.assertEqual(edges[0].tolerance()   , self.line_tol_model_edge_tol   )
        self.assertEqual(edges[0].getTolerance(), self.line_tol_model_edge_getTol)

#==============================================================================
    def test_properties(self):

        model = self.context.importModel(self.faceBody_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        faceBody = bodies[0]

        nodes = faceBody.getBodyTopos(egadslite.NODE)
        edges = faceBody.getBodyTopos(egadslite.EDGE)
        faces = faceBody.getBodyTopos(egadslite.FACE)

        (oclass, mtype, line, lim, children, sens) = edges[0].getTopology()

        (trange, periodic) = edges[0].getRange()

        # length
        self.assertEqual(line.arcLength(trange[0],trange[1]), 1.0)

        (oclass, mtype, plane, lim, children, sens) = faces[0].getTopology()

        (urange, periodic) = faces[0].getRange()

        # area
        #self.assertEqual(faces[0].getArea(), 1.0)
        #self.assertEqual(plane.getArea(urange), 1.0)

        # bounding box
        bbox = faceBody.getBoundingBox()
        self.assertEqual(bbox, [(0,0,0), (1,1,0)])

#==============================================================================
    def test_curvature_curve(self):

        # curvature

        # Only one model is allowed
        model = self.context.importModel(self.curvature_curve_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        (oclass, mtype, geom, lim, loops, sens)  = bodies[0].getTopology()
        (oclass, mtype, geom, lim, edges, sens)  = loops[0].getTopology()
        edge = edges[0]

        data = edge.curvature(0)
        self.assertEqual((1.0, (0.0, 1.0, 0.0)), data)

#==============================================================================
    def test_curvature_surface(self):

        # Only one model is allowed
        model = self.context.importModel(self.curvature_surface_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        face = bodies[0].getBodyTopos(egadslite.FACE)[0]

        data = face.curvature([0,0])
        self.assertEqual((0.0, (0.0, 0.0, 1.0), -1.0, (0.0, 1.0, 0.0)), data)

#==============================================================================
    def test_evaluate(self):

        model = self.context.importModel(self.faceBody_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        faceBody = bodies[0]

        nodes = faceBody.getBodyTopos(egadslite.NODE)
        edges = faceBody.getBodyTopos(egadslite.EDGE)
        faces = faceBody.getBodyTopos(egadslite.FACE)

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

        model = self.context.importModel(self.faceBody_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        faceBody = bodies[0]

        nodes = faceBody.getBodyTopos(egadslite.NODE)
        edges = faceBody.getBodyTopos(egadslite.EDGE)
        faces = faceBody.getBodyTopos(egadslite.FACE)

        inveval = edges[0].invEvaluate([0.5, 0., 0.])
        self.assertEqual(inveval, (0.5,          # parameter
                                  [0.5,0.0,0.])) # coordinate

        inveval = faces[0].invEvaluate([0.5, 0.5, 0.])
        self.assertEqual(inveval, ([0.5, 0.5],    # parameter
                                   [0.5,0.5,0.])) # coordinate

#==============================================================================
    def test_getEdgeUV(self):

        model = self.context.importModel(self.faceBody_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        faceBody = bodies[0]

        faces = faceBody.getBodyTopos(egadslite.FACE)
        edges = faceBody.getBodyTopos(egadslite.EDGE)

        uv = faces[0].getEdgeUV(edges[0], 0.5)
        self.assertEqual(uv, (0.5, 0.))

#==============================================================================
    def test_tessBody(self):

        model = self.context.importModel(self.box_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        box0 = bodies[0]

        tess = box0.makeTessBody([0.1, 0.1, 30.])

        xyz, t = tess.getTessEdge(1)
        lIndex = tess.getTessLoops(1)

        body, stat, icode, npts = tess.statusTessBody()

        xyz_global = [None]*npts
        for j in range(npts):
            ptype, pindex, xyz_global[j] = tess.getGlobal(j+1)

        tris_global = []
        nfaces = len(box0.getBodyTopos(egadslite.FACE))
        for fIndex in range(nfaces):
            xyz, uv, ptype, pindex, tris, tric = tess.getTessFace(fIndex+1)

            for j in range(len(tris)):
                for n in range(3):
                    tris_global.append(tess.localToGlobal(fIndex+1, tris[j][n]))

        # Build up a tesselleation (by copying)
        tess2 = box0.initTessBody()

        nedges = len(box0.getBodyTopos(egadslite.EDGE))
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
        edges = box0.getBodyTopos(egadslite.EDGE)
        facedg = box0.getBodyTopos(egadslite.FACE, edges[0])
        facedg.append( edges[0] )

        # Renamke the tessellation
        tess.remakeTess(facedg, [0.2,0.1,10])

        #-------
        # crate a quad tessellation
        quadTess = tess.quadTess()

        #-------
        # mass properties
        # volume, area, CG, I = quadTess.tessMassProps()
        # self.assertAlmostEqual(volume, 1.0, 9)
        # self.assertAlmostEqual(area, 6.0, 9)
        # for i in range(3):
        #     self.assertAlmostEqual(CG[i], (0.5, 0.5, 0.5)[i], 9)
        # Itrue = (1/6.,   0.,   0.,
        #            0., 1/6.,   0.,
        #            0.,   0., 1/6.)
        # for i in range(9):
        #     self.assertAlmostEqual(I[i], Itrue[i], 9)


#==============================================================================
    def test_tessGeom(self):

        model = self.context.importModel(self.faceBody_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        faceBody = bodies[0]

        faces = faceBody.getBodyTopos(egadslite.FACE)

        (oclass, mtype, plane, lim, loops, sens) = faces[0].getTopology()
        
        tess = plane.makeTessGeom([0,1, 0,1], [5,5])

        sizes, xyz = tess.getTessGeom()
        self.assertEqual(sizes, (5,5))

#==============================================================================
    def test_getWindingAngle(self):

        model = self.context.importModel(self.box_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        box0 = bodies[0]

        edges = box0.getBodyTopos(egadslite.EDGE)
        
        trange, per = edges[0].getRange()
        t = (trange[0] + trange[1])/2
        
        angle = edges[0].getWindingAngle(t)
        self.assertAlmostEqual(270.0, angle, 5)

#==============================================================================
    def test_inTopology(self):

        model = self.context.importModel(self.box_model)
        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        box0 = bodies[0]

        faces = box0.getBodyTopos(egadslite.FACE)
        uvrange, periodic = faces[0].getRange()
        bbox = faces[0].getBoundingBox()

        self.assertTrue( faces[0].inTopology([(bbox[0][0]+bbox[1][0])/2, (bbox[0][1]+bbox[1][1])/2, (bbox[0][2]+bbox[1][2])/2] ) )

        self.assertTrue( faces[0].inFace([(uvrange[0]+uvrange[1])/2, (uvrange[2]+uvrange[3])/2] ) )

#==============================================================================
    def test_effective(self):

        model = self.context.importModel(self.effective_model)

        (oclass, mtype, geom, lim, bodies, sens) = model.getTopology()
        body = bodies[0]
        ebody = bodies[1]

        # check the hierarchy
        (oclass, mtype, geom, lim, eshells, sens) = ebody.getTopology()
        self.assertEqual(egadslite.EBODY, oclass)
        self.assertEqual(1, len(eshells))
        self.assertEqual(body, geom)

        for eshell in eshells:
            (oclass, mtype, geom, lim, efaces, sens) = eshell.getTopology()
            self.assertEqual(oclass, egadslite.ESHELL)
            self.assertEqual(10, len(efaces))
            
            for eface in efaces:
                (oclass, mtype, geom, lim, eloops, sens) = eface.getTopology()
                self.assertEqual(oclass, egadslite.EFACE)
                self.assertEqual(1, len(eloops))
                
                for eface in efaces:
                    (oclass, mtype, geom, lim, eloops, sens) = eface.getTopology()
                    self.assertEqual(oclass, egadslite.EFACE)
                    self.assertEqual(1, len(eloops))
                    
                    for eloop in eloops:
                        (oclass, mtype, geom, lim, eedges, sens) = eloop.getTopology()
                        self.assertEqual(oclass, egadslite.ELOOPX)

        # check getting the original FACES
        faces = ebody.getBodyTopos(egadslite.FACE)
        self.assertEqual(12, len(faces))

        # get the EFACES and check the mapping
        efaces = ebody.getBodyTopos(egadslite.EFACE)
        self.assertEqual(10, len(efaces))

        uvrange, per = efaces[0].getRange()

        face, uv = efaces[0].effectiveMap([(uvrange[0]+uvrange[1])/2, (uvrange[2]+uvrange[3])/2])

        # get the EEDGES and check the mapping
        eedges = ebody.getBodyTopos(egadslite.EEDGE)

        edgeList = eedges[0].effectiveEdgeList()

        trange, per = eedges[0].getRange()

        edge, t = eedges[0].effectiveMap((trange[0]+trange[1])/2)

if __name__ == '__main__':
    unittest.main()
