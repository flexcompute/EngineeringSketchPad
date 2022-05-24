import unittest

import os
import glob
import shutil
import math

try:
    import numpy
except ImportError:
    numpy = None

from sys import version_info as pyVersion

import pyCAPS
from pyEGADS import egads

class TestGeometry(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.problemName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        cls.myGeometry = cls.myProblem.geometry

    @classmethod
    def tearDownClass(cls):
        del cls.myGeometry
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

        if os.path.isfile("myGeometry.egads"):
            os.remove("myGeometry.egads")

        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

        if os.path.isfile("myGeometry.html"):
            os.remove("myGeometry.html")

        if os.path.isfile("geomView_0.png"):
            os.remove("geomView_0.png")

        if os.path.isfile("geomView_1.png"):
            os.remove("geomView_1.png")

        if os.path.isfile("geomView_2.png"):
            os.remove("geomView_2.png")

#==============================================================================
    # Check geometry parameters
    def test_geomPmtr(self):

        myGeometry = self.myProblem.geometry

        # Get the names of all OpenCSM des/cfg/out pmtr in unitGeom.csm
        self.assertEqual(sorted(myGeometry.despmtr.keys()), sorted(['wing:dihedral', 'taper', 'wing:lesweep', 'wing:chord:root', 'despMat', 'despRow', 'despCol', 'area', 'twist', 'aspect', 'v@1:d_name', 'htail', 'htail:chord', 'vtail', 'vtail:chord', 'sphereR']))
        self.assertEqual(sorted(myGeometry.cfgpmtr.keys()), sorted(['VIEW:CFD', 'nrow', 'ncol', 'series2', 'series']))
        self.assertEqual(sorted(myGeometry.conpmtr.keys()), sorted(['nConst']))
        self.assertEqual(sorted(myGeometry.outpmtr.keys()), sorted(['dummyRow', 'dummyRow3', 'dummyRow2', 'span', 'cmean', 'dummyCol', 'dummyMat', 'sphereVol', 'boxVol', 'boxCG']))

        # Access despmtr Value Object
        taperObj = myGeometry.despmtr["taper"]
        self.assertEqual(0.5, taperObj.value)
        taperObj.value = 0.3
        self.assertEqual(0.3, taperObj.value)
        self.assertEqual([0,1], taperObj.limits)
        self.assertEqual(taperObj.props(), {"dim":0, "pmtr":0, "lfix":1, "sfix":1, "ntype":0})

        # new attributes cannot be added to the Value Object
        with self.assertRaises(AttributeError):
            taperObj.values = 0.8

        # cannot access non-existing taper2 despmtr
        with self.assertRaises(KeyError):
            taper2Obj = myGeometry.despmtr["taper2"]

        # cannot set value out of bound
        with self.assertRaises(pyCAPS.CAPSError) as e:
            taperObj.value = 1.8
        self.assertEqual(e.exception.errorName, "CAPS_RANGEERR")

        # cannot delete Value Object
        with self.assertRaises(pyCAPS.CAPSError):
            del myGeometry.despmtr["taper"]

        # cannot delete value shortcut
        with self.assertRaises(pyCAPS.CAPSError):
            del myGeometry.despmtr.taper

        # Shortcut for accessing des/cfg pmtr value
        self.assertEqual(0.3, myGeometry.despmtr.taper)
        myGeometry.despmtr.taper = 0.4
        self.assertEqual(0.4, myGeometry.despmtr.taper)
        myGeometry.despmtr.taper += 0.1
        self.assertEqual(0.5, myGeometry.despmtr.taper)


        myGeometry.despmtr.wing.dihedral = 8
        self.assertEqual( 8, myGeometry.despmtr.wing.dihedral)
        self.assertEqual( 8, myGeometry.despmtr["wing:dihedral"].value)
        self.assertEqual( 8, myGeometry.despmtr["wing"].dihedral)
        self.assertEqual( 8, myGeometry.despmtr["wing"]["dihedral"].value)

        myGeometry.despmtr.wing.lesweep = 10
        self.assertEqual(10, myGeometry.despmtr.wing.lesweep)
        self.assertEqual(sorted(myGeometry.despmtr.wing.keys()), sorted(['lesweep', 'dihedral', 'chord:root']))

        self.myGeometry.despmtr.wing.lesweep = 20.0
        self.assertEqual(self.myGeometry.despmtr.wing.lesweep, 20.0)

        self.myGeometry.despmtr.wing.lesweep = [30.0]
        self.assertEqual(self.myGeometry.despmtr.wing.lesweep, 30.0)

        self.myGeometry.despmtr.wing.lesweep = [[40.0]]
        self.assertEqual(self.myGeometry.despmtr.wing.lesweep, 40.0)

        myGeometry.despmtr["wing:chord:root"].value = 1.5
        self.assertEqual( 1.5, myGeometry.despmtr["wing:chord:root"].value)
        self.assertEqual( 1.5, myGeometry.despmtr["wing"]["chord:root"].value)
        self.assertEqual( 1.5, myGeometry.despmtr["wing"]["chord"]["root"].value)
        self.assertEqual( 1.5, myGeometry.despmtr["wing"]["chord"].root)
        self.assertEqual( 1.5, myGeometry.despmtr["wing"].chord.root)
        self.assertEqual( 1.5, myGeometry.despmtr.wing.chord.root)

        myGeometry.despmtr["htail"].value = 2
        myGeometry.despmtr["htail:chord"].value = 3
        self.assertEqual( 2, myGeometry.despmtr["htail"].value)
        self.assertEqual( 2, myGeometry.despmtr.htail)
        self.assertEqual( 3, myGeometry.despmtr["htail:chord"].value)

        myGeometry.despmtr["vtail"].value = 4
        myGeometry.despmtr["vtail:chord"].value = 5
        self.assertEqual( 4, myGeometry.despmtr["vtail"].value)
        self.assertEqual( 4, myGeometry.despmtr.vtail)
        self.assertEqual( 5, myGeometry.despmtr["vtail:chord"].value)

        self.myGeometry.cfgpmtr.VIEW.CFD = False
        self.assertEqual(self.myGeometry.cfgpmtr.VIEW.CFD, 0)
        self.myGeometry.cfgpmtr.VIEW.CFD = True
        self.assertEqual(self.myGeometry.cfgpmtr.VIEW.CFD, 1)


        self.assertEqual(8412, myGeometry.cfgpmtr.series)
        myGeometry.cfgpmtr.series = 2424
        self.assertEqual(2424, myGeometry.cfgpmtr.series)

        # cannot get non-existing des/cfg pmtr
        with self.assertRaises(AttributeError):
            series3 = myGeometry.cfgpmtr.series3

        # cannot set non-existing des/cfg pmtr
        with self.assertRaises(AttributeError):
            myGeometry.cfgpmtr.series3 = 1212


        self.assertEqual(10, myGeometry.conpmtr.nConst)
        # cannot set conpmtr
        with self.assertRaises(pyCAPS.CAPSError):
            myGeometry.conpmtr.nConst = 10

        # compute some outpmtr values consistent with equations in unitGeom.csm
        aspect = myGeometry.despmtr.aspect
        area   = myGeometry.despmtr.area
        cmean  = (area/aspect)**0.5
        span   = cmean*aspect

        # Access outpmtr class
        spanObj = myGeometry.outpmtr["span"]
        self.assertAlmostEqual(span, spanObj.value, 5)
        self.assertEqual(spanObj.props(), {"dim":0, "pmtr":0, "lfix":0, "sfix":0, "ntype":1})

        # outpmtr value cannot be modified
        with self.assertRaises(pyCAPS.CAPSError):
            myGeometry.outpmtr["cmean"].value = 5

        # Shortcut for accessing outpmtr value
        self.assertAlmostEqual(cmean, myGeometry.outpmtr.cmean, 5)
        self.assertAlmostEqual(cmean, myGeometry.outpmtr["cmean"].value, 5)

        # outpmtr value cannot be modified
        with self.assertRaises(pyCAPS.CAPSError):
            myGeometry.outpmtr.cmean = 5

        # cannot get non-existing outpmtr
        with self.assertRaises(KeyError):
            series3 = myGeometry.outpmtr["series3"]
        with self.assertRaises(AttributeError):
            series3 = myGeometry.outpmtr.series3

        # cannot set non-existing outpmtr
        with self.assertRaises(AttributeError):
            myGeometry.outpmtr.series3 = 1212

        # Check outpmtr matrix access
        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [15.0, 16.0]])
        myGeometry.despmtr.despMat = [[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]]
        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]])

        myGeometry.despmtr.despMat = [(12.0, 13.0), (13.0, 14.0), (30.0, 32.0)]
        self.assertEqual(myGeometry.despmtr.despMat, [[12.0, 13.0], [13.0, 14.0], [30.0, 32.0]])

        # Modyifying temporay variable does not change the value
        myGeometry.despmtr.despMat[0] = [2., 3.]
        self.assertEqual(myGeometry.despmtr.despMat, [[12.0, 13.0], [13.0, 14.0], [30.0, 32.0]])

        myGeometry.despmtr.despCol = [[10.],[11.],[12.]]
        self.assertEqual(myGeometry.despmtr.despCol, [10.,11.,12.])

        myGeometry.despmtr.despRow = [11.,12.,13.]
        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.,13.])

        self.assertEqual(myGeometry.outpmtr.dummyRow ,[1,2,3])
        self.assertEqual(myGeometry.outpmtr.dummyRow2,[1,2,None])
        self.assertEqual(myGeometry.outpmtr.dummyRow3,[2,3,1,1,1])
        self.assertEqual(myGeometry.outpmtr.dummyCol ,[3,2,1])
        self.assertEqual(myGeometry.outpmtr.dummyMat ,[[1,2],[3,4],[5,6]])

        # Test None values in matrix
        myGeometry.despmtr.despMat = [[None, None], [None, 14.0], [15.0, None]]
        self.assertEqual(myGeometry.despmtr.despMat, [[None, None], [None, 14.0], [15.0, None]])
        myGeometry.despmtr.despMat = [[11.0, 12.0], [13.0, 14.0], [15.0, 16.0]]
        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [15.0, 16.0]])

        # Get just the names of the output values
        valueList = self.myGeometry.outpmtr.keys()
        self.assertEqual("cmean" in valueList, True)

#==============================================================================
    # Test geometry values with numpy
    @unittest.skipIf(numpy is None, "numpy not installed")
    def test_numpy(self):

        myGeometry = self.myProblem.geometry

        myGeometry.cfgpmtr.series = numpy.int32(2425)
        self.assertEqual(2425, myGeometry.cfgpmtr.series)

        myGeometry.cfgpmtr.series = numpy.int64(2426)
        self.assertEqual(2426, myGeometry.cfgpmtr.series)

        myGeometry.cfgpmtr.series = numpy.float32(2427)
        self.assertEqual(2427, myGeometry.cfgpmtr.series)

        myGeometry.cfgpmtr.series = numpy.float64(2428)
        self.assertEqual(2428, myGeometry.cfgpmtr.series)

        myGeometry.despmtr.despMat = numpy.array([[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]], numpy.float32)
        self.assertEqual(myGeometry.despmtr.despMat, [[11.0, 12.0], [13.0, 14.0], [30.0, 32.0]])

        myGeometry.despmtr.despMat = numpy.array([(12.0, 13.0), (13.0, 14.0), (30.0, 32.0)], numpy.float64)
        self.assertEqual(myGeometry.despmtr.despMat, [[12.0, 13.0], [13.0, 14.0], [30.0, 32.0]])

        # Modyifying temporay variable does not change the value
        myGeometry.despmtr.despMat[0] = numpy.array([2., 3.], numpy.float32)
        self.assertEqual(myGeometry.despmtr.despMat, [[12.0, 13.0], [13.0, 14.0], [30.0, 32.0]])

        myGeometry.despmtr.despCol = numpy.array([[10],[11],[12]], numpy.int32)
        self.assertEqual(myGeometry.despmtr.despCol, [10.,11.,12.])

        myGeometry.despmtr.despRow = numpy.array([11,12,13], numpy.int64)
        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.,13.])

#==============================================================================
    # Retrieve geometry attributes
    def test_matShapeChange(self):

        myGeometry = self.myProblem.geometry

        #---------
        self.assertEqual(myGeometry.cfgpmtr.nrow, 3)
        self.assertEqual(myGeometry.cfgpmtr.ncol, 2)
        myGeometry.despmtr.despMat = [[11.,12.],
                                      [13.,14.],
                                      [15.,16.]]

        self.assertEqual(myGeometry.despmtr.despMat, [[11.,12.],
                                                      [13.,14.],
                                                      [15.,16.]])

        #---------
        myGeometry.cfgpmtr.nrow = 5
        self.assertEqual(myGeometry.cfgpmtr.nrow, 5)
        myGeometry.despmtr.despMat = [[11.,12.],
                                      [13.,14.],
                                      [15.,16.],
                                      [17.,18.],
                                      [19.,20.]]
        self.assertEqual(myGeometry.despmtr.despMat, [[11.,12.],
                                                      [13.,14.],
                                                      [15.,16.],
                                                      [17.,18.],
                                                      [19.,20.]])
        myGeometry.despmtr.despCol = [11.,12.,13.,14.,15.]
        self.assertEqual(myGeometry.despmtr.despCol, [11.,12.,13.,14.,15.])
        myGeometry.despmtr.despRow = [11.,12.,13.,14.,15.]
        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.,13.,14.,15.])

        #---------
        myGeometry.cfgpmtr.nrow = 1
        myGeometry.cfgpmtr.ncol = 2
        self.assertEqual(myGeometry.cfgpmtr.nrow, 1)
        self.assertEqual(myGeometry.cfgpmtr.ncol, 2)
        self.assertEqual(myGeometry.despmtr.despMat, [11.,12.])
        myGeometry.despmtr.despMat = [13.,14.]
        self.assertEqual(myGeometry.despmtr.despMat, [13.,14.])

        myGeometry.despmtr.despCol = 11.
        self.assertEqual(myGeometry.despmtr.despCol, 11.)
        myGeometry.despmtr.despRow = 11.
        self.assertEqual(myGeometry.despmtr.despRow, 11.)

        #---------
        myGeometry.cfgpmtr.nrow = 2
        myGeometry.cfgpmtr.ncol = 1
        self.assertEqual(myGeometry.cfgpmtr.nrow, 2)
        self.assertEqual(myGeometry.cfgpmtr.ncol, 1)
        self.assertEqual(myGeometry.despmtr.despMat, [13.,14.])
        myGeometry.despmtr.despMat = [11.,12.]
        self.assertEqual(myGeometry.despmtr.despMat, [11.,12.])

        myGeometry.despmtr.despCol = [11.,12.]
        self.assertEqual(myGeometry.despmtr.despCol, [11.,12.])
        myGeometry.despmtr.despRow = [11.,12.]
        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.])

        #---------
        myGeometry.cfgpmtr.nrow = 3
        myGeometry.cfgpmtr.ncol = 2
        self.assertEqual(myGeometry.cfgpmtr.nrow, 3)
        self.assertEqual(myGeometry.cfgpmtr.ncol, 2)
        self.assertEqual(myGeometry.despmtr.despMat, [[11.,12.],
                                                      [12.,12.],
                                                      [12.,12.]])
        self.assertEqual(myGeometry.despmtr.despCol, [11.,12.,12.])
        self.assertEqual(myGeometry.despmtr.despRow, [11.,12.,12.])


#==============================================================================
    # Retrieve geometry attributes
    def test_geomAttrList(self):

        # Make a new problem and load the egads file
        problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=self.file, outLevel=0 ); self.__class__.iProb += 1

        attributeList = problem.geometry.attrList("capsReferenceArea")
        self.assertEqual(attributeList, [40])

        attributeList = problem.geometry.attrList("capsGroup", attrLevel="Body")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

        attributeList = problem.geometry.attrList("capsGroup", attrLevel="Face")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

        attributeList = problem.geometry.attrList("capsGroup", attrLevel="Node")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

        # Check CSYSTEM attribute on body
        attributeList = problem.geometry.attrList("leftWingSkin", attrLevel="Body")
        self.assertEqual(len(attributeList), 1)

        attributeList = attributeList[0] # Only want the first element
        self.assertIsInstance(attributeList, egads.csystem)

        self.assertEqual(len(attributeList), 9)

        leftWingSkin = [ 3.77246051677845,  0               ,  0.004591314188236,
                        10.92812876825235, 1.941581372589092, 26.203935380022607,
                         2.13241446403669, -7.07106781186547, -0.365374367002283]
        for i in range(9):
            self.assertAlmostEqual(attributeList[i], leftWingSkin[i], 5)

#==============================================================================
    # Retrieve geometry attribute maps
    def test_geomAttrMap(self):

        attributeMap = self.myGeometry.attrMap()
        Body1 = attributeMap['Wing1']
        Body2 = attributeMap['Wing2']

        keys = list(Body1['Body'].keys())
        self.assertEqual(sorted(keys), sorted(['capsGroup', 'capsReferenceArea', '.tParam', 'capsIntent', 'capsLength', 'capsAIM', 'capsReferenceSpan', 'leftWingSkin', 'riteWingSkin', 'capsReferenceChord']))

        keys = list(Body1['Faces'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8])

        self.assertEqual(Body1['Faces'][1]['capsGroup'], 'Wing1')

        keys = list(Body1['Edges'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

        keys = list(Body1['Nodes'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8, 9])

        keys = list(Body2['Body'].keys())
        self.assertEqual(sorted(keys), sorted(['capsGroup', 'capsReferenceArea', '.tParam', 'capsIntent', 'capsLength', 'capsAIM', 'capsReferenceSpan', 'capsReferenceChord']))

        keys = list(Body2['Faces'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8])

        keys = list(Body2['Edges'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15])

        keys = list(Body2['Nodes'].keys())
        self.assertEqual(keys, [1, 2, 3, 4, 5, 6, 7, 8, 9])

        # check the reverse map
        attributeMap = self.myGeometry.attrMap(reverseMap = True)
        Body1 = attributeMap['Wing1']
        Body2 = attributeMap['Wing2']
        Body3 = attributeMap['Farfield']

        keys = list(Body1.keys())
        self.assertEqual(sorted(keys), sorted(['capsBound','capsGroup', 'capsReferenceArea', '.tParam', 'capsIntent', 'capsLength', 'capsAIM', 'capsReferenceSpan', 'leftWingSkin', 'riteWingSkin', 'capsReferenceChord']))

        Wing1 = Body1['capsGroup']['Wing1']
        self.assertEqual(Wing1, {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': [1, 2, 3, 4, 5, 6, 7, 8]})

        capsAIM = Body1['capsAIM']
        self.assertEqual(capsAIM, {'fun3dAIM;su2AIM;egadsTessAIM;aflr4AIM;tetgenAIM;aflr3AIM;masstranAIM': {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': []}})

        Wing2 = Body2['capsGroup']['Wing2']
        self.assertEqual(Wing2, {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': [1, 2, 3, 4, 5, 6, 7, 8]})

        capsAIM = Body2['capsAIM']
        self.assertEqual(capsAIM, {'fun3dAIM;su2AIM;egadsTessAIM;aflr4AIM;tetgenAIM;aflr3AIM;masstranAIM': {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': []}})

        Farfield = Body3['capsGroup']['Farfield']
        self.assertEqual(Farfield, {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': [1, 2]})

        capsAIM = Body3['capsAIM']
        self.assertEqual(capsAIM, {'fun3dAIM;su2AIM;egadsTessAIM;aflr4AIM;tetgenAIM;aflr3AIM;masstranAIM': {'Body': True, 'Nodes': [], 'Edges': [], 'Faces': []}})

#==============================================================================
    # Build geometry - force regeneration
    def test_build(self):

        self.myGeometry.despmtr.wing.lesweep = 12.0
        self.myGeometry.build()

#==============================================================================
    # View geometry
    def test_viewGeometry(self):

        try:
            self.myGeometry.view(filename = "geomView", showImage =False)

            self.assertTrue(os.path.isfile("geomView_0.png"))
            self.assertTrue(os.path.isfile("geomView_1.png"))
            self.assertTrue(os.path.isfile("geomView_2.png"))

        except ImportError as e: # The machine doesn't have matplotlib
            self.assertTrue(True)

        except pyCAPS.CAPSError as e:
            raise e
        except: # Ignore other errors  TclError being thrown up by matplotlib on some machines
            pass


#==============================================================================
    # Save geometry
    def test_saveGeometry(self):

        self.myGeometry.save("myGeometry.egads")
        self.assertTrue(os.path.isfile("myGeometry.egads"))

#==============================================================================
    # Create html tree
    def test_createTree(self):

        self.myGeometry.createTree(internetAccess = True,
                                   embedJSON = True,
                                   reverseMap = True)
        self.assertTrue(os.path.isfile("myGeometry.html"))

#==============================================================================
    # Test getting the EGADS bodies
    def test_bodies(self):

        bodies, units = self.myGeometry.bodies()

        for key in bodies:
            self.assertIsInstance(bodies[key], egads.ego)
        self.assertEqual(units, pyCAPS.Unit("ft"))

        self.assertEqual(sorted(bodies.keys()), sorted(["Farfield", "Wing1", "Wing2", "aBox"]))

        boundingBox = {}
        for name, body in bodies.items():
            box  = body.getBoundingBox()

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        self.assertAlmostEqual(boundingBox["Farfield"][0], -80, 3)
        self.assertAlmostEqual(boundingBox["Farfield"][3],  80, 3)

#==============================================================================
    # Get parameters from EGADS file
    def test_staticGeometry(self):

        # Make a new problem and load the egads file
        problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=self.file, outLevel=0 ); self.__class__.iProb += 1

        # Make sure geometry is built to dump unitGeom.egads (this should be replaced with saveGeometry)
        problem.geometry.build()

        # Make a new problem and load the egads file
        problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile="unitGeom.egads", outLevel=0 ); self.__class__.iProb += 1

        self.assertEqual(problem.geometry.despmtr.aspect, 5.0)
        self.assertEqual(problem.geometry.cfgpmtr.series, 8412)

        self.assertAlmostEqual(problem.geometry.outpmtr.span, 14.14, 2)
        self.assertAlmostEqual(problem.geometry.outpmtr.cmean, 2.82843, 5)

        # Parameters are now read only because we started from an egads file
        with self.assertRaises(pyCAPS.CAPSError) as e:
            problem.geometry.despmtr.wing.lesweep = 20.0

        self.assertEqual(e.exception.errorName, "CAPS_READONLYERR")

#=============================================================================-
    def test_parameterIO(self):
        problem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        problem.geometry.writeParameters("unitGeom.param")
        problem.geometry.readParameters("unitGeom.param")

        os.remove("unitGeom.param")

        with self.assertRaises(pyCAPS.CAPSError) as e:
            problem.geometry.readParameters("unitGeom.param")
        self.assertEqual(e.exception.errorName, "OCSM_FILE_NOT_FOUND")

#==============================================================================
    # Test value derivatives
    def test_deriv(self):

        problem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        for fac in [1.1, 1.2]:
            problem.geometry.despmtr.area *= fac
        
            aspect = problem.geometry.despmtr.aspect
            area   = problem.geometry.despmtr.area
        
            cmean            = (area/aspect)**0.5
            cmean_areaTrue   = 0.5/cmean * 1/aspect
            cmean_aspectTrue = 0.5/cmean * -area/aspect**2
        
            # Compute derivative w.r.t. area and aspect
            cmean_area   = problem.geometry.outpmtr["cmean"].deriv("area")
            cmean_aspect = problem.geometry.outpmtr["cmean"].deriv("aspect")
        
            self.assertAlmostEqual(cmean_areaTrue  , cmean_area  , 5)
            self.assertAlmostEqual(cmean_aspectTrue, cmean_aspect, 5)

        for fac in [1.1, 1.2]:
            problem.geometry.despmtr.area *= fac
        
            self.assertEqual(problem.geometry.outpmtr["dummyRow"].deriv("taper"), [0.0, 0.0, 0.0])
            self.assertEqual(problem.geometry.outpmtr["dummyRow"].deriv("twist"), [0.0, 0.0, 0.0])

        sphereR = problem.geometry.despmtr.sphereR
        sphereVol             = 4/3 * math.pi *     sphereR**3
        sphereVol_sphereRTrue = 4/3 * math.pi * 3 * sphereR**2

        # Compute derivative w.r.t. sphereR
        sphereVol_sphereR = problem.geometry.outpmtr["sphereVol"].deriv("sphereR")

        #self.assertAlmostEqual(sphereVol    , problem.geometry.outpmtr.sphereVol, 3)
        #self.assertAlmostEqual(0.94, sphereVol_sphereR/sphereVol_sphereRTrue, 2)

        despMat = problem.geometry.despmtr.despMat
        nrow = len(despMat)
        ncol = len(despMat[0])

        boxVolTrue = despMat[0][1] * despMat[1][1] * despMat[2][1]
        boxCGTrue = [despMat[0][0] + despMat[0][1]/2,
                     despMat[1][0] + despMat[1][1]/2,
                     despMat[2][0] + despMat[2][1]/2]
                  
        boxVol_despMatTrue = [0,
                                          1 * despMat[1][1] * despMat[2][1],
                              0,
                              despMat[0][1] *             1 * despMat[2][1],
                              0,
                              despMat[0][1] * despMat[1][1] *             1]

        #  despMat indexes:   00  01  10  11  20  21
        boxCG_despMatTrue = [[ 1, .5,  0,  0,  0,  0],
                             [ 0,  0,  1, .5,  0,  0],
                             [ 0,  0,  0,  0,  1, .5]]

        # Compute derivative w.r.t. despMat
        boxVol_despMat = problem.geometry.outpmtr["boxVol"].deriv("despMat")
        boxCG_despMat  = problem.geometry.outpmtr["boxCG"].deriv("despMat")

        self.assertAlmostEqual(boxVolTrue, problem.geometry.outpmtr.boxVol, 5)

        self.assertAlmostEqual(boxCGTrue[0], problem.geometry.outpmtr.boxCG[0], 5)
        self.assertAlmostEqual(boxCGTrue[1], problem.geometry.outpmtr.boxCG[1], 5)
        self.assertAlmostEqual(boxCGTrue[2], problem.geometry.outpmtr.boxCG[2], 5)

        self.assertAlmostEqual(boxVol_despMatTrue[0], boxVol_despMat[0], 5)
        self.assertAlmostEqual(boxVol_despMatTrue[1], boxVol_despMat[1], 5)
        self.assertAlmostEqual(boxVol_despMatTrue[2], boxVol_despMat[2], 5)
        self.assertAlmostEqual(boxVol_despMatTrue[3], boxVol_despMat[3], 5)
        self.assertAlmostEqual(boxVol_despMatTrue[4], boxVol_despMat[4], 5)
        self.assertAlmostEqual(boxVol_despMatTrue[5], boxVol_despMat[5], 5)

        for irow in range(nrow):
            for icol in range(ncol):
                self.assertAlmostEqual(boxCG_despMatTrue[irow][icol], boxCG_despMat[irow][icol], 5)

        for irow in range(nrow):
            for icol in range(ncol):
                boxVol_despMat = problem.geometry.outpmtr["boxVol"].deriv("despMat[{},{}]".format(irow+1, icol+1))
                self.assertAlmostEqual(boxVol_despMatTrue[ncol*irow + icol], boxVol_despMat, 5)

#==============================================================================
    # Test value derivatives bug in OpenCSM. Remove this whenthe bug is fixed.
    def off_test_deriv2(self):

        problem1 = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile="sup.csm", outLevel=0); self.__class__.iProb += 1

        inlet = problem1.geometry

        inlet.despmtr['dp_r'].value = 3.0
        inlet.despmtr['dp_n'].value = 1.5

        #self.assertAlmostEqual(inlet.outpmtr["volume"].value, 149.23536335509598, 2)

        #self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_r'), 81.59984899433215, 2)
        #self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_n'), 47.63136146266333, 2)

        #self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_n'), 47.63136146266333, 2)
        #self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_r'), 81.59984899433215, 2)


        problem2 = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile="sup.csm", outLevel=0); self.__class__.iProb += 1

        inlet = problem2.geometry

        inlet.despmtr['dp_r'].value = 3.0
        inlet.despmtr['dp_n'].value = 1.5

        self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_n'), 47.63136146266333, 2)
        self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_r'), 82.20855787328674, 2)

        self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_r'), 82.20855787328674, 2)
        self.assertAlmostEqual(inlet.outpmtr["volume"].deriv('dp_n'), 47.63136146266333, 2)

if __name__ == '__main__':
    unittest.main()
