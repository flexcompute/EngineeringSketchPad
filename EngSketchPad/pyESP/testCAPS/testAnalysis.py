import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion
from sys import version_info

import pyCAPS

class TestAnalysis(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.problemName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        m    = pyCAPS.Unit("meter")
        kg   = pyCAPS.Unit("kg")
        s    = pyCAPS.Unit("s")
        K    = pyCAPS.Unit("Kelvin")

        unitSystem={"mass":kg, "length":m, "time":s, "temperature":K}

        cls.myAnalysis = cls.myProblem.analysis.create(aim = "su2AIM",
                                                       name = "su2",
                                                       capsIntent = "CFD",
                                                       unitSystem=unitSystem)

    @classmethod
    def tearDownClass(cls):
        del cls.myAnalysis
        del cls.myProblem
        cls.cleanUp()
        pass

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

        if os.path.isfile("myAnalysisGeometry.egads"):
            os.remove("myAnalysisGeometry.egads")

        if os.path.isfile("su2.html"):
            os.remove("su2.html")

        if os.path.isfile("d3.v3.min.js"):
            os.remove("d3.v3.min.js")

#==============================================================================
    # Same AIM twice
    def test_loadAIM(self):
        self.myProblem.analysis.create(aim = "fun3dAIM", name = "fun3d")

        self.assertEqual(("fun3d" in self.myProblem.analysis.keys() and
                          "su2"   in self.myProblem.analysis.keys()), True)

#==============================================================================
    # Analysis sequence
    def test_analysis(self):

        # Cannot add custom attributes or items
        with self.assertRaises(AttributeError):
            del self.myProblem.analysis
        with self.assertRaises(AttributeError):
            self.myProblem.analysis.foo = 'bar'
        with self.assertRaises(KeyError):
            self.myProblem.analysis["foo"] = 'bar'
        with self.assertRaises(KeyError):
            foo = self.myProblem.analysis["foo"]

        # Cannot add custom attributes
        with self.assertRaises(AttributeError):
            self.myAnalysis.foo = 'bar'

        self.assertIs(self.myAnalysis, self.myProblem.analysis["fun3d"])

        self.assertEqual(sorted(self.myAnalysis.output.keys()), sorted(['CDtot_p', 'CXtot_v', 'CXtot_p',
                                                                        'CDtot_v', 'CLtot_p', 'CLtot_v',
                                                                        'Forces', 'CMXtot', 'CYtot',
                                                                        'CXtot', 'CZtot_v', 'CZtot_p',
                                                                        'CZtot', 'CYtot_p', 'CMZtot_p',
                                                                        'CMZtot_v', 'CMZtot', 'CMYtot',
                                                                        'CYtot_v', 'CLtot', 'CMYtot_v',
                                                                        'CDtot', 'CMYtot_p', 'CMXtot_v', 'CMXtot_p']))

#==============================================================================
    # AIM creates a directory
    def test_analysisDir(self):

        self.assertTrue(os.path.isdir(self.myAnalysis.analysisDir))

#==============================================================================
    # Same name
    def test_sameName(self):

        with self.assertRaises(pyCAPS.CAPSError) as e:
            myAnalysis = self.myProblem.analysis.create(aim = "su2AIM",
                                                        name = "su2")

        self.assertEqual(e.exception.errorName, "CAPS_BADNAME")

#==============================================================================
    # Copy AIM
    def test_copyAIM(self):
        myAnalysis = None
        myAnalysis = self.myProblem.analysis.copy(src = "su2",
                                                  name = "su2_copy")

        self.assertIs(myAnalysis, self.myProblem.analysis["su2_copy"])

        myAnalysis = None
        myAnalysis = self.myProblem.analysis.copy(src = self.myAnalysis,
                                                  name = "su2_copy2")

        self.assertIs(myAnalysis, self.myProblem.analysis["su2_copy2"])

#==============================================================================
    # Multiple intents
    def test_multipleIntents(self):
        myAnalysis = None
        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM",
                                                    name = "TestCFD",
                                                    capsIntent = ["Test", "CFD"])

        self.assertIs(myAnalysis, self.myProblem.analysis["TestCFD"])

#==============================================================================
    # Retrieve geometry attribute
    def test_geomAttrList(self):

        attributeList = self.myAnalysis.geometry.attrList("capsGroup", attrLevel="Body")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

        attributeList = self.myAnalysis.geometry.attrList("capsGroup", attrLevel="Face")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

        attributeList = self.myAnalysis.geometry.attrList("capsGroup", attrLevel="Node")
        self.assertEqual(sorted(attributeList), sorted(['Farfield', 'Wing2', 'Wing1']))

#==============================================================================
    # Retrieve geometry attribute maps
    def test_geomAttrMap(self):

        attributeMap = self.myAnalysis.geometry.attrMap()
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
        attributeMap = self.myAnalysis.geometry.attrMap(reverseMap = True)
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
    def test_input(self):
        ft  = pyCAPS.Unit("ft")
        rad = pyCAPS.Unit("radian")
        deg = pyCAPS.Unit("degree")

        # input value with units
        self.myAnalysis.input.Alpha = 5*deg
        self.assertEqual( 5*deg, self.myAnalysis.input.Alpha)
        self.myAnalysis.input.Alpha += 5*deg
        self.assertEqual(10*deg, self.myAnalysis.input.Alpha)

        self.myAnalysis.input.Beta = 0.174533*rad
        self.assertAlmostEqual(self.myAnalysis.input.Beta/deg, 10.0, 5) # Should be converted to degree

        # check setting incorrect units
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myAnalysis.input.Alpha = 5
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

        # Set/Unset analysis value
        self.assertEqual(self.myAnalysis.input.Mach,None)
        self.myAnalysis.input.Mach = 0.5
        self.assertEqual(self.myAnalysis.input.Mach, 0.5)
        self.myAnalysis.input.Mach = None
        self.assertEqual(self.myAnalysis.input.Mach, None)

        # Modify unitlsss input
        self.myAnalysis.input.Mach = 0.5
        self.assertEqual(0.5 , self.myAnalysis.input.Mach)
        self.myAnalysis.input.Mach += 0.25
        self.assertEqual(0.75, self.myAnalysis.input.Mach)

        # check setting incorrect units
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myAnalysis.input.Mach = 5*ft
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

        # Get analysis input values
        self.assertEqual(self.myAnalysis.input["Overwrite_CFG"].value, True)

        # Get just the names of the input values
        valueList = self.myAnalysis.input.keys()
        self.assertEqual("Proj_Name" in valueList, True)

        # Get a string and add it to another string -> Make sure byte<->str<->unicode is correct between Python 2 and 3
        string = "myTest"
        self.myAnalysis.input.Proj_Name = string

        myString = self.myAnalysis.input.Proj_Name
        self.assertEqual(myString+"_TEST", "myTest_TEST")

#==============================================================================
    # Get analysis out values
    def test_output(self):

        # Get all output values - analysis should be dirty since pre hasn't been run
        with self.assertRaises(pyCAPS.CAPSError) as e:
            valueDict = self.myAnalysis.output.CLtot

        self.assertEqual(e.exception.errorName, "CAPS_DIRTY")

        # Get just the names of the output values
        valueList = self.myAnalysis.output.keys()
        self.assertEqual("CLtot" in valueList, True)


#==============================================================================
    # Analysis - object consistent
    def test_analysis(self):

        self.assertEqual(self.myAnalysis, self.myProblem.analysis["su2"])

#==============================================================================
    # Create html tree
    def test_createTree(self):

        inviscid = {"bcType" : "Inviscid", "wallTemperature" : 1.1}
        self.myAnalysis.input.Boundary_Condition = {"Skin": inviscid,
                                                    "SymmPlane": "SymmetryY",
                                                    "Farfield": "farfield"}

        # Create tree
        self.myAnalysis.createTree(internetAccess = False,
                                   embedJSON = True,
                                   analysisGeom = True,
                                   reverseMap = True)
        self.assertTrue(os.path.isfile("su2.html"))
        self.assertTrue(os.path.isfile("d3.v3.min.js"))

#==============================================================================
    # Get analysis info
    def test_info(self):

        # Make sure it returns a integer
        state = self.myAnalysis.info()
        self.assertEqual(isinstance(state, int), True)

        # Make sure it returns a dictionary
        infoDict = self.myAnalysis.info(printinfo = False,  infoDict = True)
        self.assertEqual(infoDict["analysisDir"], self.myAnalysis.analysisDir)

#==============================================================================
    # Get analysis info with a NULL intent
    def test_intent(self):

        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM"
                                                    ) # Use default capsIntent

        infoDict = myAnalysis.info(printInfo = False,  infoDict = True)
        self.assertEqual(infoDict["intent"], None)

        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM",
                                                    capsIntent = None) # Use no capsIntent

        infoDict = myAnalysis.info(printinfo = False,  infoDict = True)
        self.assertEqual(infoDict["intent"], None)

        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM",
                                                    capsIntent = "") # Use no capsIntent

        infoDict = myAnalysis.info(printinfo = False,  infoDict = True)
        self.assertEqual(infoDict["intent"], None)

        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM",
                                                    capsIntent = "STRUCTURE") # Use different capsIntent

        infoDict = myAnalysis.info(printinfo = False,  infoDict = True)
        self.assertEqual(infoDict["intent"], "STRUCTURE")

        myAnalysis = self.myProblem.analysis.create(aim = "fun3dAIM",
                                                    capsIntent = ["CFD","STRUCTURE"]) # Use multiple capsIntent

        infoDict = myAnalysis.info(printinfo = False,  infoDict = True)
        self.assertEqual(infoDict["intent"], ['CFD', 'STRUCTURE'])

#==============================================================================
    # Adding/getting attributes
    def test_attributes(self):

        # Check list attributes
        self.myAnalysis.attr.create("testAttr", [1, 2, 3])
        self.assertEqual(self.myAnalysis.attr.testAttr, [1,2,3])

        # Check float attributes
        self.myAnalysis.attr.create("testAttr_2", 10.0)
        self.assertEqual(self.myAnalysis.attr.testAttr_2, 10.0)

        # Check string attributes
        self.myAnalysis.attr.create("testAttr_3", "anotherAttribute")
        self.assertEqual(self.myAnalysis.attr.testAttr_3, "anotherAttribute")

        # Check over writing attribute
        self.myAnalysis.attr.create("testAttr_2", 30.0, overwrite=True)
        self.assertEqual(self.myAnalysis.attr.testAttr_2, 30.0)

        self.myAnalysis.attr.create("arrayOfStrings", ["1", "2", "3", "4"])
        self.assertEqual(self.myAnalysis.attr.arrayOfStrings, ["1", "2", "3", "4"])

#==============================================================================
    # Save geometry
    def test_saveGeometry(self):

        self.myAnalysis.geometry.save("myAnalysisGeometry")
        self.assertTrue(os.path.isfile("myAnalysisGeometry.egads"))

#==============================================================================
    # Test getting the EGADS bodies
    def test_bodies(self):

        bodies, units = self.myAnalysis.geometry.bodies()

        self.assertEqual(sorted(bodies.keys()), sorted(["Farfield", "Wing1", "Wing2"]))

        boundingBox = {}
        for name, body in bodies.items():
            box  = body.getBoundingBox()

            boundingBox[name] = [box[0][0], box[0][1], box[0][2],
                                 box[1][0], box[1][1], box[1][2]]

        self.assertAlmostEqual(boundingBox["Farfield"][0], -80, 3)
        self.assertAlmostEqual(boundingBox["Farfield"][3],  80, 3)

#==============================================================================
    # Analysis list of tuples with a mixture of types
    def test_tuple(self):

        problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        nastran = problem.analysis.create(aim="nastranAIM",
                                          name="nastran",
                                          capsIntent = "OML")

        x = {'DESMAX': 200, 'CONV1': 1e-15, 'CONVDV':1e-15, 'CONVPR':1e-15, 'P1':'a', 'P2':'b'}
        nastran.input.Design_Opt_Param = x
        self.assertEqual(nastran.input.Design_Opt_Param, x)

#==============================================================================
    # Circularly linked analysis with auto execution
    def test_circularAutoExec(self):

        self.problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile="airfoilSection.csm", outLevel=0); self.__class__.iProb += 1

        self.xfoil0 = self.problem.analysis.create(aim="xfoilAIM",
                                                   name="xfoil0")

        self.xfoil1 = self.problem.analysis.create(aim="xfoilAIM",
                                                   name="xfoil1")

        self.xfoil2 = self.problem.analysis.create(aim="xfoilAIM",
                                                   name="xfoil2")
        
        self.xfoil1.input["Alpha"].link(self.xfoil0.output["Alpha"])
        
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.xfoil0.input["CL"].link(self.xfoil1.output["CL"])
        self.assertEqual(e.exception.errorName, "CAPS_CIRCULARLINK")

        self.xfoil2.input["Alpha"].link(self.xfoil1.output["Alpha"])

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.xfoil0.input["CL"].link(self.xfoil2.output["CL"])
        self.assertEqual(e.exception.errorName, "CAPS_CIRCULARLINK")
        
        del self.problem
        del self.xfoil0
        del self.xfoil1
        del self.xfoil2

#==============================================================================
    # Analysis exectution and outputs
    def test_analysisExec(self):

        inch = pyCAPS.Unit("inch")
        m    = pyCAPS.Unit("meter")
        kg   = pyCAPS.Unit("kg")

        problem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        masstran = problem.analysis.create(aim="masstranAIM",
                                           name="masstran",
                                           unitSystem={"mass" : kg, "length": inch},
                                           capsIntent = "OML")

        # Min/Max number of points on an edge for quadding
        masstran.input.Edge_Point_Min = 2
        masstran.input.Edge_Point_Max = 10

        madeupium    = {"materialType" : "isotropic",
                        "density"      : 10*kg/m**3}
        unobtainium  = {"materialType" : "isotropic",
                        "density"      : (20, "kilogram/meter**3")}

        masstran.input.Material = {"madeupium" : madeupium, "unobtainium" : unobtainium}

        # Set properties
        shell1 = {"propertyType"      : "Shell",
                  "membraneThickness" : 2.0*inch,
                  "material"          : "madeupium"}

        shell2 = {"propertyType"      : "Shell",
                  "membraneThickness" : (3.0, "mm"),
                  "material"          : "unobtainium"}

        masstran.input.Property = {"Wing1" : shell1, "Wing2" : shell2}

        # get properties (computed automatically)
        Area     = masstran.output.Area
        Mass     = masstran.output.Mass
        Centroid = masstran.output.Centroid
        CG       = masstran.output.CG
        I        = masstran.output.I_Vector
        II       = masstran.output.I_Tensor

        # Check that unit conversion works with the output
        Area.convert("m**2")
        Mass.convert("g")
        Centroid.convert("inch")
        CG.convert("hm")
        I.convert("g*mm**2")

#==============================================================================
    # Analysis coupling
    def off_test_analysisCoupling(self):

        problem = pyCAPS.experiment.Problem(self.problemName, capsFile=self.file, outLevel=0)

        # craete an aims
        masstranWing = problem.analysis.create(aim="masstranAIM",
                                               name="masstranWing",
                                               capsIntent = "Wing")

        masstranFuse = problem.analysis.create(aim="masstranAIM",
                                               name="masstranFuse",
                                               capsIntent = "Fuselage")

        # Material
        madeupium    = {"materialType" : "isotropic",
                        "density"      : 10*lb}
        # Properties
        shell = {"propertyType"      : "Shell",
                 "membraneThickness" : 2.0*mm,
                 "material"          : "madeupium"}

        masstranWing.input.Material = ("madeupium", madeupium)
        masstranWing.input.Property = ("Wing", shell)

        masstranFuse.input.Material = ("madeupium", madeupium)
        masstranFuse.input.Property = ("Fuselage", shell)

        # compute Mass properties
        masstranWing.pre()
        masstranWing.post()

        masstranFuse.pre()
        masstranFuse.post()

        # get mass properties
        wing_mass = masstranWing.Mass
        wing_CG   = masstranWing.CG
        wing_I    = masstranWing.I_Vector

        wing_skin = {"mass":wing_mass, "CG":wing_CG, "massInertia":wing_I}

        fuse_mass = masstranFuse.Mass
        fuse_CG   = masstranFuse.CG
        fuse_I    = masstranFuse.I_Vector

        fuse_skin = {"mass":fuse_mass, "CG":fuse_CG, "massInertia":fuse_I}

        cockpit   = {"mass":3000*lb, "CG":[  8,  0.,  5]*ft, "massInertia":[0.0  , 0.0, 0.0  ]*lb*ft**2}
        Main_gear = {"mass":4500*lb, "CG":[ 86,  0., -4]*ft, "massInertia":[0.5e6, 0.0, 0.5e6]*lb*ft**2}
        Nose_gear = {"mass":1250*lb, "CG":[ 26,  0., -5]*ft, "massInertia":[0.0  , 0.0, 0.0  ]*lb*ft**2}

        avl.input.MassProp = [("wing_skin", wing_skin),
                              ("fuse_skin", fuse_skin),
                              ("cockpit"  , cockpit  ),
                              ("Main_gear", Main_gear),
                              ("Nose_gear", Nose_gear)]

        avl.input["MassProp"].link(masstranWing.output["MassProp"], irow=0)
        avl.input["MassProp"].link(masstranFuse.output["MassProp"], irow=1)

        avl.preAnalysis()
        ####### Run avl ####################
        avl.system("avl caps < avlInput.txt > avlOutput.txt");
        ######################################
        avl.postAnalysis()

        # Get neutral point, CG and MAC
        Xnp = avl.output.Xnp
        Xcg = avl.output.Xcg
        MAC = problem.geometry.despmtr.wing.mac

        StaticMargin = (Xnp - Xcg)/MAC
        print ("--> Xcg = ", Xcg)
        print ("--> Xnp = ", Xnp)
        print ("--> StaticMargin = ", StaticMargin)


if __name__ == '__main__':
    unittest.main()
