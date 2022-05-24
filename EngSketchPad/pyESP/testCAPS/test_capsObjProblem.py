import unittest

import os
import glob
import shutil
import __main__

from pyCAPS import caps
from pyEGADS import egads

class TestProblem(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.projectName = "testProblem"
        cls.analysisDir1 = "UnitTest1"
        cls.analysisDir2 = "UnitTest2"
        cls.iProb = 1

        cls.cleanUp()

    @classmethod
    def tearDownClass(cls):
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.projectName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove default projectName
        base = os.path.basename(__main__.__file__)
        projectName = os.path.splitext(base)[0]
        if os.path.isdir(projectName):
            shutil.rmtree(projectName)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

        if os.path.isfile("unitGeom2.egads"):
            os.remove("unitGeom2.egads")

#=============================================================================-
    def test_openSingleParametric(self):
        problem = caps.open(self.projectName, None, caps.oFlag.oFileName, self.file, 0)

        #check the pointer equality operator
        self.assertEqual(problem, problem)

        name, otype, stype, link, parent, last = problem.info()

        self.assertEqual(self.projectName     , name)
        self.assertEqual(caps.oType.PROBLEM   , otype)
        self.assertEqual(caps.sType.PARAMETRIC, stype)
        self.assertIsNone(link)
        self.assertIsNone(parent)
        self.assertIsInstance(last, caps.capsOwn)

#=============================================================================-
    # Multiple problems
    def test_openMultiStaic(self):

        problem1 = caps.open(self.projectName+"_1", None, caps.oFlag.oFileName, "unitGeom.csm", 0)

        problem1.writeGeometry("unitGeom2.egads")

        self.assertEqual(4, problem1.size(caps.oType.BODIES, 0))

        problem2 = caps.open(self.projectName+"_2", None, caps.oFlag.oFileName, "unitGeom2.egads", 0)

        self.assertNotEqual(problem1, problem2)

        name, otype, stype, link, parent, last = problem2.info()

        self.assertEqual(self.projectName+"_2", name)
        self.assertEqual(caps.oType.PROBLEM  , otype)
        self.assertEqual(caps.sType.STATIC   , stype)
        self.assertIsNone(link)
        self.assertIsNone(parent)
        self.assertIsInstance(last, caps.capsOwn)

        body, units = problem1.bodyByIndex(1)
        self.assertIsInstance(body, egads.ego)
        body, units = problem1.bodyByIndex(2)
        self.assertIsInstance(body, egads.ego)
        body, units = problem1.bodyByIndex(3)
        self.assertIsInstance(body, egads.ego)

        with self.assertRaises(caps.CAPSError) as e:
            problem1.bodyByIndex(0)
        self.assertEqual(e.exception.errorName, "CAPS_RANGEERR")
        with self.assertRaises(caps.CAPSError) as e:
            problem1.bodyByIndex(5)
        self.assertEqual(e.exception.errorName, "CAPS_RANGEERR")

#=============================================================================-
    def test_childByIndex(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        nGeomIn = problem.size(caps.oType.VALUE, caps.sType.GEOMETRYIN)
        self.assertGreater(nGeomIn, 0)
        for i in range(nGeomIn):
            obj = problem.childByIndex(caps.oType.VALUE, caps.sType.GEOMETRYIN, i+1)

#=============================================================================-
    def test_geometryInOut(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        # get the aspaect ratio and area
        aspectObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYIN, "aspect")
        aspect = aspectObj.getValue()
        self.assertEqual(5.0, aspect)

        areaObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYIN, "area")
        area = areaObj.getValue()
        self.assertEqual(40.0, area)

        cmean = (area/aspect)**0.5
        cmeanObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "cmean")
        data = cmeanObj.getValue()
        self.assertAlmostEqual(cmean, data, 5)

        span = cmean*aspect
        spanObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "span")
        data = spanObj.getValue()
        self.assertAlmostEqual(span, data, 5)

        # modify the aspect ratio and check cmean and span
        aspect *= 1.5
        area *= 1.1
        aspectObj.setValue(aspect)
        areaObj.setValue(area)

        cmean = (area/aspect)**0.5
        data = cmeanObj.getValue()
        self.assertAlmostEqual(cmean, data, 5)

        span = cmean*aspect
        data = spanObj.getValue()
        self.assertAlmostEqual(span, data, 5)

        # test array and matrix
        dummyRowObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "dummyRow")
        data = dummyRowObj.getValue()
        self.assertEqual(data, [1,2,3])

        dummyRow2Obj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "dummyRow2")
        data = dummyRow2Obj.getValue()
        self.assertEqual(data, [1,2,None])

        dummyRow3Obj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "dummyRow3")
        data = dummyRow3Obj.getValue()
        self.assertEqual(data, [2,3,1,1,1])

        dummyColObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "dummyCol")
        data = dummyColObj.getValue()
        self.assertEqual(data, [3,2,1])

        dummyMatObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "dummyMat")
        data = dummyMatObj.getValue()
        self.assertEqual(data, [[1,2],[3,4], [5,6]])

#=============================================================================-
    def test_parameterIO(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        problem.writeParameters("unitGeom.param")
        problem.readParameters("unitGeom.param")

        os.remove("unitGeom.param")

        with self.assertRaises(caps.CAPSError) as e:
            problem.readParameters("unitGeom.param")
        self.assertEqual(e.exception.errorName, "OCSM_FILE_NOT_FOUND")

#=============================================================================-
    def test_queryAnalysis(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1
        nIn, nOut, execute = problem.queryAnalysis("fun3dAIM")

#=============================================================================-
    def test_makeAnalysis(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        self.assertEqual(0, problem.size(caps.oType.ANALYSIS, 0))

        analysis1 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir1)
        analysis2 = problem.makeAnalysis("fun3dAIM", name = self.analysisDir2)

        self.assertEqual(2, problem.size(caps.oType.ANALYSIS, 0))

#=============================================================================-
    def test_makeValue(self):
        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        data = 3
        value = problem.makeValue("scalar_Int", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = [10]*3
        value = problem.makeValue("vector_Int", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = [[10]*3]*2
        value = problem.makeValue("array2D", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = 3.0
        value = problem.makeValue("scalar_Double", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = [10.0]*3
        value = problem.makeValue("vector_Double", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = "foo"
        value = problem.makeValue("string", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        # vector of strings are not supported...
        #data = ["a", "b", "c"]
        #value = problem.makeValue("vector_String", caps.sType.USER, data, units=None)
        #data_out = value.getValue()
        #self.assertEqual(data, data_out)

        data = ["Hello", "World"]
        value = problem.makeValue("list", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = [["Hello", "World"], ["World", "Hello"]]
        value = problem.makeValue("listlist", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

        data = [10]*3
        value = problem.makeValue("param", caps.sType.PARAMETER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)

#=============================================================================-
    # Set verbosity
    def test_setOutLevel(self):

        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        problem.setOutLevel(0)
        problem.setOutLevel(1)
        problem.setOutLevel(2)

        with self.assertRaises(caps.CAPSError) as e:
            problem.setOutLevel(-1)
        self.assertEqual(e.exception.errorName, "CAPS_RANGEERR")

        with self.assertRaises(caps.CAPSError) as e:
            problem.setOutLevel(3)
        self.assertEqual(e.exception.errorName, "CAPS_RANGEERR")

#=============================================================================-
    # Adding/getting attributes
    def test_attributes(self):

        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        # Check list attribute
        data = [1, 2, 3]
        value = problem.makeValue("vector_Int", caps.sType.USER, data)

        problem.setAttr(value)
        del value
        self.assertEqual(problem.attrByName("vector_Int").getValue(), data)
        self.assertEqual(problem.attrByIndex(1).getValue(), data)

        # Check float attribute
        data = 10.0
        value = problem.makeValue("float", caps.sType.USER, data)

        problem.setAttr(value)
        del value
        self.assertEqual(problem.attrByName("float").getValue(), data)
        self.assertEqual(problem.attrByIndex(2).getValue(), data)

        # Check string attributes
        data = "Hello World!"
        value = problem.makeValue("string", caps.sType.USER, data)

        problem.setAttr(value)
        del value
        self.assertEqual(problem.attrByName("string").getValue(), data)
        self.assertEqual(problem.attrByIndex(3).getValue(), data)

        # Check over writing attribute
        data = 20.0
        value = problem.makeValue("string", caps.sType.USER, data)

        problem.setAttr(value)
        del value
        self.assertEqual(problem.attrByName("string").getValue(), data)
        self.assertEqual(problem.attrByIndex(3).getValue(), data)

        self.assertEqual(problem.size(caps.oType.ATTRIBUTES, 0), 3)

        problem.deleteAttr("string")

        self.assertEqual(problem.size(caps.oType.ATTRIBUTES, 0), 2)

        with self.assertRaises(caps.CAPSError) as e:
            problem.attrByName("string")
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

#=============================================================================-
    # Get geometry out values - if geometry has not been set
    def test_getGeometryOut_NotSet(self):

        problem = caps.open(self.projectName+str(self.iProb), None, caps.oFlag.oFileName, self.file, 0); self.__class__.iProb += 1

        # Not set with OutParam
        with self.assertRaises(caps.CAPSError) as e:
            xtipObj = problem.childByName(caps.oType.VALUE, caps.sType.GEOMETRYOUT, "xtip")

        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

if __name__ == '__main__':
    unittest.main()
