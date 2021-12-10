import unittest

import os
import glob
import shutil
import __main__

import pyCAPS

class TestValue(unittest.TestCase):

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


        ft = pyCAPS.Unit("ft")

        cls.myParam = cls.myProblem.parameter.create("Altitude", [0.0, 30000.0, 60000.0]*ft, limits=[0.0, 70000]*ft)

    @classmethod
    def tearDownClass(cls):
        del cls.myAnalysis
        del cls.myParam
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove default problemName
        base = os.path.basename(__main__.__file__)
        problemName = os.path.splitext(base)[0]
        if os.path.isdir(problemName):
            shutil.rmtree(problemName)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

#==============================================================================
    # Value - object consistent
    def test_parameter(self):
        self.assertEqual(self.myParam, self.myProblem.parameter["Altitude"])

        ft = pyCAPS.Unit("ft")

        self.assertEqual(self.myParam.name, "Altitude")
        self.assertEqual(self.myParam.value, [0.0, 30000.0, 60000.0]*ft)
        self.assertEqual(self.myParam.limits, [0.0, 70000]*ft)

#==============================================================================
    # Set values
    def test_value(self):

        ft = pyCAPS.Unit("ft")

        myParam = self.myProblem.parameter.create("Altitude2", [0.0, 30000.0, 60000.0]*ft, limits=[0.0, 70000]*ft)

        self.assertEqual(myParam.value, [0.0, 30000.0, 60000]*ft)

        myParam.value = [0.0, 40000.0, 50000]*ft
        self.assertEqual(myParam.value, [0.0, 40000.0, 50000]*ft)

        myParam.transferValue(pyCAPS.tMethod.Copy, self.myParam)
        self.assertEqual(self.myParam.value, [0.0, 30000.0, 60000.0]*ft)

#==============================================================================
    # Set limits
    def test_limits(self):

        ft = pyCAPS.Unit("ft")

        myParam = self.myProblem.parameter.create("Altitude3", [0.0, 30000.0, 60000.0]*ft, limits=[0.0, 70000]*ft)

        self.assertEqual(myParam.limits, [0.0, 70000.0]*ft)

        myParam.limits = [0.0, 80000.0]*ft
        self.assertEqual(myParam.limits, [0.0, 80000.0]*ft)

#==============================================================================
    # Array of strings
    def test_arrayString(self):

        myParam = self.myProblem.parameter.create("Strings", ["test1", "test2", "test3"])
        self.assertEqual(myParam.value, ["test1", "test2", "test3"])

        myParam2 = self.myProblem.parameter.create("Strings2", ["foo1", "foo2", "foo3"])
        myParam2.transferValue(pyCAPS.tMethod.Copy, myParam)
        self.assertEqual(myParam2.value, ["test1", "test2", "test3"])

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myParam.transferValue(pyCAPS.tMethod.Copy, myParam)
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

#==============================================================================
    # Change the shape/length of fixed value
    def test_changefixedShapeAndLength(self):

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myParam.value = 0.0
        self.assertEqual(e.exception.errorName, "CAPS_SHAPEERR")

#==============================================================================
    # Change the length of fixed shape value
    def test_changefixedShape(self):

        myParam = self.myProblem.parameter.create("FixedShape", [0.0, 30000.0, 60000.0], fixedLength=False)

        myParam.value = [0.0, 1.0]
        self.assertEqual(myParam.value, [0.0, 1.0])

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myParam.transferValue(pyCAPS.tMethod.Copy, myParam)
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

#==============================================================================
    # Change the shape shape value
    def test_changefixedShape2(self):

        myParam = self.myProblem.parameter.create("FixedShape2", [0.0, 30000.0, 60000.0], fixedLength=False)
        with self.assertRaises(pyCAPS.CAPSError) as e:
            myParam.value = [[0.0, 1.0], [2.0, 3.0]]

        self.assertEqual(e.exception.errorName, "CAPS_SHAPEERR")

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myParam.transferValue(pyCAPS.tMethod.Copy, myParam)
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

#==============================================================================
    # Change the length and shape
    def test_changeLengthAndShape(self):

        myParam = self.myProblem.parameter.create("FreeShape", 0.0, fixedLength = False, fixedShape=False)

        myParam.value = [0.0, 1.0]
        self.assertEqual(myParam.value, [0.0, 1.0])

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myParam.transferValue(pyCAPS.tMethod.Copy, myParam)
        self.assertEqual(e.exception.errorName, "CAPS_UNITERR")

#==============================================================================
    # Autolink
    def test_autolink(self):
        
        deg = pyCAPS.Unit("degree")

        self.myAnalysis.input.Beta = 0.0*deg
        self.myAnalysis.input.Alpha = 1.0*deg

        myParam = self.myProblem.parameter.create("Beta", 1.0*deg)
        self.myProblem.autoLinkParameter("Beta")

        myParam = self.myProblem.parameter.create("Alpha", 10.0*deg)
        self.myProblem.autoLinkParameter(myParam)

        self.assertEqual(self.myAnalysis.input.Beta, 1.0*deg)
        self.assertEqual(self.myAnalysis.input.Alpha, 10.0*deg)
        
        # Modify the parameter
        self.myProblem.parameter.Beta = 3.0*deg
        self.myProblem.parameter.Alpha = 5.0*deg
        
        # Check the analysis input is modified
        self.assertEqual(self.myAnalysis.input.Beta, 3.0*deg)
        self.assertEqual(self.myAnalysis.input.Alpha, 5.0*deg)

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myAnalysis.input.Beta = 1.0*deg
            
        self.myAnalysis.input["Beta"].unlink()
        self.myAnalysis.input["Alpha"].unlink()

        self.assertEqual(self.myAnalysis.input.Beta, 0.0*deg)
        self.assertEqual(self.myAnalysis.input.Alpha, 1.0*deg)

        # Auto link all parameters
        self.myProblem.autoLinkParameter()

        self.assertEqual(self.myAnalysis.input.Beta, 3.0*deg)
        self.assertEqual(self.myAnalysis.input.Alpha, 5.0*deg)

#==============================================================================
    def test_matrixValue(self):

        param = self.myProblem.parameter.create("matrixValue", [[1,2], [3,4], [5,6]])
        self.assertEqual(param.value, [[1,2], [3,4], [5,6]])

#==============================================================================
    def test_rowValue(self):

        param = self.myProblem.parameter.create("rowValue", [1,2,3])
        self.assertEqual(param.value, [1,2,3])

#==============================================================================
    def test_colValue(self):

        # Should get convert back to a row
        param = self.myProblem.parameter.create("colValue", [[1,2,3]])
        self.assertEqual(param.value, [1,2,3])

#==============================================================================
    def test_inconsistentValueShape(self):

        with self.assertRaises(pyCAPS.CAPSError):
            value = self.myProblem.parameter.create("inconValue", [[1,2,3], [0]])


if __name__ == '__main__':
    unittest.main()
