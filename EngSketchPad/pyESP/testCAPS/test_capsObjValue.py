import unittest

import os
import glob
import shutil
import __main__

from pyCAPS import caps

class TestValue(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.analysisDir = "UnitTest"
        cls.projectName = "testProblem"
        cls.cleanUp()

        cls.problemObj = caps.open(cls.projectName, None, caps.oFlag.oFileName, cls.file, 0)

    @classmethod
    def tearDownClass(cls):
        del cls.problemObj
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

#=============================================================================-
    def test_makeValue(self):

        ft = caps.Unit("ft")

        data = [0.0, 30000.0, 60000.0]*ft
        value = self.problemObj.makeValue("Altitude1", caps.sType.USER, data)
        data_out = value.getValue()
        self.assertEqual(data, data_out)
        self.assertIsInstance(data_out, caps.Quantity)
        self.assertEqual(caps.Unit("ft"), data_out._units)

#=============================================================================-
    # Set/Get values
    def test_setgetValue(self):

        ft = caps.Unit("ft")

        data = [0.0, 30000.0, 60000.0]*ft
        value = self.problemObj.makeValue("Altitude1", caps.sType.USER, data)

        data = [0.0, 40000.0, 50000.0]*ft
        value.setValue(data)
        data_out = value.getValue()
        self.assertIsInstance(data_out, caps.Quantity)
        for i in range(len(data)):
            self.assertAlmostEqual(data[i], data_out[i], 5)
        self.assertEqual(caps.Unit("ft"), data_out._units)

#=============================================================================-
    # Set limits
    def test_limits(self):

        ft = caps.Unit("ft")

        value = self.problemObj.makeValue("Altitude2", caps.sType.USER, [0.0, 30000.0, 60000.0]*ft)

        limits = [0.0, 80000.0]*ft
        value.setLimits(limits)
        limits_out = value.getLimits()
        self.assertEqual(limits, limits_out)

#=============================================================================-
    # Convert units
    def test_convertUnits(self):

        ft = caps.Unit("ft")

        value = self.problemObj.makeValue("Altitude3", caps.sType.USER, [0.0, 30000.0, 60000.0]*ft)

        alt_m = value.convertValue(60000.0, "m")
        self.assertAlmostEqual(196850.4, alt_m, 1)

#=============================================================================-
#     # Array of strings
#     def test_arrayString(self):
#
#         myValue = self.myProblem.createValue("Strings", ["test1", "test2", "test3"])
#         self.assertEqual(myValue.value, ["test1", "test2", "test3"])

#=============================================================================-
    def test_changefixedShapeAndLength(self):

        ft = caps.Unit("ft")
        m  = caps.Unit("m")

        value = self.problemObj.makeValue("Altitude4", caps.sType.USER, [0.0, 30000.0, 60000.0]*ft)

        # Change the shape/length of fixed value
        with self.assertRaises(caps.CAPSError) as e:
            value.setValue(0.0*ft)

        self.assertEqual(e.exception.errorName, "CAPS_SHAPEERR")

        e = None
        # Change the length of fixed shape value
        with self.assertRaises(caps.CAPSError) as e:
            value.setValue([0.0, 1.0]*ft)

        self.assertEqual(e.exception.errorName, "CAPS_SHAPEERR")

#=============================================================================-
    # Change the shape shape value
    def test_changefixedShape2(self):

        ft = caps.Unit("ft")
        m  = caps.Unit("m")

        value = self.problemObj.makeValue("FixedShape2", caps.sType.USER, 0.0*ft)
        value.setValueProps(dim=caps.vDim.Vector, lfix=caps.Fixed.Change, sfix=caps.Fixed.Fixed, ntype=caps.Null.NotNull)

        value.setValue( [0.0, 1.0]*m )

        with self.assertRaises(caps.CAPSError) as e:
            value.setValue( [[0.0, 1.0], [2.0, 3.0]]*ft )

        self.assertEqual(e.exception.errorName, "CAPS_SHAPEERR")

        value.setValueProps(dim=caps.vDim.Vector, lfix=caps.Fixed.Change, sfix=caps.Fixed.Change, ntype=caps.Null.NotNull)
        value.setValue( [[0.0, 1.0], [2.0, 3.0]]*ft )

if __name__ == '__main__':
    unittest.main()
