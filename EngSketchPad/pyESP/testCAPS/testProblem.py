import unittest

import os
import glob
import shutil
import __main__

from sys import version_info as pyVersion

import pyCAPS

class TestProblem(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = "unitGeom.csm"
        cls.probemName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.myProblem = pyCAPS.Problem(cls.probemName, capsFile=cls.file, outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove project directories
        dirs = glob.glob( cls.probemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        dirs = glob.glob( os.path.join("..", cls.probemName + '*') )
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        # Remove created files
        if os.path.isfile("unitGeom.egads"):
            os.remove("unitGeom.egads")

        if os.path.isfile("myProblem.html"):
            os.remove("myProblem.html")

#==============================================================================
    def test_init(self):
        # Tets all variants of capsFile strings
        myProblem = pyCAPS.Problem(self.probemName + str(self.iProb), capsFile= "unitGeom.csm", outLevel=0); self.__class__.iProb += 1
        myProblem = pyCAPS.Problem(self.probemName + str(self.iProb), capsFile=b"unitGeom.csm", outLevel=0); self.__class__.iProb += 1
        myProblem = pyCAPS.Problem(self.probemName + str(self.iProb), capsFile=u"unitGeom.csm", outLevel=0); self.__class__.iProb += 1

#==============================================================================
    def test_close(self):
        myProblem = pyCAPS.Problem(self.probemName + str(self.iProb), capsFile= "unitGeom.csm", outLevel=0); self.__class__.iProb += 1

        myProblem.close()
        
        # Closing a second time raises and error
        with self.assertRaises(pyCAPS.CAPSError):
            myProblem.close()

#==============================================================================
    # Cannot add Python attributes
    def test_static(self):

        # cannot set non-existing attributes on Problem Object
        with self.assertRaises(AttributeError):
            self.myProblem.foo = "bar"

        # cannot del attributes on Problem Object
        with self.assertRaises(AttributeError):
            del self.myProblem.geometry

#==============================================================================
    # Multiple problems
    def test_multiProblem(self):

        name = self.probemName + str(self.iProb); self.__class__.iProb += 1

        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)
        self.assertNotEqual(self.myProblem, myProblemNew)
        
        self.assertEqual(myProblemNew.name, name)

        # Test a problemName with a relative path
        name = os.path.join("..", self.probemName + str(self.iProb)); self.__class__.iProb += 1
        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)
        self.assertNotEqual(self.myProblem, myProblemNew)

        del myProblemNew # close the problem
        shutil.rmtree(name)


        # Create a problem, close it, and re-create with the same name
        name = self.probemName + str(self.iProb); self.__class__.iProb += 1
        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)

        del myProblemNew # close the problem

        # create the problem again with the same name
        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)


#==============================================================================
    # Set outLevel
    def test_setOutLevel(self):

        self.myProblem.setOutLevel("minimal")
        self.myProblem.setOutLevel("standard")
        self.myProblem.setOutLevel("debug")

        self.myProblem.setOutLevel(2)
        self.myProblem.setOutLevel(1)
        self.myProblem.setOutLevel(0)

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.myProblem.setOutLevel(10)

        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

#==============================================================================
    # Create html tree
    def test_createTree(self):

        self.myProblem.createTree()

        self.assertTrue(os.path.isfile("myProblem.html"))

#==============================================================================
    # Add/Get attributes
    def test_attributes(self):

        ft = pyCAPS.Unit("ft")

        problem = pyCAPS.Problem(self.probemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        # Check list attributes
        problem.attr.create("testAttr", [1, 2, 3])
        self.assertEqual(problem.attr.testAttr, [1,2,3])
        self.assertEqual(problem.attr["testAttr"].value, [1,2,3])

        # Check float attributes
        problem.attr.create("testAttr_2", 10.0*ft)
        self.assertEqual(problem.attr.testAttr_2, 10.0*ft)
        self.assertEqual(problem.attr["testAttr_2"].value, 10.0*ft)

        # Check string attributes
        problem.attr.create("testAttr_3", ["str0", "str1"])
        self.assertEqual(problem.attr.testAttr_3, ["str0", "str1"])
        self.assertEqual(problem.attr["testAttr_3"].value, ["str0", "str1"])

        # Check over writing attribute
        problem.attr.create("testAttr_2", 30.0, overwrite=True)
        self.assertEqual(problem.attr.testAttr_2, 30.0)
        self.assertEqual(problem.attr["testAttr_2"].value, 30.0)

        # Check modifying attribute
        attrObj = problem.attr["testAttr_2"]
        problem.attr.testAttr_2 = 20.0
        self.assertEqual(attrObj.value, 20.0)
        self.assertEqual(problem.attr.testAttr_2, 20.0)

        problem.attr.testAttr_3 = ["CAPS", "attribute"]
        self.assertEqual(problem.attr.testAttr_3, ["CAPS", "attribute"])

        # try to set incorrect units on attribute
        with self.assertRaises(pyCAPS.CAPSError):
            problem.attr.testAttr_2 = 5*ft

        # Delete an attribute
        del problem.attr.testAttr_2
        del problem.attr["testAttr_3"]

        # try to set non-existing attribute
        with self.assertRaises(AttributeError):
            problem.attr.testAttr_2 = 5
        with self.assertRaises(AttributeError):
            problem.attr.testAttr_3 = 6

        # try to get non-existing attribute
        with self.assertRaises(AttributeError):
            val = problem.attr.testAttr_2
        with self.assertRaises(KeyError):
            val = problem.attr["testAttr_2"]

        # try to delete non-existing attribute
        with self.assertRaises(AttributeError):
            del problem.attr.testAttr_2
        with self.assertRaises(KeyError):
            del problem.attr["testAttr_2"]
            
        # Explicitly close the problem object to check for memory leaks
        problem.close()

#==============================================================================
    # Create parameters
    def test_parameters(self):

        deg = pyCAPS.Unit("degree")

        problem = pyCAPS.Problem(self.probemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        problem.parameter.create("Mach", 0.5)

        self.assertEqual(0.5, problem.parameter.Mach)

        avl = problem.analysis.create(aim="avlAIM", name="avl")
        problem.analysis.create(aim="xfoilAIM", name="xfoil")
        xfoil = problem.analysis["xfoil"]

        avl.input.Mach = 0.2
        self.assertEqual(0.2, avl.input.Mach)
        xfoil.input.Mach = 0.3
        self.assertEqual(0.3, xfoil.input.Mach)

        avl.input["Mach"].link(problem.parameter["Mach"])
        xfoil.input["Mach"].link(problem.parameter["Mach"])

        self.assertEqual(0.5, avl.input.Mach)
        self.assertEqual(0.5, xfoil.input.Mach)

        problem.parameter.Mach = 0.4
        self.assertEqual(0.4, avl.input.Mach)
        self.assertEqual(0.4, xfoil.input.Mach)

        xfoil.input["Mach"].unlink()
        problem.parameter.Mach = 0.6
        self.assertEqual(0.6, avl.input.Mach)
        # returns to the originally stored value
        self.assertEqual(0.3, xfoil.input.Mach)
        

if __name__ == '__main__':
    unittest.main()
