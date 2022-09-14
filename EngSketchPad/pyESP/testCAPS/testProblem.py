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
        cls.problemName = "basicTest"
        cls.iProb = 1

        cls.cleanUp()

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

    @classmethod
    def tearDownClass(cls):
        del cls.myProblem
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove project directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

        dirs = glob.glob( os.path.join("..", cls.problemName + '*') )
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
        myProblem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile= "unitGeom.csm", outLevel=0); self.__class__.iProb += 1
        myProblem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=b"unitGeom.csm", outLevel=0); self.__class__.iProb += 1
        myProblem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile=u"unitGeom.csm", outLevel=0); self.__class__.iProb += 1

#==============================================================================
    def test_close(self):
        myProblem = pyCAPS.Problem(self.problemName + str(self.iProb), capsFile= "unitGeom.csm", outLevel=0); self.__class__.iProb += 1

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

        name = self.problemName + str(self.iProb); self.__class__.iProb += 1

        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)
        self.assertNotEqual(self.myProblem, myProblemNew)
        
        self.assertEqual(myProblemNew.name, name)

        # Test a problemName with a relative path
        name = os.path.join("..", self.problemName + str(self.iProb)); self.__class__.iProb += 1
        myProblemNew = pyCAPS.Problem(name, capsFile=self.file, outLevel=0)
        self.assertNotEqual(self.myProblem, myProblemNew)

        del myProblemNew # close the problem
        shutil.rmtree(name)


        # Create a problem, close it, and re-create with the same name
        name = self.problemName + str(self.iProb); self.__class__.iProb += 1
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

        problem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

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

        problem = pyCAPS.Problem(self.problemName+str(self.iProb), capsFile=self.file, outLevel=0); self.__class__.iProb += 1

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
        
#==============================================================================
    def test_intentPhrase(self):

        problem = pyCAPS.Problem(self.problemName+str(self.iProb), phaseName="Phase0", capsFile=self.file, outLevel=0); self.__class__.iProb += 1

        problem.intentPhrase(["Phase does amazing things!"])

        problem.parameter.create("Mach", 0.5)
        
        hist = problem.parameter["Mach"].history()
        
        self.assertEqual(1, len(hist))
        self.assertEqual(["Phase does amazing things!"], hist[0].intentPhrase)

        problem.intentPhrase(["but I have moved on", "to greener pastures"])

        problem.parameter.Mach = 0.75

        hist = problem.parameter["Mach"].history()

        self.assertEqual(2, len(hist))
        self.assertEqual(["Phase does amazing things!"], hist[0].intentPhrase)
        self.assertEqual(["but I have moved on", "to greener pastures"], hist[1].intentPhrase)

#==============================================================================
    def test_phase(self):

        #--------------------------
        # Initialize Problem object
        # This would be the first phasing script: phase0.py
        #--------------------------
        problemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile="airfoilSection.csm", outLevel=0)

        # Load xfoil aim
        xfoil = myProblem.analysis.create(aim = "xfoilAIM", name = "xfoil")

        # Set new Mach/Alt parameters
        xfoil.input.Mach  = 0.5
        xfoil.input.Alpha = 1.0

        # Retrieve results
        CL = xfoil.output.CL
        CD = xfoil.output.CD

        myProblem.closePhase()

        #--------------------------
        # Initialize Problem from the last phase and make a new phase
        # This would be the second phasing script: phase1.py
        #--------------------------
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0, phaseContinuation = True)

        xfoil = myProblem.analysis["xfoil"]

        # Retrieve results
        self.assertEqual(CL, xfoil.output.CL)
        self.assertEqual(CD, xfoil.output.CD)

        xfoil.input.Alpha = 3.0
        self.assertAlmostEqual(0.4396, xfoil.output.CL, 3)

        myProblem.close()

        #--------------------------
        # Emulate continued work on a phase
        # This would be the second phasing script: phase1.py
        # Add a line at the end and run the script again
        #--------------------------
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0, phaseContinuation = True)
        self.assertTrue(myProblem.journaling())

        xfoil = myProblem.analysis["xfoil"]
        self.assertTrue(myProblem.journaling())

        # Retrieve results as was done in the first phase development
        self.assertEqual(CL, xfoil.output.CL)
        self.assertTrue(myProblem.journaling())
        self.assertEqual(CD, xfoil.output.CD)
        self.assertTrue(myProblem.journaling())

        # Run Alpha=3 as in the first phase development
        xfoil.input.Alpha = 3.0
        self.assertTrue(myProblem.journaling())
        self.assertAlmostEqual(0.4396, xfoil.output.CL, 3)
        self.assertTrue(myProblem.journaling())

        # Add a new step to the phase development
        xfoil.input.Alpha = 2.0
        self.assertAlmostEqual(0.2918, xfoil.output.CL, 3)

        myProblem.close()

        #--------------------------
        # Emulate continued work on a phase
        # This would be the second phasing script: phase1.py
        # Add a line at the beginning, remove some lines and run the script again
        #--------------------------
    
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0, phaseContinuation = False)
        self.assertFalse(myProblem.journaling())

        xfoil = myProblem.analysis["xfoil"]

        # Run Alpha=4 which is new
        xfoil.input.Alpha = 4.0
        self.assertAlmostEqual(0.5897, xfoil.output.CL, 3)

        # Retrieve results as was done in the first phase development
        xfoil.input.Alpha = 1.0
        self.assertEqual(CL, xfoil.output.CL)
        self.assertEqual(CD, xfoil.output.CD)

        # Run Alpha=3 as in the first phase development
        xfoil.input.Alpha = 3.0
        self.assertAlmostEqual(0.4396, xfoil.output.CL, 3)

        # Run Alpha=2 as in the second phase development
        xfoil.input.Alpha = 2.0
        self.assertAlmostEqual(0.2918, xfoil.output.CL, 3)

        myProblem.close()

#==============================================================================
    def test_Scratch_Continuation(self):

        #--------------------------
        # Initialize Problem object
        # This would be the first phasing script: phase0.py
        #--------------------------
        problemName = self.problemName + "_Scratch"
        myProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile="airfoilSection.csm", outLevel=0)

        # Load xfoil aim
        xfoil = myProblem.analysis.create(aim = "xfoilAIM", name = "xfoil")

        # Set new Mach/Alt parameters
        xfoil.input.Mach  = 0.5
        xfoil.input.Alpha = 1.0

        # Retrieve results
        CL = xfoil.output.CL
        CD = xfoil.output.CD

        myProblem.close()

        #--------------------------
        # Use continuation with the Scratch phase
        #--------------------------
        myProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile="airfoilSection.csm", outLevel=0, phaseContinuation = True)
        self.assertTrue(myProblem.journaling())

        # Load xfoil aim
        xfoil = myProblem.analysis.create(aim = "xfoilAIM", name = "xfoil")
        self.assertTrue(myProblem.journaling())

        # Set new Mach/Alt parameters
        xfoil.input.Mach  = 0.5
        self.assertTrue(myProblem.journaling())
        xfoil.input.Alpha = 1.0
        self.assertTrue(myProblem.journaling())

        # Retrieve results
        self.assertEqual(CL, xfoil.output.CL)
        self.assertTrue(myProblem.journaling())
        self.assertEqual(CD, xfoil.output.CD)
        self.assertTrue(myProblem.journaling())
        
        
if __name__ == '__main__':
    unittest.main()
