import unittest

import os
import glob
import shutil

import pyCAPS

 
class Testxfoil_NACA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_xfoilNACATest"
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Initialize Problem object
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","airfoilSection.csm"),
                                       outLevel = 0)

        # Change a design parameter - area in the geometry
        cls.myProblem.geometry.despmtr.camber = 0.02

        # Load xfoil aim
        cls.xfoil = cls.myProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        cls.xfoil.input.Mach = 0.2
        cls.xfoil.input.Re   = 5.0e5
        cls.xfoil.input.Viscous_Iteration = 100

    @classmethod
    def tearDownClass(cls):
        del cls.xfoil
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_alpha_custom_increment(self):

        AlphaTrue = [0.0, 3.0, 5.0, 9.0, 11.0, 13.0, 14.0, 15.0]
        ClTrue    = [0.2236, 0.6711, 0.8515, 1.0934, 1.1727, 1.2227, 1.23, 1.2291]
        CdTrue    = [0.00637, 0.00822, 0.00957, 0.01718, 0.0238, 0.03607, 0.04548, 0.05711]
        TranXTrue = [0.7239, 0.5064, 0.3542, 0.0542, 0.034, 0.0263, 0.0242, 0.0226]

        # Set custom AoA
        self.xfoil.input.Alpha = AlphaTrue

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i])

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i])

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i])

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i])
            
        # Unset custom AoA
        self.xfoil.input.Alpha = None

#==============================================================================
    def test_alpha_uniform_inrement(self):

        AlphaTrue = [1.0 + x * 0.1 for x in range(0, 11)]
        ClTrue    = [0.3544, 0.3704, 0.3831, 0.3984, 0.4158, 0.4319, 0.4462, 0.4588, 0.476, 0.4931, 0.5091]
        CdTrue    = [0.00675, 0.00684, 0.0069, 0.00698, 0.00705, 0.00713, 0.00722, 0.00729, 0.00738, 0.00744, 0.00754]
        TranXTrue = [0.6539, 0.6467, 0.6395, 0.6321, 0.6247, 0.6171, 0.6101, 0.6029, 0.5957, 0.588, 0.5808]

        # Set AoA seq
        self.xfoil.input.Alpha_Increment = [1.0, 2.0, 0.10]

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i])

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i])

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i])

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i])

        # Unset AoA seq
        self.xfoil.input.Alpha_Increment = None

#==============================================================================
    def test_Cl(self):

        AlphaTrue = -1.308
        ClTrue = 0.1
        CdTrue = 0.00697
        TranXTrue = 0.8055

        # Set custom Cl
        self.xfoil.input.CL = 0.1

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertAlmostEqual(AlphaTrue, Alpha)
        self.assertAlmostEqual(ClTrue, Cl)
        self.assertAlmostEqual(CdTrue, Cd)
        self.assertAlmostEqual(TranXTrue, TranX)
        
        # Unset custom Cl
        self.xfoil.input.CL = None

#==============================================================================
    def test_CL_uniform_increment(self):

        AlphaTrue = [-2.325, 0.244, 1.943, 3.674, 7.233]
        ClTrue    = [0.0 + (x * 0.25) for x in range(0, int(1.0/0.25)+1)]
        CdTrue    = [0.00786, 0.00638, 0.00748, 0.00858, 0.01337]
        TranXTrue = [0.8596, 0.7076, 0.5848, 0.4581, 0.1206]

        # Set Cl seq
        self.xfoil.input.CL_Increment = [0.0, 1.0, .25]

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i])

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i])

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i])

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i])

        # Unset Cl seq
        self.xfoil.input.CL_Increment = None

#==============================================================================
    def test_append(self):

        AlphaTrue = [0.0, 3.0]
        ClTrue    = [0.2236, 0.6711]
        CdTrue    = [0.00637, 0.00822]
        TranXTrue = [0.7239, 0.5064]

        # Set custom AoA
        self.xfoil.input.Alpha = 0.0

        # run xfoil
        self.xfoil.runAnalysis()

        # Append the polar file if it already exists - otherwise the AIM will delete the file
        self.xfoil.input.Append_PolarFile = True

        # Set custom AoA
        self.xfoil.input.Alpha = 3.0

        # run xfoil
        self.xfoil.runAnalysis()

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i])

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i])

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i])

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i])

        # Unset custom AoA
        self.xfoil.input.Alpha = None
        self.xfoil.input.Append_PolarFile = False


class Testxfoil_Kulfan(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_xfoilKulfanTest"
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Load CSM file
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","kulfanSection.csm"),
                                       outLevel = 0)

        # Load xfoil aim
        cls.xfoil = cls.myProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        cls.xfoil.input.Mach = 0.5
        cls.xfoil.input.Re   = 1.0e6
        cls.xfoil.input.Viscous_Iteration = 100

    @classmethod
    def tearDownClass(cls):
        del cls.xfoil
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    def test_alpha_custom_increment(self):

        AlphaTrue = [0.0, 3.0, 5.0, 9.0, 11.0, 12.0]
        ClTrue    = [0.5875, 0.9729, 1.2117, 1.4685, 1.494, 1.4821]
        CdTrue    = [0.00758, 0.00771, 0.00966, 0.02614, 0.04777, 0.06245]
        TranXTrue = [0.577, 0.4657, 0.3656, 0.0727, 0.0458, 0.0412]

        # Set custom AoA
        self.xfoil.input.Alpha = AlphaTrue

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i])

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i])

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i])

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i])
            
        # Unset custom AoA
        self.xfoil.input.Alpha = None

#==============================================================================
    def test_Cl(self):

        AlphaTrue = -3.705
        ClTrue = 0.1
        CdTrue = 0.00908
        TranXTrue = 0.7186

        # Set custom Cl
        self.xfoil.input.CL = 0.1

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        #print("Alpha = ", Alpha)
        #print("Cl = ", Cl)
        #print("Cd = ", Cd)
        #print("Transition location = ", TranX)

        self.assertAlmostEqual(AlphaTrue, Alpha)
        self.assertAlmostEqual(ClTrue, Cl)
        self.assertAlmostEqual(CdTrue, Cd)
        self.assertAlmostEqual(TranXTrue, TranX)
        
        # Unset custom Cl
        self.xfoil.input.CL = None
        
if __name__ == '__main__':
    unittest.main()
