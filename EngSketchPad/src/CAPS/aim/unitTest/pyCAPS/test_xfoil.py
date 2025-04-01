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
        cls.capsFile = capsFile = os.path.join("..","csmData","airfoilSection.csm")

        cls.cleanUp()

        # Initialize Problem object
        cls.capsProblem = pyCAPS.Problem(problemName = cls.problemName,
                                         capsFile = cls.capsFile,
                                         outLevel = 0)

        # Change a design parameter - camber in the geometry
        cls.capsProblem.geometry.despmtr.camber = 0.02

        # Load xfoil aim
        cls.xfoil = cls.capsProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        cls.xfoil.input.Mach = 0.2
        cls.xfoil.input.Re   = 5.0e5
        cls.xfoil.input.Viscous_Iteration = 100

    @classmethod
    def tearDownClass(cls):
        del cls.xfoil
        del cls.capsProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)


#==============================================================================
    def test_normalize(self):

        # Initialize Problem object
        capsProblem = pyCAPS.Problem(problemName = self.problemName + "_normalize",
                                     capsFile = os.path.join("..","csmData","airfoilSection.csm"),
                                     outLevel = 0)

        # Change a design parameter - camber, scale and rotate
        capsProblem.geometry.despmtr.camber = 0.02
        capsProblem.geometry.despmtr.chord = 0.09
        capsProblem.geometry.despmtr.alpha = 10.00
        capsProblem.geometry.despmtr.offset = [0.1, 0.2, 0.3]

        # Load xfoil aim
        xfoil = capsProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        xfoil.input.Mach = 0.2
        xfoil.input.Re   = 5.0e5
        xfoil.input.Viscous_Iteration = 100

        Alphas =  [0.0, 2.0, 5.0, 7.0, 8.0, 10.0]

        # Compute using non-scaled or rotated
        self.xfoil.input.Alpha = Alphas

        AlphaTrue = self.xfoil.output.Alpha
        ClTrue    = self.xfoil.output.CL
        CdTrue    = self.xfoil.output.CD
        TranXTrue = self.xfoil.output.Transition_Top

        self.xfoil.input.Alpha = None

        # Set custom AoA
        xfoil.input.Alpha = Alphas

        # Retrieve results, the arifoil should be normalized and the same
        Alpha = xfoil.output.Alpha
        Cl    = xfoil.output.CL
        Cd    = xfoil.output.CD
        TranX = xfoil.output.Transition_Top

        # Check
        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 6)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 6)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 6)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 6)

#==============================================================================
    def test_findTE(self):

        # Initialize Problem object
        capsProblem = pyCAPS.Problem(problemName = self.problemName + "_camber",
                                     capsFile = os.path.join("..","csmData","airfoilSection.csm"),
                                     outLevel = 0)

        # Check finding the TE with a very large camber
        capsProblem.geometry.despmtr.camber = 0.6
        capsProblem.geometry.cfgpmtr.sharpTE = 0

        # Load xfoil aim
        xfoil = capsProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number
        xfoil.input.Mach = 0.2

        Alphas =  [0.0, 2.0, 5.0, 7.0, 8.0, 10.0]

        # Set custom AoA
        xfoil.input.Alpha = Alphas

        # Retrieve results
        Cl    = xfoil.output.CL
        Cd_p  = xfoil.output.CD_p
        # print()
        # print("Cl", Cl)
        # print("Cd_p", Cd_p)
        # print()
        ClTrue   = [6.43239337898359, 6.75340744254531, 7.22292325613111, 7.52720472349114, 7.67645735051033, 7.96861490777808]
        Cd_pTrue = [-0.0276050055602635, -0.0242341126499494, -0.0171414876853362, -0.0112335579038935, -0.00803327501699627, -0.00145106077886209]

        # Check
        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 4)

        self.assertEqual(len(Cd_pTrue), len(Cd_p))
        for i in range(len(Cd_pTrue)):
            self.assertAlmostEqual(Cd_pTrue[i], Cd_p[i], 4)

#==============================================================================
    def test_alpha_custom_increment(self):

        AlphaTrue =  [0.0, 2.0, 5.0, 7.0, 8.0, 10.0]
        ClTrue    =  [0.2242, 0.5071, 0.846, 0.9851, 1.0486, 1.1502]
        CdTrue    =  [0.00637, 0.00749, 0.0095, 0.01287, 0.01501, 0.01991]
        TranXTrue =  [0.724, 0.5816, 0.3571, 0.14, 0.0785, 0.0424]

        # Set custom AoA
        self.xfoil.input.Alpha = AlphaTrue

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

        # Unset custom AoA
        self.xfoil.input.Alpha = None

#==============================================================================
    def test_alpha_uniform_inrement(self):

        AlphaTrue =  [1.0 + x * 0.1 for x in range(0, 6)]
        ClTrue    =  [0.3529, 0.3679, 0.3815, 0.3978, 0.4119, 0.4274]
        CdTrue    =  [0.00674, 0.00681, 0.00688, 0.00696, 0.00704, 0.0071]
        TranXTrue =  [0.6544, 0.6473, 0.64, 0.6326, 0.6255, 0.6181]

        # Set AoA seq
        self.xfoil.input.Alpha_Increment = [1.0, 1.5, 0.10]

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

        # Unset AoA seq
        self.xfoil.input.Alpha_Increment = None

#==============================================================================
    def test_Cl(self):

        AlphaTrue =  -1.31
        ClTrue    =  0.1
        CdTrue    =  0.00697
        TranXTrue =  0.8059

        # Set custom Cl
        self.xfoil.input.CL = 0.1

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue, Cl, 3)
        self.assertAlmostEqual(CdTrue, Cd, 3)
        self.assertAlmostEqual(TranXTrue, TranX, 3)

        # Unset custom Cl
        self.xfoil.input.CL = None

#==============================================================================
    def test_CL_uniform_increment(self):

        AlphaTrue =  [-2.322, 0.239, 1.957, 3.799, 7.23]
        ClTrue    =  [0.0 + (x * 0.25) for x in range(0, int(1.0/0.25)+1)]
        CdTrue    =  [0.00786, 0.00638, 0.00746, 0.00858, 0.01336]
        TranXTrue =  [0.8597, 0.708, 0.5847, 0.4523, 0.1209]

        # Set Cl seq
        self.xfoil.input.CL_Increment = [0.0, 1.0, .25]

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

        # Unset Cl seq
        self.xfoil.input.CL_Increment = None

#==============================================================================
    def test_append(self):

        AlphaTrue =  [0.0, 3.0]
        ClTrue    =  [0.2242, 0.6723]
        CdTrue    =  [0.00637, 0.00815]
        TranXTrue =  [0.724, 0.5065]

        # Set custom AoA
        self.xfoil.input.Alpha = 0.0
        self.xfoil.input.Append_PolarFile = False

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
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

        # Unset custom AoA
        self.xfoil.input.Alpha = None
        self.xfoil.input.Append_PolarFile = False

#==============================================================================
    def test_phase(self):

        # Initialize Problem object
        problemName = self.problemName + "_Phase"
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=self.capsFile, outLevel=0)

        # Load xfoil aim
        xfoil = capsProblem.analysis.create(aim = "xfoilAIM", name = "xfoil")

        # Set new Mach/Alt parameters
        xfoil.input.Mach  = 0.5
        xfoil.input.Alpha = 1.0

        # Retrieve results
        CL = xfoil.output.CL
        CD = xfoil.output.CD

        capsProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        capsProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        xfoil = capsProblem.analysis["xfoil"]

        # Retrieve results
        self.assertEqual(CL, xfoil.output.CL)
        self.assertEqual(CD, xfoil.output.CD)

        xfoil.input.Alpha = 3.0
        self.assertAlmostEqual(0.4396, xfoil.output.CL, 3)

class Testxfoil_Kulfan(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_xfoilKulfanTest"
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Load CSM file
        cls.capsProblem = pyCAPS.Problem(problemName = cls.problemName,
                                         capsFile = os.path.join("..","csmData","kulfanSection.csm"),
                                         outLevel = 0)

        # Load xfoil aim
        cls.xfoil = cls.capsProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        cls.xfoil.input.Mach = 0.5
        cls.xfoil.input.Re   = 1.0e6
        cls.xfoil.input.Viscous_Iteration = 100

    @classmethod
    def tearDownClass(cls):
        del cls.xfoil
        del cls.capsProblem
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

        AlphaTrue =  [0.0, 1.0, 1.75, 2.5, 3.0, 4.0]
        ClTrue    =  [0.5878, 0.7116, 0.8169, 0.9101, 0.9724, 1.0947]
        CdTrue    =  [0.00758, 0.00693, 0.00693, 0.00739, 0.00772, 0.00855]
        TranXTrue =  [0.5771, 0.5414, 0.5134, 0.4854, 0.466, 0.4218]

        # Set custom AoA
        self.xfoil.input.Alpha = AlphaTrue

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

        # Unset custom AoA
        self.xfoil.input.Alpha = None

#==============================================================================
    def test_Cl(self):

        AlphaTrue =  -3.7
        ClTrue    =  0.1
        CdTrue    =  0.00889
        TranXTrue =  0.7186

        # Set custom Cl
        self.xfoil.input.CL = 0.1

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue, Cl, 3)
        self.assertAlmostEqual(CdTrue, Cd, 3)
        self.assertAlmostEqual(TranXTrue, TranX, 3)

        # Unset custom Cl
        self.xfoil.input.CL = None


class Testxfoil_Parsec(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_xfoilParsecTest"

        cls.cleanUp()

        # Load CSM file
        cls.capsProblem = pyCAPS.Problem(problemName = cls.problemName,
                                         capsFile = os.path.join("..","csmData","parsecSection.csm"),
                                         outLevel = 0)

        # Load xfoil aim
        cls.xfoil = cls.capsProblem.analysis.create(aim = "xfoilAIM")

        # Set Mach number, Reynolds number
        cls.xfoil.input.Mach = 0.5
        cls.xfoil.input.Re   = 1.0e6
        cls.xfoil.input.Viscous_Iteration = 100

    @classmethod
    def tearDownClass(cls):
        del cls.xfoil
        del cls.capsProblem
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

        AlphaTrue =  [0.0, 1.0, 1.75, 2.5, 3.0, 4.0]
        ClTrue    =  [0.0, 0.116551006161837, 0.202697964696352, 0.287150860648956, 0.343494900507778, 0.472305881135068]
        CdTrue    =  [0.00620798257714116, 0.00631820696421362, 0.00654607156464654, 0.00693180688295738, 0.00727050781339858, 0.00829599074999285]
        TranXTrue =  [0.593028893427919, 0.496005232275858, 0.412049671019148, 0.320575792075857, 0.261365920500577, 0.160530241137239]

        # Set custom AoA
        self.xfoil.input.Alpha = AlphaTrue
        self.xfoil.input.CL = None

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertEqual(len(AlphaTrue), len(Alpha))
        for i in range(len(AlphaTrue)):
            self.assertAlmostEqual(AlphaTrue[i], Alpha[i], 3)

        self.assertEqual(len(ClTrue), len(Cl))
        for i in range(len(ClTrue)):
            self.assertAlmostEqual(ClTrue[i], Cl[i], 3)

        self.assertEqual(len(CdTrue), len(Cd))
        for i in range(len(CdTrue)):
            self.assertAlmostEqual(CdTrue[i], Cd[i], 3)

        self.assertEqual(len(TranXTrue), len(TranX))
        for i in range(len(TranXTrue)):
            self.assertAlmostEqual(TranXTrue[i], TranX[i], 3)

#==============================================================================
    def test_Cl(self):

        AlphaTrue =  0.85640
        ClTrue    =  0.1
        CdTrue    =  0.0062884
        TranXTrue =  0.510978571240592

        # Set custom Cl
        self.xfoil.input.CL = 0.1
        self.xfoil.input.Alpha = None

        # Retrieve results
        Alpha = self.xfoil.output.Alpha
        Cl    = self.xfoil.output.CL
        Cd    = self.xfoil.output.CD
        TranX = self.xfoil.output.Transition_Top
        # print()
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("TranXTrue = ", TranX)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue, Cl, 3)
        self.assertAlmostEqual(CdTrue, Cd, 3)
        self.assertAlmostEqual(TranXTrue, TranX, 3)

if __name__ == '__main__':
    unittest.main()
