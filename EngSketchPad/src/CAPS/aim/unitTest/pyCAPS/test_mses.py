import unittest

import os
import glob
import shutil

import pyCAPS

class TestMSES_Kulfan(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_msesKulfanTest"
        cls.iProb = 1
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Load CSM file
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","kulfanSection.csm"),
                                       outLevel = 0)

        # Load mses aim
        cls.mses = cls.myProblem.analysis.create(aim = "msesAIM")

        # Set Mach number and GridAlpha
        cls.mses.input.Mach = 0.5
        cls.mses.input.GridAlpha = 3.0

    @classmethod
    def tearDownClass(cls):
        del cls.mses
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
    def test_allInputs(self):
        mses = self.myProblem.analysis.create(aim = "msesAIM")

        mses.input.Mach = 0.5
        mses.input.Re = 1e6
        mses.input.Alpha = None
        mses.input.CL = 0.1
        mses.input.Acrit = 5
        mses.input.xTransition_Upper = 0.1
        mses.input.xTransition_Lower = 0.5
        mses.input.Mcrit = 0.9
        mses.input.MuCon = 1.0
        mses.input.ISMOM = 4
        mses.input.IFFBC = 2
        mses.input.GridAlpha = 1.0
        mses.input.Airfoil_Points = 201
        mses.input.Inlet_Points = 41
        mses.input.Outlet_Points = 41
        mses.input.Upper_Stremlines = 21
        mses.input.Lower_Stremlines = 21
        mses.input.Coarse_Iteration = 100
        mses.input.xGridRange = [-2, 3]
        mses.input.yGridRange = [-2, 3]
        mses.input.Design_Variable = {"auppervar": {},
                                      "alowervar": {},
                                      "classvar": {},
                                      "ztailvar": {}}
        mses.input.Cheby_Modes = [0]*20

#==============================================================================
    def test_execute(self):

        AlphaTrue =  3.0
        ClTrue    =  0.8394186474518254
        CdTrue    =  0.006555029283599085
        CdpTrue   =  0.002216345232296256
        CdvTrue   =  0.006555029283599085
        CdwTrue   =  0.0
        CmTrue    =  -0.10463276571925233
        McritTrue =  0.6281087595405522

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 1e6
        self.mses.input.Mach  = 0.2
        self.mses.input.xTransition_Upper = 1.0
        self.mses.input.xTransition_Lower = 1.0

        # Retrieve results
        Alpha = self.mses.output.Alpha
        Cl    = self.mses.output.CL
        Cd    = self.mses.output.CD
        Cdp   = self.mses.output.CD_p
        Cdv   = self.mses.output.CD_v
        Cdw   = self.mses.output.CD_w
        Cm    = self.mses.output.CM
        Mcrit = self.mses.output.Mcrit
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("CdpTrue   = ", Cdp)
        # print("CdvTrue   = ", Cdv)
        # print("CdwTrue   = ", Cdw)
        # print("CmTrue    = ", Cm)
        # print("McritTrue = ", Mcrit)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue   , Cl   , 3)
        self.assertAlmostEqual(CdTrue   , Cd   , 3)
        self.assertAlmostEqual(CdpTrue  , Cdp  , 3)
        self.assertAlmostEqual(CdvTrue  , Cdv  , 3)
        self.assertAlmostEqual(CdwTrue  , Cdw  , 3)
        self.assertAlmostEqual(CmTrue   , Cm   , 3)
        self.assertAlmostEqual(McritTrue, Mcrit, 2)


        AlphaTrue =  3.601263317256423
        ClTrue    =  0.9
        CdTrue    =  0.007219418142688761
        CdpTrue   =  0.0024180255576456413
        CdvTrue   =  0.007219418142688761
        CdwTrue   =  0.0
        CmTrue    =  -0.10344222164778835
        McritTrue =  0.6133316825117032

        # Set custom Cl
        self.mses.input.GridAlpha = 3.0
        self.mses.input.Alpha = None
        self.mses.input.CL    = ClTrue
        self.mses.input.Re    = 1e6
        self.mses.input.xTransition_Upper = 0.7
        self.mses.input.xTransition_Lower = 0.8

        # Retrieve results
        Alpha = self.mses.output.Alpha
        Cl    = self.mses.output.CL
        Cd    = self.mses.output.CD
        Cdp   = self.mses.output.CD_p
        Cdv   = self.mses.output.CD_v
        Cdw   = self.mses.output.CD_w
        Cm    = self.mses.output.CM
        Mcrit = self.mses.output.Mcrit
        # print("AlphaTrue = ", Alpha)
        # print("ClTrue    = ", Cl)
        # print("CdTrue    = ", Cd)
        # print("CdpTrue   = ", Cdp)
        # print("CdvTrue   = ", Cdv)
        # print("CdwTrue   = ", Cdw)
        # print("CmTrue    = ", Cm)
        # print("McritTrue = ", Mcrit)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue   , Cl   , 3)
        self.assertAlmostEqual(CdTrue   , Cd   , 3)
        self.assertAlmostEqual(CdpTrue  , Cdp  , 3)
        self.assertAlmostEqual(CdvTrue  , Cdv  , 3)
        self.assertAlmostEqual(CdwTrue  , Cdw  , 3)
        self.assertAlmostEqual(CmTrue   , Cm   , 3)
        self.assertAlmostEqual(McritTrue, Mcrit, 3)

#==============================================================================
    def test_Cheby_Modes(self):

        AlphaTrue = 3.0

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 0.0
        self.mses.input.Mach  = 0.2

        self.mses.input.Cheby_Modes = [1e-4]*40
        Cl_mod   = self.mses.output["CL"].deriv("Cheby_Modes")
        self.assertEqual(len(Cl_mod), 40)

        self.mses.input.Cheby_Modes = [1e-4]*20
        Cl_mod   = self.mses.output["CL"].deriv("Cheby_Modes")
        self.assertEqual(len(Cl_mod), 20)

        self.mses.input.Cheby_Modes = None
        self.mses.input.Design_Variable = {"auppervar": {},
                                           "alowervar": {},
                                           "classvar": {},
                                           "ztailvar": {}}

        Cheby_Modes = self.mses.output.Cheby_Modes
        Cheby_Modes_auppervar = self.mses.output["Cheby_Modes"].deriv("auppervar")
        Cheby_Modes_alowervar = self.mses.output["Cheby_Modes"].deriv("alowervar")
        Cheby_Modes_classvar  = self.mses.output["Cheby_Modes"].deriv("classvar")
        Cheby_Modes_ztailvar  = self.mses.output["Cheby_Modes"].deriv("ztailvar")

        self.mses.input.Design_Variable = None

#==============================================================================
    def test_Cheby_Modes_Symmetric(self):
        # For a symmetric airfoil, the Cheby_Mode derivatives w.r.t. thickness
        # should be identical top and bottom

        myProblem = pyCAPS.Problem(problemName = self.problemName+str(self.iProb),
                                   capsFile = os.path.join("..","csmData","kulfanSymmetricSection.csm"),
                                   outLevel = 0); self.iProb += 1

        # Load mses aim
        mses = myProblem.analysis.create(aim = "msesAIM")

        AlphaTrue = 0.0

        # Set custom AoA and Re
        mses.input.Alpha = AlphaTrue
        mses.input.CL    = None
        mses.input.Re    = 1e6
        mses.input.Mach  = 0.2
        mses.input.xTransition_Upper = 0.1
        mses.input.xTransition_Lower = 0.1

        mses.input.Cheby_Modes = [0]*40
        mses.input.Design_Variable = None

        Cheby_Modes = mses.output.Cheby_Modes
        CD_Cheby_Modes = mses.output["CD"].deriv("Cheby_Modes")

        nCheby = int(len(CD_Cheby_Modes)/2)
        for i in range(nCheby):
            self.assertAlmostEqual(CD_Cheby_Modes[i], CD_Cheby_Modes[i+nCheby], 2)


        mses.input.Cheby_Modes = None
        mses.input.Design_Variable = {"avar": {}}

        Cheby_Modes = mses.output.Cheby_Modes
        Cheby_Modes_avar = mses.output["Cheby_Modes"].deriv("avar")

        nCheby = int(len(Cheby_Modes_avar)/2)
        for i in range(1,nCheby):
            for j in range(len(Cheby_Modes_avar[i])):
                self.assertAlmostEqual(Cheby_Modes_avar[i][j], Cheby_Modes_avar[i+nCheby][j], 5)

        mses.input.Design_Variable = None

#==============================================================================
    def test_sensitivity_AnalysisIn(self):

        eps = 1e-4
        AlphaTrue = 3.0
        MachTrue = 0.2
        ReTrue = 0 #1e6

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = ReTrue
        self.mses.input.Mach  = MachTrue

        Cl_Alpha  = self.mses.output["CL"].deriv("Alpha")
        Cd_Alpha  = self.mses.output["CD"].deriv("Alpha")
        Mc_Alpha  = self.mses.output["Mcrit"].deriv("Alpha")

        AlphaTrue += eps
        self.mses.input.Alpha = AlphaTrue
        Clp = self.mses.output.CL
        Cdp = self.mses.output.CD
        Mcp = self.mses.output.Mcrit

        AlphaTrue -= 2*eps
        self.mses.input.Alpha = AlphaTrue
        Clm = self.mses.output.CL
        Cdm = self.mses.output.CD
        Mcm = self.mses.output.Mcrit

        # remove the step change
        AlphaTrue += eps
        self.mses.input.Alpha = AlphaTrue

        # print("******", (Clp - Clm)/(2*eps), Cl_Alpha)
        # print("******", (Cdp - Cdm)/(2*eps), Cd_Alpha)
        # print("******", (Mcp - Mcm)/(2*eps), Mc_Alpha)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_Alpha, 2)
        self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_Alpha, 2)
        self.assertAlmostEqual((Mcp - Mcm)/(2*eps), Mc_Alpha, 2)


        Cl_Mach  = self.mses.output["CL"].deriv("Mach")
        Cd_Mach  = self.mses.output["CD"].deriv("Mach")
        Mc_Mach  = self.mses.output["Mcrit"].deriv("Mach")

        MachTrue += eps
        self.mses.input.Mach = MachTrue
        Clp = self.mses.output.CL
        Cdp = self.mses.output.CD
        Mcp = self.mses.output.Mcrit

        MachTrue -= 2*eps
        self.mses.input.Mach = MachTrue
        Clm = self.mses.output.CL
        Cdm = self.mses.output.CD
        Mcm = self.mses.output.Mcrit

        # remove the step change
        MachTrue += eps
        self.mses.input.Mach = MachTrue

        # print("******", (Clp - Clm)/(2*eps), Cl_Mach)
        # print("******", (Cdp - Cdm)/(2*eps), Cd_Mach)
        # print("******", (Mcp - Mcm)/(2*eps), Mc_Mach)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_Mach, 2)
        self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_Mach, 2)
        self.assertAlmostEqual((Mcp - Mcm)/(2*eps), Mc_Mach, 2)

        Cl_Re = self.mses.output["CL"].deriv("Re")

        return
        eps = 1000
        ReTrue += eps
        self.mses.input.Re = ReTrue
        Clp = self.mses.output.CL

        ReTrue -= 2*eps
        self.mses.input.Re = ReTrue
        Clm = self.mses.output.CL

        # remove the step change
        ReTrue += eps
        self.mses.input.Re = ReTrue

        #print("******", (Clp - Clm)/(2*eps), Cl_Re)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_Re, 2)

#==============================================================================
    def off_test_geom_sensitivity_CL(self):

        eps = 1e-4
        AlphaTrue = 3.0

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 0

        # self.mses.input.Design_Variable = {"auppervar": {}}
        #
        # Cl = self.mses.output.CL
        # Cl_auppervar  = self.mses.output["CL"].deriv("auppervar")
        #
        # self.mses.input.Design_Variable = None
        # auppervar  = self.myProblem.geometry.despmtr.auppervar
        #
        # print("auppervar")
        # for i in range(len(auppervar)):
        #     auppervar[i] += eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #     Clp = self.mses.output.CL
        #
        #     auppervar[i] -= 2*eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #     Clm = self.mses.output.CL
        #
        #     # remove the step change
        #     auppervar[i] += eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #
        #     print("******", (Clp - Clm)/(2*eps), Cl_auppervar[i])
        #     #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_auppervar[i], 3)
        #
        # return

        self.mses.input.Design_Variable = {"auppervar": {},
                                           "alowervar": {},
                                           "classvar": {}}


        Cl = self.mses.output.CL
        Cl_auppervar = self.mses.output["CL"].deriv("auppervar")
        Cl_alowervar = self.mses.output["CL"].deriv("alowervar")
        Cl_classvar  = self.mses.output["CL"].deriv("classvar")

        self.mses.input.Design_Variable = None
        auppervar = self.myProblem.geometry.despmtr.auppervar
        alowervar = self.myProblem.geometry.despmtr.alowervar
        classvar  = self.myProblem.geometry.despmtr.classvar

        print("auppervar")
        for i in range(len(auppervar)):
            auppervar[i] += eps
            self.myProblem.geometry.despmtr.auppervar = auppervar
            Clp = self.mses.output.CL

            auppervar[i] -= 2*eps
            self.myProblem.geometry.despmtr.auppervar = auppervar
            Clm = self.mses.output.CL

            # remove the step change
            auppervar[i] += eps
            self.myProblem.geometry.despmtr.auppervar = auppervar

            print("******", (Clp - Clm)/(2*eps), Cl_auppervar[i])
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_auppervar[i], 3)

        print("alowervar")
        for i in range(len(alowervar)):
            alowervar[i] += eps
            self.myProblem.geometry.despmtr.alowervar = alowervar
            Clp = self.mses.output.CL

            alowervar[i] -= 2*eps
            self.myProblem.geometry.despmtr.alowervar = alowervar
            Clm = self.mses.output.CL

            # remove the step change
            alowervar[i] += eps
            self.myProblem.geometry.despmtr.alowervar = alowervar

            print("******", (Clp - Clm)/(2*eps), Cl_alowervar[i])
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_alowervar[i], 3)

        print("classvar")
        for i in range(len(classvar)):
            classvar[i] += eps
            self.myProblem.geometry.despmtr.classvar = classvar
            Clp = self.mses.output.CL

            classvar[i] -= 2*eps
            self.myProblem.geometry.despmtr.classvar = classvar
            Clm = self.mses.output.CL

            # remove the step change
            classvar[i] += eps
            self.myProblem.geometry.despmtr.classvar = classvar

            print("******", (Clp - Clm)/(2*eps), Cl_classvar[i])
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_classvar[i], 3)

#==============================================================================
    def off_test_geom_sensitivity_CD(self):

        eps = 1e-4
        AlphaTrue = 3.0

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 1e7

        self.mses.input.Design_Variable = {"auppervar": {}}

        Cd = self.mses.output.CD
        Cd_auppervar  = self.mses.output["CD"].deriv("auppervar")

        self.mses.input.Design_Variable = None
        auppervar  = self.myProblem.geometry.despmtr.auppervar

        # print("auppervar")
        # for i in range(len(auppervar)):
        #     auppervar[i] += eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #     Cdp = self.mses.output.CD
        #
        #     auppervar[i] -= 2*eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #     Cdm = self.mses.output.CD
        #
        #     # remove the step change
        #     auppervar[i] += eps
        #     self.myProblem.geometry.despmtr.auppervar = auppervar
        #
        #     print("******", (Cdp - Cdm)/(2*eps), Cd_auppervar[i])
        #     #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_auppervar[i], 3)
        #
        # return

        self.mses.input.Design_Variable = {"auppervar": {},
                                           "alowervar": {},
                                           "classvar": {}}


        Cd = self.mses.output.CD
        Cd_auppervar = self.mses.output["CD"].deriv("auppervar")
        Cd_alowervar = self.mses.output["CD"].deriv("alowervar")
        Cd_classvar  = self.mses.output["CD"].deriv("classvar")

        self.mses.input.Design_Variable = None
        auppervar = self.myProblem.geometry.despmtr.auppervar
        alowervar = self.myProblem.geometry.despmtr.alowervar
        classvar  = self.myProblem.geometry.despmtr.classvar

        print("auppervar")
        for i in range(len(auppervar)):
            auppervar[i] += eps
            self.myProblem.geometry.despmtr.auppervar = auppervar
            Cdp = self.mses.output.CD

            auppervar[i] -= 2*eps
            self.myProblem.geometry.despmtr.auppervar = auppervar
            Cdm = self.mses.output.CD

            # remove the step change
            auppervar[i] += eps
            self.myProblem.geometry.despmtr.auppervar = auppervar

            print("******", (Cdp - Cdm)/(2*eps), Cd_auppervar[i])
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_auppervar[i], 3)

        print("alowervar")
        for i in range(len(alowervar)):
            alowervar[i] += eps
            self.myProblem.geometry.despmtr.alowervar = alowervar
            Cdp = self.mses.output.CD

            alowervar[i] -= 2*eps
            self.myProblem.geometry.despmtr.alowervar = alowervar
            Cdm = self.mses.output.CD

            # remove the step change
            alowervar[i] += eps
            self.myProblem.geometry.despmtr.alowervar = alowervar

            print("******", (Cdp - Cdm)/(2*eps), Cd_alowervar[i])
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_alowervar[i], 3)

        print("classvar")
        for i in range(len(classvar)):
            classvar[i] += eps
            self.myProblem.geometry.despmtr.classvar = classvar
            Cdp = self.mses.output.CD

            classvar[i] -= 2*eps
            self.myProblem.geometry.despmtr.classvar = classvar
            Cdm = self.mses.output.CD

            # remove the step change
            classvar[i] += eps
            self.myProblem.geometry.despmtr.classvar = classvar

            print("******", (Cdp - Cdm)/(2*eps), Cd_classvar[i])
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_classvar[i], 3)


###############################################################################
class TestMSES_NACA(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        # Create working directory variable
        cls.problemName = "workDir_msesNACATest"
        cls.iProb = 1
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Load CSM file
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","airfoilSection.csm"),
                                       outLevel = 0)

        # Load mses aim
        cls.mses = cls.myProblem.analysis.create(aim = "msesAIM")

        # Set Mach number
        cls.mses.input.Mach = 0.5

    @classmethod
    def tearDownClass(cls):
        del cls.mses
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
    def test_Cheby_Modes_Symmetric(self):
        # For a symmetric airfoil, the Cheby_Mode derivatives w.r.t. thickness
        # should be identical top and bottom

        AlphaTrue = 0.0

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 0.0
        self.mses.input.Mach  = 0.2

        self.mses.input.Cheby_Modes = None
        self.mses.input.Design_Variable = {"thick": {}}

        Cheby_Modes = self.mses.output.Cheby_Modes
        Cheby_Modes_thick = self.mses.output["Cheby_Modes"].deriv("thick")

        nCheby = int(len(Cheby_Modes_thick)/2)
        for i in range(nCheby):
            self.assertAlmostEqual(Cheby_Modes_thick[i], Cheby_Modes_thick[i+nCheby], 5)

        self.mses.input.Design_Variable = None


#==============================================================================
    def test_geom_sensitivity(self):

        eps = 1e-4
        AlphaTrue = 3.0

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 0

        self.mses.input.Design_Variable = {"thick": {},
                                           "camber": {}}


        Cl_thick  = self.mses.output["CL"].deriv("thick")
        Cl_camber = self.mses.output["CL"].deriv("camber")
        Mc_thick  = self.mses.output["Mcrit"].deriv("thick")
        Mc_camber = self.mses.output["Mcrit"].deriv("camber")

        self.mses.input.Design_Variable = None
        thick = self.myProblem.geometry.despmtr.thick
        camber = self.myProblem.geometry.despmtr.camber

        #-------------------------------------------
        thick += eps
        self.myProblem.geometry.despmtr.thick = thick
        Clp = self.mses.output.CL
        Mcp = self.mses.output.Mcrit

        thick -= 2*eps
        self.myProblem.geometry.despmtr.thick = thick
        Clm = self.mses.output.CL
        Mcm = self.mses.output.Mcrit

        # remove the step change
        thick += eps
        self.myProblem.geometry.despmtr.thick = thick

        # print("thick")
        # print("******", (Clp - Clm)/(2*eps), Cl_thick)
        # print("******", (Mcp - Mcm)/(2*eps), Mc_thick)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_thick, 2)
        self.assertAlmostEqual((Mcp - Mcm)/(2*eps), Mc_thick, 1)
        #-------------------------------------------

        #-------------------------------------------
        camber += eps
        self.myProblem.geometry.despmtr.camber = camber
        Clp = self.mses.output.CL
        Mcp = self.mses.output.Mcrit

        camber -= 2*eps
        self.myProblem.geometry.despmtr.camber = camber
        Clm = self.mses.output.CL
        Mcm = self.mses.output.Mcrit

        # remove the step change
        camber += eps
        self.myProblem.geometry.despmtr.camber = camber

        # print("camber")
        # print("******", (Clp - Clm)/(2*eps), Cl_camber)
        # print("******", (Mcp - Mcm)/(2*eps), Mc_camber)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_camber, 2)
        self.assertAlmostEqual((Mcp - Mcm)/(2*eps), Mc_camber, 1)
        #-------------------------------------------

#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load AIM
        if verbose: print(6*"-","Load msesAIM", line)
        mses = myProblem.analysis.create(aim = "msesAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())


        # Set custom AoA and Re
        if verbose: print(6*"-","Modify Alpha", line)
        mses.input.Alpha = 0.0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Re", line)
        mses.input.Re    = 1.0e6; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Mach", line)
        mses.input.Mach  = 0.2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Design_Variable", line)
        mses.input.Design_Variable = {"thick": {}}; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Get Cheby_Modes_thick", line)
        Cheby_Modes_thick = mses.output["Cheby_Modes"].deriv("thick"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        nCheby = int(len(Cheby_Modes_thick)/2)
        for i in range(nCheby):
            self.assertAlmostEqual(Cheby_Modes_thick[i], Cheby_Modes_thick[i+nCheby], 4)

        if verbose: print(6*"-","Get CL", line)
        CL1 = mses.output.CL; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Modify Alpha", line)
        mses.input.Alpha = 2.0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-","Get CL", line)
        CL2 = mses.output.CL; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        self.assertGreater(CL2, CL1, 3)

        if verbose: print(6*"-","Get Cheby_Modes_thick", line)
        Cheby_Modes_thick = mses.output["Cheby_Modes"].deriv("thick"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        nCheby = int(len(Cheby_Modes_thick)/2)
        for i in range(nCheby):
            self.assertAlmostEqual(Cheby_Modes_thick[i], Cheby_Modes_thick[i+nCheby], 4)

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        capsFile = os.path.join("..","csmData","airfoilSection.csm")
        problemName = self.problemName+str(self.iProb)

        myProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_journal(myProblem, -1)

        myProblem.close()
        shutil.rmtree(problemName)

        #print(80*"=")
        #print(80*"=")
        # Create the problem to start journaling
        myProblem = pyCAPS.Problem(problemName, capsFile=capsFile, outLevel=0)
        myProblem.close()

        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile=capsFile, outLevel=0)
            self.run_journal(myProblem, line_exit)
            myProblem.close()

        self.__class__.iProb += 1

if __name__ == '__main__':
    unittest.main()
