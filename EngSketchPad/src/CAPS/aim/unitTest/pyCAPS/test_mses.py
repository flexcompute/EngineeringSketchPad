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
        #workDir = os.path.join(str(args.workDir[0]), workDir)

        cls.cleanUp()

        # Load CSM file
        cls.myProblem = pyCAPS.Problem(problemName = cls.problemName,
                                       capsFile = os.path.join("..","csmData","kulfanSection.csm"),
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
    def test_execute(self):

        AlphaTrue = 3.0
        ClTrue    = 0.8404364788053338
        CdTrue    = 0.006555756386418418
        CdpTrue   = 0.002221088763744491
        CdvTrue   = 0.006555756386418418
        CdwTrue   = 0.0
        CmTrue    = -0.10472886188421521

        # Set custom AoA and Re
        self.mses.input.Alpha = AlphaTrue
        self.mses.input.CL    = None
        self.mses.input.Re    = 1e6
        self.mses.input.Mach  = 0.2

        # Retrieve results
        Alpha = self.mses.output.Alpha
        Cl    = self.mses.output.CL
        Cd    = self.mses.output.CD
        Cdp   = self.mses.output.CD_p
        Cdv   = self.mses.output.CD_v
        Cdw   = self.mses.output.CD_w
        Cm    = self.mses.output.CM
        # print("Alpha = ", Alpha)
        # print("Cl  = ", Cl)
        # print("Cd  = ", Cd)
        # print("Cdp = ", Cdp)
        # print("Cdv = ", Cdv)
        # print("Cdw = ", Cdw)
        # print("Cm  = ", Cm)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue , Cl , 3)
        self.assertAlmostEqual(CdTrue , Cd , 3)
        self.assertAlmostEqual(CdpTrue, Cdp, 3)
        self.assertAlmostEqual(CdvTrue, Cdv, 3)
        self.assertAlmostEqual(CdwTrue, Cdw, 3)
        self.assertAlmostEqual(CmTrue , Cm , 3)


        AlphaTrue = 3.5586937597499455
        ClTrue    = 0.9
        CdTrue    = 0.0068401834989933906
        CdpTrue   = 0.002417738651453858
        CdvTrue   = 0.0068401834989933906
        CdwTrue   = 0.0
        CmTrue    = -0.10404788841883214

        # Set custom AoA and Re
        self.mses.input.Alpha = None
        self.mses.input.CL    = ClTrue
        self.mses.input.Re    = 1e6

        # Retrieve results
        Alpha = self.mses.output.Alpha
        Cl    = self.mses.output.CL
        Cd    = self.mses.output.CD
        Cdp   = self.mses.output.CD_p
        Cdv   = self.mses.output.CD_v
        Cdw   = self.mses.output.CD_w
        Cm    = self.mses.output.CM
        # print("Alpha = ", Alpha)
        # print("Cl  = ", Cl)
        # print("Cd  = ", Cd)
        # print("Cdp = ", Cdp)
        # print("Cdv = ", Cdv)
        # print("Cdw = ", Cdw)
        # print("Cm  = ", Cm)

        self.assertAlmostEqual(AlphaTrue, Alpha, 3)
        self.assertAlmostEqual(ClTrue , Cl , 3)
        self.assertAlmostEqual(CdTrue , Cd , 3)
        self.assertAlmostEqual(CdpTrue, Cdp, 3)
        self.assertAlmostEqual(CdvTrue, Cdv, 3)
        self.assertAlmostEqual(CdwTrue, Cdw, 3)
        self.assertAlmostEqual(CmTrue , Cm , 3)


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

        Cl = self.mses.output.CL
        Cl_Alpha  = self.mses.output["CL"].deriv("Alpha")

        AlphaTrue += eps
        self.mses.input.Alpha = AlphaTrue
        Clp = self.mses.output.CL

        AlphaTrue -= 2*eps
        self.mses.input.Alpha = AlphaTrue
        Clm = self.mses.output.CL

        # remove the step change
        AlphaTrue += eps
        self.mses.input.Alpha = AlphaTrue

        #print("******", (Clp - Clm)/(2*eps), Cl_Alpha)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_Alpha, 2)


        Cl_Mach  = self.mses.output["CL"].deriv("Mach")

        MachTrue += eps
        self.mses.input.Mach = MachTrue
        Clp = self.mses.output.CL

        MachTrue -= 2*eps
        self.mses.input.Mach = MachTrue
        Clm = self.mses.output.CL

        # remove the step change
        MachTrue += eps
        self.mses.input.Mach = MachTrue

        #print("******", (Clp - Clm)/(2*eps), Cl_Mach)
        self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_Mach, 2)

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
        #     #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_auppervar[i], 4)
        #
        # return

        self.mses.input.Design_Variable = {"auppervar": {},
                                           "alowervar": {},
                                           "classvar": {},
                                           "ztailvar": {}}


        Cl = self.mses.output.CL
        Cl_auppervar = self.mses.output["CL"].deriv("auppervar")
        Cl_alowervar = self.mses.output["CL"].deriv("alowervar")
        Cl_classvar  = self.mses.output["CL"].deriv("classvar")
        Cl_ztailvar  = self.mses.output["CL"].deriv("ztailvar")

        self.mses.input.Design_Variable = None
        auppervar = self.myProblem.geometry.despmtr.auppervar
        alowervar = self.myProblem.geometry.despmtr.alowervar
        classvar  = self.myProblem.geometry.despmtr.classvar
        ztailvar  = self.myProblem.geometry.despmtr.ztailvar

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
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_auppervar[i], 4)

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
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_alowervar[i], 4)

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
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_classvar[i], 4)

        print("ztailvar")
        for i in range(len(ztailvar)):
            ztailvar[i] += eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar
            Clp = self.mses.output.CL

            ztailvar[i] -= 2*eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar
            Clm = self.mses.output.CL

            # remove the step change
            ztailvar[i] += eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar

            print("******", (Clp - Clm)/(2*eps), Cl_ztailvar[i])
            #self.assertAlmostEqual((Clp - Clm)/(2*eps), Cl_ztailvar[i], 4)


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
        #     #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_auppervar[i], 4)
        #
        # return

        self.mses.input.Design_Variable = {"auppervar": {},
                                           "alowervar": {},
                                           "classvar": {},
                                           "ztailvar": {}}


        Cd = self.mses.output.CD
        Cd_auppervar = self.mses.output["CD"].deriv("auppervar")
        Cd_alowervar = self.mses.output["CD"].deriv("alowervar")
        Cd_classvar  = self.mses.output["CD"].deriv("classvar")
        Cd_ztailvar  = self.mses.output["CD"].deriv("ztailvar")

        self.mses.input.Design_Variable = None
        auppervar = self.myProblem.geometry.despmtr.auppervar
        alowervar = self.myProblem.geometry.despmtr.alowervar
        classvar  = self.myProblem.geometry.despmtr.classvar
        ztailvar  = self.myProblem.geometry.despmtr.ztailvar

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
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_auppervar[i], 4)

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
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_alowervar[i], 4)

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
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_classvar[i], 4)

        print("ztailvar")
        for i in range(len(ztailvar)):
            ztailvar[i] += eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar
            Cdp = self.mses.output.CD

            ztailvar[i] -= 2*eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar
            Cdm = self.mses.output.CD

            # remove the step change
            ztailvar[i] += eps
            self.myProblem.geometry.despmtr.ztailvar = ztailvar

            print("******", (Cdp - Cdm)/(2*eps), Cd_ztailvar[i])
            #self.assertAlmostEqual((Cdp - Cdm)/(2*eps), Cd_ztailvar[i], 4)

if __name__ == '__main__':
    unittest.main()
