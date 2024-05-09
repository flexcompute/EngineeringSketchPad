import unittest

import os
import glob
import shutil

import sys

import pyCAPS

class TestCBAERO(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.problemName = "workDir_cbaero"
        cls.iProb = 1
        cls.cleanUp()

        cornerFile = os.path.join("..","csmData","cfdSingleBody.csm")
        cls.capsProblem = pyCAPS.Problem(cls.problemName, capsFile=cornerFile)
        
        cls.capsProblem.geometry.cfgpmtr.farfield = 0

        cls.aflr4 = cls.capsProblem.analysis.create(aim = "aflr4AIM",
                                              name = "aflr4")

        # Modify local mesh sizing parameters
        cls.aflr4.input.Mesh_Quiet_Flag = True
        cls.aflr4.input.Mesh_Length_Factor = 20

    @classmethod
    def tearDownClass(cls):
        del cls.capsProblem
        del cls.aflr4
        cls.cleanUp()
        pass

    @classmethod
    def cleanUp(cls):

        # Remove analysis directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
    # Test inputs
    def test_inputs(self):

        m   = pyCAPS.Unit("meter")
        deg = pyCAPS.Unit("degree")
        psi = pyCAPS.Unit("psi")

        cbaero = self.capsProblem.analysis.create(aim = "cbaeroAIM",autoExec=False)

        cbaero.input.Proj_Name = "cbaero_CAPS"
        cbaero.input.Mach = 0.5
        cbaero.input.Alpha = 10 * deg
        cbaero.input.Beta = 5 * deg
        cbaero.input.ReferenceArea = 2 * m*m
        cbaero.input.ReferenceChord = 0.6 * m
        cbaero.input.ReferenceSpan = 6 * m
        cbaero.input.Moment_Center = [1,0.1,0.2] * m
        cbaero.input.Flow_Type = "Laminar"
        cbaero.input.Critical_Transition = 9
        cbaero.input.Planet = "Jupiter"
        cbaero.input.Default_Body_Method = "TangentCone"
        cbaero.input.Default_Body_Method = "TangentWedge"
        cbaero.input.Default_Low_Speed_Method = "LowAR"
        cbaero.input.Leading_Edge_Suction = 0.8
        cbaero.input.Aero_Surface = {"Wing1":"Wing"}
        cbaero.input.Mesh_Morph = False
        cbaero.input["Surface_Mesh"].link(self.aflr4.output["Surface_Mesh"])

        # Write input files
        cbaero.preAnalysis()

#==============================================================================
    def test_Output(self):

        m      = pyCAPS.Unit("meter")
        cm     = pyCAPS.Unit("centimeter")
        watt   = pyCAPS.Unit("watt")
        kelvin = pyCAPS.Unit("kelvin")
        deg    = pyCAPS.Unit("degree")
        bar    = pyCAPS.Unit("bar")

        cbaero = self.capsProblem.analysis.create(aim = "cbaeroAIM", autoExec=False)

        cbaero.input["Surface_Mesh"].link(self.aflr4.output["Surface_Mesh"])

        # Write input files
        cbaero.preAnalysis()

        # Create a dummy plt file
        filename = os.path.join(cbaero.analysisDir, cbaero.input.Proj_Name+".plt")
        with open(filename, "w") as f:
            f.write("  Beta       Alpha        Bars      Mach      PerTrb       CL          CD        Cm        LoD        CLp        CDp         CLf        CDf     StagTemp   QdotConv    QdotRad  StagRadius CL_Trefftz CD_Trefftz\n")
            f.write("BLOCK B0.000000_Q0.689476_M20.000000\n")
            f.write("   0.00001    0.00000  0.6894757   20.00000   15.71943    0.00004    1.65778   -0.00000    0.00002    0.00004    1.65693    0.00000    0.00086    6638.99       0.10       0.01       0.10    0.10000    4.00000\n")
            f.write("   0.00020    5.00000  0.6894757   20.00000   11.42210   -0.11959    1.64205   -0.02117   -0.07283   -0.11980    1.64120    0.00021    0.00086    6638.23       0.20       0.02       0.01    0.02000    0.30000\n")
            f.write("   0.00300    7.00000  0.6894757   20.00000   12.70497   -0.16611    1.62714   -0.02950   -0.10209   -0.16641    1.62622    0.00029    0.00092    6639.11       0.30       0.03       0.20    0.00300    0.00200\n")
            f.write("   0.04000   10.00000  0.6894757   20.00000   13.25314   -0.23335    1.59578   -0.04172   -0.14623   -0.23374    1.59478    0.00040    0.00100    6635.25       0.40       0.04       0.02    0.00040    0.00010\n")
            f.write("\n")
            f.write("\n")
            f.write("  Beta       Alpha        Bars      Mach      PerTrb       CL          CD        Cm        LoD        CLp        CDp         CLf        CDf     StagTemp   QdotConv    QdotRad  StagRadius CL_Trefftz CD_Trefftz\n")
            f.write("BLOCK B0.000000_Q0.689476_M25.000000\n")
            f.write("   0.00000    0.00000  0.6894757   25.00000    1.92223    0.00004    1.66538   -0.00000    0.00002    0.00004    1.66467    0.00000    0.00070    7685.48       0.00       0.00       0.00    0.00000    0.00000\n")
            f.write("   0.00000    5.00000  0.6894757   25.00000    1.67473   -0.12021    1.64958   -0.02130   -0.07287   -0.12035    1.64887    0.00014    0.00071    7684.17       0.00       0.00       0.00    0.00000    0.00000\n")
            f.write("   0.00000    7.00000  0.6894757   25.00000    3.20561   -0.16694    1.63458   -0.02966   -0.10213   -0.16717    1.63382    0.00022    0.00076    7685.68       0.00       0.00       0.00    0.00000    0.00000\n")
            f.write("   0.00000   10.00000  0.6894757   25.00000    4.44776   -0.23449    1.60305   -0.04195   -0.14628   -0.23481    1.60222    0.00031    0.00083    7679.05       0.00       0.00       0.00    0.00000    0.00000\n")

        cbaero.postAnalysis()

        TrueValues = {"Beta": [1e-05, 0.0002, 0.003, 0.04, 0.0, 0.0, 0.0, 0.0]*deg,
                      "Alpha": [0.0, 5.0, 7.0, 10.0, 0.0, 5.0, 7.0, 10.0]*deg,
                      "Dynamic_Pressure": [0.6894757, 0.6894757, 0.6894757, 0.6894757, 0.6894757, 0.6894757, 0.6894757, 0.6894757]*bar,
                      "Mach": [20.0, 20.0, 20.0, 20.0, 25.0, 25.0, 25.0, 25.0],
                      "PerTrb": [15.71943, 11.4221, 12.70497, 13.25314, 1.92223, 1.67473, 3.20561, 4.44776],
                      "CLtot": [4e-05, -0.11959, -0.16611, -0.23335, 4e-05, -0.12021, -0.16694, -0.23449],
                      "CDtot": [1.65778, 1.64205, 1.62714, 1.59578, 1.66538, 1.64958, 1.63458, 1.60305],
                      "CMYtot": [-0.0, -0.02117, -0.0295, -0.04172, -0.0, -0.0213, -0.02966, -0.04195],
                      "LoDtot": [2e-05, -0.07283, -0.10209, -0.14623, 2e-05, -0.07287, -0.10213, -0.14628],
                      "CL_p": [4e-05, -0.1198, -0.16641, -0.23374, 4e-05, -0.12035, -0.16717, -0.23481],
                      "CD_p": [1.65693, 1.6412, 1.62622, 1.59478, 1.66467, 1.64887, 1.63382, 1.60222],
                      "CL_v": [0.0, 0.00021, 0.00029, 0.0004, 0.0, 0.00014, 0.00022, 0.00031],
                      "CD_v": [0.00086, 0.00086, 0.00092, 0.001, 0.0007, 0.00071, 0.00076, 0.00083],
                      "Stagnation_Temperature": [6638.99, 6638.23, 6639.11, 6635.25, 7685.48, 7684.17, 7685.68, 7679.05]*kelvin,
                      "Stagnation_Radius": [0.1, 0.01, 0.2, 0.02, 0.0, 0.0, 0.0, 0.0]*m,
                      "Convective_Flux": [0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.0]*watt / cm**2,
                      "Radiative_Flux": [0.01, 0.02, 0.03, 0.04, 0.0, 0.0, 0.0, 0.0]*watt / cm**2,
                      "CL_Trefftz": [0.1, 0.02, 0.003, 0.0004, 0.0, 0.0, 0.0, 0.0],
                      "CD_Trefftz": [4.0, 0.3, 0.002, 0.0001, 0.0, 0.0, 0.0, 0.0]}

        self.assertEqual(TrueValues.keys(), cbaero.output.keys())

        for key in cbaero.output.keys():
            self.assertAlmostEqual(TrueValues[key], cbaero.output[key].value, 1e-5)

if __name__ == '__main__':
    unittest.main()
