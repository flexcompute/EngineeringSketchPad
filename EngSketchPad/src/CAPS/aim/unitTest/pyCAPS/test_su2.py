import unittest

import os
import glob
import shutil
import __main__

import pyCAPS

class TestSU2(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = os.path.join("..","csmData","cfdMultiBody.csm")
        cls.problemName = "su2_problem"

        cls.cleanUp()

        cls.configFile = "su2_CAPS.cfg"
        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file)

        cls.myGeometry = cls.myProblem.geometry

        egadsTess = cls.myProblem.analysis.create(aim = "egadsTessAIM")

        # Modify local mesh sizing parameters
        egadsTess.input.Tess_Params = [0.5, 0.01, 15.0]
        Mesh_Sizing = {"Farfield": {"tessParams" : [7.,  2., 20.0]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        tetgen = cls.myProblem.analysis.create(aim = "tetgenAIM",
                                               name = "tetgen")

        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        cls.myAnalysis = cls.myProblem.analysis.create(aim = "su2AIM")

        cls.myAnalysis.input["Mesh"].link(tetgen.output["Volume_Mesh"])

    @classmethod
    def tearDownClass(cls):
        del cls.myAnalysis
        del cls.myGeometry
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):
        # Remove default problemName directory
        if os.path.isdir(cls.problemName):
            shutil.rmtree(cls.problemName)

    # Test inputs
    def test_inputs(self):
        
        self.myAnalysis.input.Proj_Name = "su2_CAPS"
        self.myAnalysis.input.Mach = 0.5
        self.myAnalysis.input.Re = 1e6
        self.myAnalysis.input.Physical_Problem = "Navier_Stokes"
        self.myAnalysis.input.Equation_Type = "Incompressible"
        self.myAnalysis.input.Turbulence_Model = "SA_NEG"
        self.myAnalysis.input.Alpha = 10
        self.myAnalysis.input.Beta = 4
        self.myAnalysis.input.Init_Option = "TD_CONDITIONS"
        self.myAnalysis.input.Overwrite_CFG = True
        self.myAnalysis.input.Num_Iter = 5
        self.myAnalysis.input.CFL_Number = 4
        self.myAnalysis.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"},
                                                    "Wing2": {"bcType" : "Inviscid"},
                                                    "Farfield":"farfield"}
        self.myAnalysis.input.MultiGrid_Level = 3
        self.myAnalysis.input.Residual_Reduction = 10
        self.myAnalysis.input.Unit_System = "SI"
        self.myAnalysis.input.Reference_Dimensionalization = "Dimensional"
        self.myAnalysis.input.Freestream_Pressure = 101325
        self.myAnalysis.input.Freestream_Temperature = 300
        self.myAnalysis.input.Freestream_Density = 1.225
        self.myAnalysis.input.Freestream_Velocity = 100
        self.myAnalysis.input.Freestream_Viscosity = 1e-5
        self.myAnalysis.input.Moment_Center = [10,2,4]
        self.myAnalysis.input.Moment_Length = 2
        self.myAnalysis.input.Reference_Area = 10
        self.myAnalysis.input.Pressure_Scale_Factor = 5
        self.myAnalysis.input.Pressure_Scale_Offset = 9
        self.myAnalysis.input.Output_Format = "Tecplot"
        self.myAnalysis.input.Two_Dimensional = False
        self.myAnalysis.input.Convective_Flux = "JST"
        self.myAnalysis.input.Surface_Monitor = ["Wing1", "Wing2"]
        self.myAnalysis.input.Surface_Deform = ["Wing1", "Wing2"]
        self.myAnalysis.input.Input_String = ["value1 = 1", "value2 = 2"]
        
        for version in ["Cardinal", "Raven", "Falcon", "Blackbird"]:
            
            self.myAnalysis.input.SU2_Version = version
            
            self.myAnalysis.preAnalysis()
            self.myAnalysis.postAnalysis()
        
            self.assertEqual(os.path.isfile(os.path.join(self.myAnalysis.analysisDir, self.configFile)), True)
            os.remove(os.path.join(self.myAnalysis.analysisDir, self.configFile))

    # Try put an invalid boundary name
    def test_invalidBoundaryName(self):
        # Create a new instance
        self.su2 = self.myProblem.analysis.create(aim = "su2AIM")
        
        self.su2.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"}, # No capsGroup 'X'
                                             "X": {"bcType" : "Inviscid"},
                                             "Farfield":"farfield"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.su2.preAnalysis()
 
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")
        del self.su2

    # Try an invalid boundary type
    def test_invalidBoundary(self):
        # Create a new instance
        self.su2 = self.myProblem.analysis.create(aim = "su2AIM")
        
        self.su2.input.Boundary_Condition = {"Wing1": {"bcType" : "X"}}
        
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.su2.preAnalysis()
 
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")
        del self.su2
        
    # Test re-enter 
    def test_reenter(self):
        
        self.myAnalysis.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"},
                                                    "Wing2": {"bcType" : "Inviscid"},
                                                    "Farfield":"farfield"}
        self.myAnalysis.preAnalysis()
        self.myAnalysis.postAnalysis()
        
        self.assertEqual(os.path.isfile(os.path.join(self.myAnalysis.analysisDir, self.configFile)), True)
        
        os.remove(os.path.join(self.myAnalysis.analysisDir, self.configFile))
        
        self.myAnalysis.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                                    "Wing2": {"bcType" : "Inviscid"},
                                                    "Farfield":"farfield"}
        self.myAnalysis.preAnalysis()
        self.myAnalysis.postAnalysis()
        
        self.assertEqual(os.path.isfile(os.path.join(self.myAnalysis.analysisDir, self.configFile)), True)
 
    def test_units(self):

        m  = pyCAPS.Unit("meter")
        kg = pyCAPS.Unit("kg")
        s  = pyCAPS.Unit("s")
        K  = pyCAPS.Unit("K")
        kPa  = pyCAPS.Unit("kPa")

        unitSystem={"mass" : kg, "length": m, "time":s, "temperature":K}

        myAnalysis = self.myProblem.analysis.create(aim = "su2AIM", unitSystem=unitSystem)

        myAnalysis.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        myAnalysis.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        myAnalysis.input.Freestream_Density  = 1.225*kg/m**3
        myAnalysis.input.Freestream_Velocity = 200*m/s
        myAnalysis.input.Freestream_Pressure = 101.325*kPa

        myAnalysis.preAnalysis()
 
if __name__ == '__main__':
    unittest.main()
