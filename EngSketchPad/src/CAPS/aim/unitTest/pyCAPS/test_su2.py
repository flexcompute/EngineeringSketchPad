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
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)

#==============================================================================
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

#==============================================================================
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

#==============================================================================
    # Try an invalid boundary type
    def test_invalidBoundary(self):
        # Create a new instance
        self.su2 = self.myProblem.analysis.create(aim = "su2AIM")
        
        self.su2.input.Boundary_Condition = {"Wing1": {"bcType" : "X"}}
        
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.su2.preAnalysis()
 
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")
        del self.su2
        
#==============================================================================
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
 
#==============================================================================
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

#==============================================================================
    def test_phase(self):
        file = os.path.join("..","csmData","cfdMultiBody.csm")
        
        probemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase0", capsFile=file, outLevel=0)

        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM",
                                              name = "egadsTess")

        egadsTess.input.Mesh_Quiet_Flag = True

        # Modify local mesh sizing parameters
        Mesh_Sizing = {"Farfield": {"tessParams" : [0.3*80, 0.2*80, 30]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        tetgen = myProblem.analysis.create(aim = "tetgenAIM",
                                           name = "tetgen")

        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        tetgen.input.Mesh_Quiet_Flag = True

        su2 = myProblem.analysis.create(aim = "su2AIM",
                                        name = "su2")

        su2.input["Mesh"].link(tetgen.output["Volume_Mesh"])

        su2.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"},
                                        "Wing2": {"bcType" : "Inviscid"},
                                        "Farfield":"farfield"}
        su2.preAnalysis()
        su2.postAnalysis()

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(probemName, phaseName="Phase1", prevPhaseName="Phase0", outLevel=0)

        egadsTess = myProblem.analysis["egadsTess"]
        tetgen    = myProblem.analysis["tetgen"]
        su2       = myProblem.analysis["su2"]
        
        su2.input.Alpha = 5
        
        su2.preAnalysis()
        su2.postAnalysis()

#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Load egadsAIM")
        egadsTess = myProblem.analysis.create(aim = "egadsTessAIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess Mesh_Quiet_Flag")
        egadsTess.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Modify local mesh sizing parameters
        if verbose: print(6*"-", line,"egadsTess Tess_Params")
        egadsTess.input.Tess_Params = [0.5, 0.01, 15.0]; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"egadsTess Mesh_Sizing")
        Mesh_Sizing = {"Farfield": {"tessParams" : [7.,  2., 20.0]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        #if verbose: print(6*"-", line,"egadsTess runAnalysis")
        #egadsTess.runAnalysis(); line += 1
        #if line == line_exit: return line
        #if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Load tetgenAIM")
        tetgen = myProblem.analysis.create(aim = "tetgenAIM",
                                           name = "tetgen"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"tetgen Mesh_Quiet_Flag")
        tetgen.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Link Surface_Mesh")
        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        #if verbose: print(6*"-", line,"tetgen runAnalysis")
        #tetgen.runAnalysis(); line += 1
        #if line == line_exit: return line
        #if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Load su2AIM")
        su2 = myProblem.analysis.create(aim = "su2AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Link Mesh")
        su2.input["Mesh"].link(tetgen.output["Volume_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify Boundary_Condition")
        su2.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                        "Wing2": {"bcType" : "Inviscid"},
                                        "Farfield":"farfield"}; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Run the first time
        if verbose: print(6*"-", line,"su2 preAnalysis")
        su2.preAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"su2 system")
        su2.system("echo 'Running SU2 1'"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"su2 postAnalysis")
        su2.postAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Change angle of attach
        if verbose: print(6*"-", line,"Modify Alpha")
        su2.input.Alpha = 10; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())
        
        # Run a 2nd time
        if verbose: print(6*"-", line,"su2 preAnalysis")
        su2.preAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"su2 system")
        su2.system("echo 'Running SU2 2'"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"su2 postAnalysis")
        su2.postAnalysis(); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # make sure the last call journals everything
        return line+2

#==============================================================================
    def test_journal(self):

        problemName = self.problemName+"2"
        
        myProblem = pyCAPS.Problem(problemName, capsFile=self.file, outLevel=0)

        # Run once to get the total line count
        line_total = self.run_journal(myProblem, -1)
        
        myProblem.close()
        shutil.rmtree(problemName)
        
        #print(80*"=")
        #print(80*"=")
        # Create the problem to start journaling
        myProblem = pyCAPS.Problem(problemName, capsFile=self.file, outLevel=0)
        myProblem.close()
        
        for line_exit in range(line_total):
            #print(80*"=")
            myProblem = pyCAPS.Problem(problemName, capsFile=self.file, outLevel=0, useJournal=True)
            self.run_journal(myProblem, line_exit)
            myProblem.close()

if __name__ == '__main__':
    unittest.main()
