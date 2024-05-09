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

        cls.myProblem.geometry.cfgpmtr.wake = 0

        cls.myGeometry = cls.myProblem.geometry

        aflr4 = cls.myProblem.analysis.create(aim = "aflr4AIM",
                                              name = "aflr4")

        # Modify local mesh sizing parameters
        aflr4.input.Mesh_Quiet_Flag = True
        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr = 1.1

        aflr3 = cls.myProblem.analysis.create(aim = "aflr3AIM",
                                              name = "aflr3")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        cls.myAnalysis = cls.myProblem.analysis.create(aim = "su2AIM")

        cls.myAnalysis.input["Mesh"].link(aflr3.output["Volume_Mesh"])

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

        self.aflr3 = self.myProblem.analysis["aflr3"]

        self.su2.input["Mesh"].link(self.aflr3.output["Volume_Mesh"])

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

        self.aflr3 = self.myProblem.analysis["aflr3"]

        self.su2.input["Mesh"].link(self.aflr3.output["Volume_Mesh"])

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

        mm = pyCAPS.Unit("mm")
        m  = pyCAPS.Unit("meter")
        kg = pyCAPS.Unit("kg")
        s  = pyCAPS.Unit("s")
        K  = pyCAPS.Unit("K")
        kPa  = pyCAPS.Unit("kPa")

        unitSystem={"mass" : kg, "length": m, "time":s, "temperature":K}

        myAnalysis = self.myProblem.analysis.create(aim = "su2AIM", unitSystem=unitSystem)

        myAnalysis.input["Mesh"].link(self.myProblem.analysis["aflr3"].output["Volume_Mesh"])

        myAnalysis.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        myAnalysis.input.Freestream_Density  = 1.225*kg/m**3
        myAnalysis.input.Freestream_Velocity = 200*m/s
        myAnalysis.input.Freestream_Pressure = 101.325*kPa

        myAnalysis.input.Moment_Length = 20 * mm
        myAnalysis.input.Reference_Area = 100 * mm**2

        myAnalysis.preAnalysis()

#==============================================================================
    def test_phase(self):
        file = os.path.join("..","csmData","cfdMultiBody.csm")

        problemName = self.problemName + "_Phase"
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase0", capsFile=file, outLevel=0)

        myProblem.geometry.cfgpmtr.wake = 0

        aflr4 = myProblem.analysis.create(aim = "aflr4AIM",
                                          name = "aflr4")

        aflr4.input.Mesh_Quiet_Flag = True
        aflr4.input.Mesh_Length_Factor = 2
        aflr4.input.min_scale = 0.1
        aflr4.input.ff_cdfr = 1.1

        aflr3 = myProblem.analysis.create(aim = "aflr3AIM",
                                          name = "aflr3")

        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

        aflr3.input.Mesh_Quiet_Flag = True

        su2 = myProblem.analysis.create(aim = "su2AIM",
                                        name = "su2")

        su2.input["Mesh"].link(aflr3.output["Volume_Mesh"])

        su2.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"},
                                        "Wing2": {"bcType" : "Inviscid"},
                                        "Farfield":"farfield"}
        su2.preAnalysis()
        su2.postAnalysis()

        myProblem.closePhase()

        # Initialize Problem from the last phase and make a new phase
        myProblem = pyCAPS.Problem(problemName, phaseName="Phase1", phaseStart="Phase0", outLevel=0)

        aflr4 = myProblem.analysis["aflr4"]
        aflr3 = myProblem.analysis["aflr3"]
        su2   = myProblem.analysis["su2"]

        su2.input.Alpha = 5

        su2.preAnalysis()
        su2.postAnalysis()

#==============================================================================
    def run_journal(self, myProblem, line_exit):

        verbose = False

        line = 0
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Dissable wake")
        myProblem.geometry.cfgpmtr.wake = 0; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())


        # Load aflr4 AIM
        if verbose: print(6*"-", line,"Load aflr4AIM")
        aflr4 = myProblem.analysis.create(aim = "aflr4AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        aflr4.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Length_Factor")
        aflr4.input.Mesh_Length_Factor = 2; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify min_scale")
        aflr4.input.min_scale = 0.1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify ff_cdfr")
        aflr4.input.ff_cdfr = 1.1; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())


        # Create the aflr3 AIM
        if verbose: print(6*"-", line,"Load aflr3AIM")
        aflr3 = myProblem.analysis.create(aim = "aflr3AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Modify Mesh_Quiet_Flag")
        aflr3.input.Mesh_Quiet_Flag = True; line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Link the surface mesh
        if verbose: print(6*"-", line,"Link Surface_Mesh")
        aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"]); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        # Load SU2 AIM
        if verbose: print(6*"-", line,"Load su2AIM")
        su2 = myProblem.analysis.create(aim = "su2AIM"); line += 1
        if line == line_exit: return line
        if line_exit > 0: self.assertTrue(myProblem.journaling())

        if verbose: print(6*"-", line,"Link Mesh")
        su2.input["Mesh"].link(aflr3.output["Volume_Mesh"]); line += 1
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
            myProblem = pyCAPS.Problem(problemName, phaseName="Scratch", capsFile=self.file, outLevel=0)
            self.run_journal(myProblem, line_exit)
            myProblem.close()

if __name__ == '__main__':
    unittest.main()
