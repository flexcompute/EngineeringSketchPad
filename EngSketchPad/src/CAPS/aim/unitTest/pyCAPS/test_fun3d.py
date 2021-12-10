import unittest

import os
import glob
import shutil
import __main__
import gc


# f90nml is used to write fun3d inputs not available in the aim
import f90nml

import pyCAPS

class TestFUN3D(unittest.TestCase):

    @classmethod
    def setUpClass(cls):

        cls.file = os.path.join("..","csmData","cfdMultiBody.csm")
        cls.problemName = "workDir_fun3dTest"
        cls.iProb = 1
        cls.cleanUp()

        cls.configFile = "fun3d.nml"

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        egadsTess = cls.myProblem.analysis.create(aim = "egadsTessAIM")

        # Modify local mesh sizing parameters
        egadsTess.input.Tess_Params = [0.5, 0.01, 15.0]
        Mesh_Sizing = {"Farfield": {"tessParams" : [7.,  2., 20.0]}}
        egadsTess.input.Mesh_Sizing = Mesh_Sizing

        tetgen = cls.myProblem.analysis.create(aim = "tetgenAIM",
                                               name = "tetgen")

        tetgen.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

        cls.myAnalysis = cls.myProblem.analysis.create(aim = "fun3dAIM")

        cls.myAnalysis.input["Mesh"].link(tetgen.output["Volume_Mesh"])

        cls.myAnalysis.input.Overwrite_NML = True

    @classmethod
    def tearDownClass(cls):
        del cls.myAnalysis
        del cls.myProblem
        cls.cleanUp()

    @classmethod
    def cleanUp(cls):

        # Remove problem directories
        dirs = glob.glob( cls.problemName + '*')
        for dir in dirs:
            if os.path.isdir(dir):
                shutil.rmtree(dir)


    # Try put an invalid boundary name
    def test_invalidBoundaryName(self):
        # Create a new instance
        # Use self.fun3d to ensure fun3d is not tracked in the traceback frames
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"}, # No capsGroup 'X'
                                               "X": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}
        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.fun3d.preAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        del self.fun3d

    # Try an invalid boundary type
    def test_invalidBoundary(self):
        # Create a new instance
        # Use self.fun3d to ensure fun3d is not tracked in the traceback frames
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "X"}}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.fun3d.preAnalysis()

        self.assertEqual(e.exception.errorName, "CAPS_NOTFOUND")

        del self.fun3d

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

    # Turn off Overwrite_NML
    def test_overwriteNML(self):

        # Create a new instance
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        self.fun3d.input.Overwrite_NML = False

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        self.fun3d.preAnalysis() # Don't except a config file because Overwrite_NML = False

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, self.configFile)), False)

        del self.fun3d

    # Create sensitvities
    def test_sensitivity(self):

        # Create a new instance
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        self.fun3d.input.Design_Functional = {"Composite": [{"function":"ClCd", "weight": 1.0, "target": 2.7, "power": 1.0, "capsGroup":"Wing1"},
                                                            {"function":"Cd"  , "weight": 1.0, "target": 0.0, "power": 2.0}],
                                              "Lift^2": {"function":"Cl", "weight": 1.0, "target": 0.0, "power": 2.0}}

        self.fun3d.input.Design_Sensitivity = False

        self.fun3d.input.Alpha = 1
        self.fun3d.input.Use_Python_NML = False
        self.fun3d.input.Overwrite_NML = True

        self.fun3d.preAnalysis()
        self.fun3d.postAnalysis()

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, os.path.join("Flow", self.configFile))), True)
        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, "rubber.data")), True)
        self.assertEqual(os.path.isdir(os.path.join(self.fun3d.analysisDir, "Rubberize")), True)
        self.assertEqual(os.path.isdir(os.path.join(self.fun3d.analysisDir, "Adjoint")), True)

        # Turn on the python nml writer, sensitvity, and remove the old files
        self.fun3d.input.Use_Python_NML = True
    
        self.fun3d.input.Design_Sensitivity = True
        self.fun3d.input.Design_Variable = {"Alpha":"",
                                            "Beta":{},
                                            "area":{},
                                            "aspect":{},
                                            "taper":"",
                                            "twist":"",
                                           }

        for file in os.listdir(self.fun3d.analysisDir):
            path = os.path.join(self.fun3d.analysisDir, file)
            try:
                shutil.rmtree(path)
            except OSError:
                os.remove(path)

        # Write the inputs again
        self.fun3d.preAnalysis()

        rubber = os.path.join(self.fun3d.analysisDir, "rubber.data")
        rubber_bak = os.path.join(self.fun3d.analysisDir, "rubber.data.bak")
        shutil.copy2(rubber, rubber_bak)

        writeValue = False
        varDeriv = 0
        varDerivMax = 0
        with open(rubber_bak, 'r') as bak:
            with open(rubber, 'w') as f:
                for line in bak:
                    if 'Current value of function' in line:
                        writeValue = True
                    elif writeValue:
                        line = "10.0\n"
                        writeValue = False
                    elif 'Current derivatives of function wrt global design variables' in line:
                        varDeriv = 1
                        varDerivMax = 6
                    elif 'Current derivatives of function wrt shape design variables' in line:
                        varDeriv += 1
                        varDerivMax = 8
                    elif varDeriv > 0:
                        line = str(float(varDeriv)) + '\n'
                        varDeriv += 1
                        if varDeriv > varDerivMax: varDeriv = 0

                    f.write(line)

        self.fun3d.postAnalysis()

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, os.path.join("Flow", self.configFile))), True)
        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, "rubber.data")), True)
        self.assertEqual(os.path.isdir(os.path.join(self.fun3d.analysisDir, "Rubberize")), True)
        self.assertEqual(os.path.isdir(os.path.join(self.fun3d.analysisDir, "Adjoint")), True)

        Objective_Value  = self.fun3d.dynout["Composite"].value
        Objective_Alpha  = self.fun3d.dynout["Composite"].deriv("Alpha")
        Objective_Beta   = self.fun3d.dynout["Composite"].deriv("Beta")
        Objective_area   = self.fun3d.dynout["Composite"].deriv("area")
        Objective_aspect = self.fun3d.dynout["Composite"].deriv("aspect")
        Objective_taper  = self.fun3d.dynout["Composite"].deriv("taper")
        Objective_twist  = self.fun3d.dynout["Composite"].deriv("twist")

        self.assertEqual(Objective_Value   , 10.0)
        self.assertEqual(Objective_Alpha   , 2.0)
        self.assertEqual(Objective_Beta    , 3.0)
        self.assertEqual(Objective_area    , 1.0 *3)
        self.assertEqual(Objective_aspect  , 2.0 *3)
        self.assertEqual(Objective_taper[0], 3.0 *3)
        self.assertEqual(Objective_taper[1], 4.0 *3)
        self.assertEqual(Objective_taper[2], 5.0 *3)
        self.assertEqual(Objective_taper[3], 6.0 *3)
        self.assertEqual(Objective_twist[0], 7.0 *3)
        self.assertEqual(Objective_twist[1], 8.0 *3)

        # Check without shape design parameters
        self.fun3d.input.Design_Variable = {"Alpha": {"upperBound": 10.0}}

        self.fun3d.preAnalysis()
        self.fun3d.postAnalysis()

        Objective_Value = self.fun3d.dynout["Composite"].value
        Objective_Alpha = self.fun3d.dynout["Composite"].deriv("Alpha")

        del self.fun3d

    # Test using Cython to write and modify the *.nml file
    def test_cythonNML(self):

        # Create a new instance
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        fun3dnml = f90nml.Namelist()
        fun3dnml['boundary_output_variables'] = f90nml.Namelist()
        fun3dnml['boundary_output_variables']['mach'] = True
        fun3dnml['boundary_output_variables']['cp'] = True
        fun3dnml['boundary_output_variables']['average_velocity'] = True

        fun3dnml.write(os.path.join(self.fun3d.analysisDir,self.configFile), force=True)

        self.fun3d.input.Use_Python_NML = True
        self.fun3d.input.Overwrite_NML = False # append

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous", "wallTemperature": 1},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        self.fun3d.preAnalysis()

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, self.configFile)), True)

        del self.fun3d

     # Test using Cython to write and modify the *.nml file with reentrance into the AIM
    def test_cythonNMLReentrance(self):

        # Create a new instance
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        self.fun3d.input["Mesh"].link(self.myProblem.analysis["tetgen"].output["Volume_Mesh"])

        self.fun3d.input.Use_Python_NML = True
        self.fun3d.input.Overwrite_NML = False # append

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Viscous"},
                                               "Wing2": {"bcType" : "Inviscid"},
                                               "Farfield":"farfield"}

        self.fun3d.preAnalysis()

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, self.configFile)), True)

        self.fun3d.postAnalysis()

        self.fun3d.input.Mach = 0.8

        self.fun3d.preAnalysis()

        self.assertEqual(os.path.isfile(os.path.join(self.fun3d.analysisDir, self.configFile)), True)

        self.fun3d.postAnalysis()

        del self.fun3d

    # Test using Cython to write and modify the *.nml file - catch an error
    def test_cythonNMLError(self):

        # Create a new instance
        # Use self.fun3d to ensure fun3d is not tracked in the traceback frames
        self.fun3d = self.myProblem.analysis.create(aim = "fun3dAIM")

        # Create a bad nml file
        f = open(os.path.join(self.fun3d.analysisDir,self.configFile), "w")
        f.write("&badNamelist")
        f.close()

        self.fun3d.input.Use_Python_NML = True
        self.fun3d.input.Overwrite_NML = False # append

        self.fun3d.input.Boundary_Condition = {"Wing1": {"bcType" : "Inviscid"}}

        with self.assertRaises(pyCAPS.CAPSError) as e:
            self.fun3d.preAnalysis()
        self.assertEqual(e.exception.errorName, "CAPS_BADVALUE")

        del self.fun3d

if __name__ == '__main__':
    unittest.main()
