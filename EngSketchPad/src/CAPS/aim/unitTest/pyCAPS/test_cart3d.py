import unittest

import os
import glob
import shutil
import __main__
import gc
import platform

import pyCAPS

# Helper function to check if an executable exists
def which(program):
    import os
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    fpath, fname = os.path.split(program)
    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None


class TestCart3D(unittest.TestCase):

    @classmethod
    def setUpClass(cls):
        if which("flowCart") == None:
            return

        cls.file = os.path.join("..","csmData","cfdMultiBody.csm")
        cls.problemName = "workDir_cart3dTest"
        cls.iProb = 1
        cls.cleanUp()

        cls.configFile = "Components.i.tri"

        cls.myProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        cls.myAnalysis = cls.myProblem.analysis.create(aim = "cart3dAIM")

    @classmethod
    def tearDownClass(cls):
        if which("flowCart") == None:
            return

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

    # This is executed prior to each test
    def setUp(self):
        if platform.system() == "Windows":
            self.skipTest("Cart3D does not exist for Windows")
        else:
            if which("flowCart") == None:
                self.skipTest("No flowCart executable")

    # Test available outputs
    def test_outputs(self):

        keys = list(self.myAnalysis.output.keys())
        self.assertEqual(sorted(keys), sorted(["C_A", "C_Y", "C_N",   "C_D",   "C_S", "C_L", "C_l",
                                               "C_m", "C_n", "C_M_x", "C_M_y", "C_M_z"])  )
    # Test re-enter
    def test_reenter(self):

        self.myAnalysis.input.Mach = 0.5
        self.myAnalysis.input.outer_box = 80

        self.myAnalysis.input.maxCycles = 1
        self.myAnalysis.input.nOrders = 1

        self.myAnalysis.input.aerocsh = ["set it_fc = 20", # Number of fine-grid flowCart iterations on initial mesh
                                         "set it_ad = 10", # Number of fine-grid adjointCart iterations on each mesh
                                         "set mg_init = 1",# Number of multigrid levels for initial mesh (default 2, ramps up to mg_fc/ad)
                                         "set mg_fc = 1",  # Number of flowCart multigrid levels (default=3)
                                         "set mg_ad = 1",  # Number of adjointCart multigrid levels (usually same as flowCart)
                                        ]

        self.myAnalysis.input.Restart = True

        self.myAnalysis.input.Adapt_Functional = {"Drag": {"function":"C_D"}}

        # Generate Components.i.tri
        self.myAnalysis.input.Mach = 0.5
        self.myAnalysis.runAnalysis()

        # Try exatracting values for each output
        for name in self.myAnalysis.output:
            val = self.myAnalysis.output[name].value
            print(name, val)
            #self.assertNotEqual(val, 0)

        # Doe not change grid, so no new Components.i.tri
        self.myAnalysis.input.Mach = 0.4
        self.myAnalysis.runAnalysis()

        # Try exatracting values for each output
        for name in self.myAnalysis.output:
            val = self.myAnalysis.output[name].value
            print(name, val)
            #self.assertNotEqual(val, 0)

        # Changes surface grid, makes new Components.i.tri
        self.myAnalysis.input.Tess_Params = self.myAnalysis.input.Tess_Params
        self.myAnalysis.runAnalysis()

        # Try exatracting values for each output
        for name in self.myAnalysis.output:
            val = self.myAnalysis.output[name].value
            print(name, val)
            #self.assertNotEqual(val, 0)

    # Create sensitvities
    def test_sensitivity_AnalysisIn(self):

        # Create a new instance
        self.cart3d = self.myProblem.analysis.create(aim = "cart3dAIM", autoExec=True)

        self.cart3d.input.Mach = 0.5
        self.cart3d.input.outer_box = 80

        self.cart3d.input.maxCycles = 1
        self.cart3d.input.nOrders = 4

        self.cart3d.input.alpha = 1

        self.cart3d.input.aerocsh = ["set it_fc = 10", # Number of fine-grid flowCart iterations on initial mesh
                                     "set it_ad = 10", # Number of fine-grid adjointCart iterations on each mesh
                                     "set mg_init = 1",# Number of multigrid levels for initial mesh (default 2, ramps up to mg_fc/ad)
                                     "set mg_fc = 1",  # Number of flowCart multigrid levels (default=3)
                                     "set mg_ad = 1",  # Number of adjointCart multigrid levels (usually same as flowCart)
                                    ]

        self.cart3d.input.Design_Variable = {"Mach":""}

        funs = [["C_A"   , "C_Y"   , "C_N"],
                ["C_D"   , "C_S"   , "C_L"],
                ["C_l"]  ,["C_m"]  ,["C_n"],
                ["C_M_x"],["C_M_y"], ["C_M_z"]]

        for names in funs:
            # Create a design functional for each regular output function
            Design_Functional = {}

            for name in names: #self.cart3d.output:
                Design_Functional[name] = {"function":name, 'capsGroup':"Wing1"}
            #Design_Functional["LoD"] = {"function":"LoD"}

            self.cart3d.input.Design_Functional = Design_Functional
            #self.cart3d.input.Design_Adapt = "C_D"

            self.cart3d.input.Design_Sensitivity = False
            self.cart3d.runAnalysis()

            # Try exatracting values for each output amd checm the dynamic outputs match
            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = True
            self.cart3d.runAnalysis()

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = False
            self.cart3d.runAnalysis()

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

        del self.cart3d


    # Create sensitvities
    def test_sensitivity_GeometryIn(self):

        # Create a new instance
        self.cart3d = self.myProblem.analysis.create(aim = "cart3dAIM", autoExec=True)

        self.cart3d.input.Mach = 0.5
        self.cart3d.input.outer_box = 80

        self.cart3d.input.maxCycles = 1
        self.cart3d.input.nOrders = 4

        self.cart3d.input.alpha = 1

        self.cart3d.input.aerocsh = ["set it_fc = 10", # Number of fine-grid flowCart iterations on initial mesh
                                     "set it_ad = 10", # Number of fine-grid adjointCart iterations on each mesh
                                     "set mg_init = 1",# Number of multigrid levels for initial mesh (default 2, ramps up to mg_fc/ad)
                                     "set mg_fc = 1",  # Number of flowCart multigrid levels (default=3)
                                     "set mg_ad = 1",  # Number of adjointCart multigrid levels (usually same as flowCart)
                                    ]

        self.cart3d.input.Design_Variable = {"area":"", "Mach":""}

        funs = [["C_A"   , "C_Y"   , "C_N"],
                #["C_D"   , "C_S"   , "C_L"],
                #["C_l"]  ,["C_m"]  ,["C_n"],
                #["C_M_x"],["C_M_y"], ["C_M_z"]]
               ]

        for names in funs:
            # Create a design functional for each regular output function
            Design_Functional = {}

            for name in names: #self.cart3d.output:
                Design_Functional[name] = {"function":name, 'capsGroup':"Wing1"}
            #Design_Functional["LoD"] = {"function":"LoD"}

            self.cart3d.input.Design_Functional = Design_Functional
            #self.cart3d.input.Design_Adapt = "C_D"

            self.cart3d.input.Design_Sensitivity = False
            self.cart3d.runAnalysis()

            # Try exatracting values for each output amd checm the dynamic outputs match
            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = True
            self.cart3d.runAnalysis()

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = False
            self.cart3d.runAnalysis()

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

        del self.cart3d

if __name__ == '__main__':
    unittest.main()
