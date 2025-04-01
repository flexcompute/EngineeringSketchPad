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

        cls.capsProblem = pyCAPS.Problem(cls.problemName, capsFile=cls.file, outLevel=0)

        cls.capsProblem.geometry.cfgpmtr.series     = 12
        cls.capsProblem.geometry.cfgpmtr.series2    = 12
        cls.capsProblem.geometry.cfgpmtr.sharpte    = 1
        cls.capsProblem.geometry.despmtr.twist      = [0, 0]
        cls.capsProblem.geometry.despmtr.dihedral   = 0
        cls.capsProblem.geometry.cfgpmtr.aerosystem = 0

        cls.cart3d = cls.capsProblem.analysis.create(aim = "cart3dAIM")

    @classmethod
    def tearDownClass(cls):
        if which("flowCart") == None:
            return

        del cls.cart3d
        del cls.capsProblem
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

#==============================================================================
    def test_setInput(self):

        cart3d = self.capsProblem.analysis.create(aim = "cart3dAIM")

        cart3d.input.Tess_Params = [0.025,0.001,15]
        cart3d.input.mesh2d = False
        cart3d.input.outer_box = 30.
        cart3d.input.nDiv = 5
        cart3d.input.maxR = 7
        cart3d.input.preSpec = False
        cart3d.input.symmX = False
        cart3d.input.symmY = False
        cart3d.input.symmZ = False
        cart3d.input.halfBody = False
        
        # inlet boundary conditions       
        # Cart3d BC opt. 1: "SurfBC 1  1.0 0.3 0. 0. 0.714285  # compID   Rho    Xvel     Yvel    Zvel   Press"
        inletbc = {"bcType" : "Riemann", "staticDensity" : 1.0, "uVelocity": 0.3, "vVelocity":0, "wVelocity":0, "staticPressure":0.714285}
        
        # exhaust boundary condition
        # Cart3d BC opt. 1: "SurfBC 1  2.0 1.5 0. 0. 1.5  # compID   Rho    Xvel     Yvel    Zvel   Press"
        exhaustbc = {"bcType" : "Riemann", "staticDensity" : 2.0, "uVelocity": 1.5, "vVelocity":0, "wVelocity":0, "staticPressure":1.5}
        
        # Apply the selected BCs
        cart3d.input.Boundary_Condition = {"inlet"   : inletbc,
                                           "exhaust" : exhaustbc}

        cart3d.input.Mach = 0.76
        cart3d.input.alpha = 0.
        cart3d.input.beta = 0.
        cart3d.input.gamma = 1.4
        cart3d.input.Pressure_Scale_Factor = 1.0

        cart3d.input.TargetCL = 0.1
        cart3d.input.TargetCL_Tol = 0.01
        cart3d.input.TargetCL_Start_Iter = 60
        cart3d.input.TargetCL_Freq = 5
        cart3d.input.TargetCL_NominalAlphaStep = 0.2

        cart3d.input.maxCycles = 1000
        cart3d.input.nMultiGridLevels = 1
        cart3d.input.MultiGridCycleType = 2
        cart3d.input.MultiGridPreSmoothing = 1
        cart3d.input.MultiGridPostSmoothing = 1

        cart3d.input.RK = [[0.0695, 1],[0.1602,0], [0.2898,0], [0.5060,0], [1.0,0]]
        cart3d.input.RK = [[0.0695, 1],[0.1602,0], [0.2898,0], [1.0,0]]

        cart3d.input.CFL = 1.2
        cart3d.input.Limiter = 2
        cart3d.input.FluxFun = 0
        cart3d.input.iForce = 1
        cart3d.input.iHist = 1
        cart3d.input.nOrders = 1
        cart3d.input.nAdaptCycles = 1
        
        cart3d.input.nAdaptCycles = 5
        cart3d.input.Adapt_Functional = {"CD" : {"function" : "C_D", "power" : 1}}

        self.cart3d.input.Design_Variable = {"Mach":""}
        self.cart3d.input.Design_Functional =  {"CL" : {"function":"C_L", 'capsGroup':"Wing1"}}
        self.cart3d.input.Design_Sensitivity = False
        self.cart3d.input.Design_Adapt = "CL"
        self.cart3d.input.Design_Run_Config = "production"
        self.cart3d.input.Design_Gradient_Memory_Budget = 32
        
        self.cart3d.input.Xslices = [0, 0.5, 1.0]
        self.cart3d.input.Yslices = [0, 0.5, 1.0]
        self.cart3d.input.Zslices = [0, 0.5, 1.0]
        
        self.cart3d.input.y_is_spanwise = False

        self.cart3d.input.Model_X_axis = "-Xb"
        self.cart3d.input.Model_Y_axis = "Yb"
        self.cart3d.input.Model_Z_axis = "-Zb"

        self.cart3d.input.Restart = False

        self.cart3d.input.aerocsh = ["set it_fc = 10", # Number of fine-grid flowCart iterations on initial mesh
                                     "set it_ad = 10", # Number of fine-grid adjointCart iterations on each mesh
                                    ]
        
#==============================================================================
    # Test available outputs
    def test_outputs(self):

        keys = list(self.cart3d.output.keys())
        self.assertEqual(sorted(keys), sorted(["C_A", "C_Y", "C_N",   "C_D",   "C_S", "C_L", "C_l",
                                               "C_m", "C_n", "C_M_x", "C_M_y", "C_M_z", "alpha"])  )

#==============================================================================
    # Test re-enter
    def test_reenter(self):

        self.cart3d.input.Mach = 0.5
        self.cart3d.input.outer_box = 80

        self.cart3d.input.maxCycles = 1
        self.cart3d.input.nOrders = 1

        self.cart3d.input.aerocsh = ["set it_fc = 20", # Number of fine-grid flowCart iterations on initial mesh
                                     "set it_ad = 10", # Number of fine-grid adjointCart iterations on each mesh
                                     "set mg_init = 1",# Number of multigrid levels for initial mesh (default 2, ramps up to mg_fc/ad)
                                     "set mg_fc = 1",  # Number of flowCart multigrid levels (default=3)
                                     "set mg_ad = 1",  # Number of adjointCart multigrid levels (usually same as flowCart)
                                    ]

        self.cart3d.input.Restart = True

        self.cart3d.input.nAdaptCycles = 0
        self.cart3d.input.Adapt_Functional = {"Drag": {"function":"C_D"}}

        # Generate Components.i.tri
        self.cart3d.input.Mach = 0.5

        vals = {}

        # Try exatracting values for each output
        for name in self.cart3d.output:
            vals[name] = self.cart3d.output[name].value
            #print(name, val)
            #self.assertNotEqual(val, 0)

        # Doe not change grid, so no new Components.i.tri
        self.cart3d.input.Mach = 0.4

        # Try exatracting values for each output agan and make sure they have changed
        for name in self.cart3d.output:
            if name == "alpha": continue
            val = self.cart3d.output[name].value
            #print(name, val, vals[name])
            self.assertNotEqual(val, vals[name])
            vals[name] = val

        # Changes surface grid, makes new Components.i.tri
        Tess_Params = self.cart3d.input.Tess_Params
        Tess_Params[0] *= 0.8
        self.cart3d.input.Tess_Params = Tess_Params

        # Try exatracting values for each output and make sure they have changed
        for name in self.cart3d.output:
            if name == "alpha": continue
            val = self.cart3d.output[name].value
            #print(name, val)
            self.assertNotEqual(val, vals[name])

#==============================================================================
    # Test maching a target CL
    def test_TargetCL(self):

        # Create a new instance
        self.cart3d = self.capsProblem.analysis.create(aim = "cart3dAIM", autoExec=True)

        TargetCL = 0.1

        self.cart3d.input.mesh2d = True
        self.cart3d.input.nDiv = 9
        self.cart3d.input.maxR = 7

        self.cart3d.input.Mach = 0.2
        self.cart3d.input.outer_box = 10*self.capsProblem.geometry.outpmtr.cmean

        self.cart3d.input.Model_X_axis = "-Xb"
        self.cart3d.input.Model_Y_axis = "-Zb"
        self.cart3d.input.Model_Z_axis = "-Yb"
        #self.cart3d.input.Yslices = [0]
        self.cart3d.input.Zslices = [0]

        self.cart3d.input.alpha = 1
        self.cart3d.input.TargetCL = TargetCL

        # print(self.cart3d.output.C_L)
        # print(self.cart3d.output.C_S)
        # print(self.cart3d.output.C_D)
        # print(self.cart3d.output.alpha)
        self.assertAlmostEqual( TargetCL, self.cart3d.output.C_L, 2 )
        #self.assertAlmostEqual( 0.45018301428, self.cart3d.output.alpha, 2 )

        del self.cart3d

#==============================================================================
    # Test analysis sensitvities
    def test_sensitivity_AnalysisIn(self):

        # Create a new instance
        self.cart3d = self.capsProblem.analysis.create(aim = "cart3dAIM", autoExec=True)

        self.cart3d.input.Model_X_axis = "-Xb"
        self.cart3d.input.Model_Y_axis = "-Zb"
        self.cart3d.input.Model_Z_axis = "-Yb"

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

            # Try exatracting values for each output and check the dynamic outputs match
            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = True

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)
                self.assertNotEqual(0.0, self.cart3d.dynout[name].deriv("Mach"))

            self.cart3d.input.Design_Sensitivity = False

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

        del self.cart3d


#==============================================================================
    # Test geometry sensitvities
    def test_sensitivity_GeometryIn(self):

        # Create a new instance
        self.cart3d = self.capsProblem.analysis.create(aim = "cart3dAIM", autoExec=True)

        self.cart3d.input.Model_X_axis = "-Xb"
        self.cart3d.input.Model_Y_axis = "-Zb"
        self.cart3d.input.Model_Z_axis = "-Yb"

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

            # Try exatracting values for each output and check the dynamic outputs match
            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

            self.cart3d.input.Design_Sensitivity = True

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)
                self.assertNotEqual(0.0, self.cart3d.dynout[name].deriv("Mach"))
                self.assertNotEqual(0.0, self.cart3d.dynout[name].deriv("area"))

            self.cart3d.input.Design_Sensitivity = False

            for name in names: #self.cart3d.output:
                valOut = self.cart3d.output[name].value
                valDyn = self.cart3d.dynout[name].value
                #print(name, valOut, valDyn)
                self.assertAlmostEqual(valOut, valDyn, 1)

        del self.cart3d

if __name__ == '__main__':
    unittest.main()
