# Import pyCAPS module
import pyCAPS

# Import other modules
import os
import argparse

# Prevent openmdao from importing mpi4py
# For OpenMPI, importing mpi4py causes this enviornment variable to be set:
# OMPI_MCA_ess=singleton
# which in turms prevents a system call with mpirun to function
# Using os.gettenv, os.unsetenv, os.setenv may be another workaround...
os.environ['OPENMDAO_USE_MPI'] = "0"

# Import OpenMDAO v3 module
import openmdao
import openmdao.api as om

if openmdao.__version__[:openmdao.__version__.find(".")] != "3":
    print( "This example is designed for OpenMDAO version 3")
    print( "Installed OpenMDAO :", openmdao.__version__ )
    exit()

class Fun3dAnalysis(om.ExplicitComponent):

    def initialize(self):
    # Setup variables so that arguments can be provided to the component.

        # CAPS Problem input
        self.options.declare('capsProblem', types=object)
        self.options.declare('cores', types=int)

    def setup(self):
    # Assign the actual values to the variables.

        capsProblem = self.options['capsProblem']
        fun3d = capsProblem.analysis['fun3d']

        # attach free design parameters to OpenMDAO fluid object
        self.add_input('Alpha' , val=fun3d.input.Alpha)
        self.add_input('Beta'  , val=fun3d.input.Beta)
        self.add_input('area'  , val=capsProblem.geometry.despmtr.area)
        self.add_input('aspect', val=capsProblem.geometry.despmtr.aspect)

        # Add outputs for both objective and constraints
        self.add_output('Cd2')
        self.add_output('Const')

    def compute(self, inputs, outputs):

        cores = self.options['cores']
        capsProblem = self.options['capsProblem']
        fun3d = capsProblem.analysis['fun3d']

        # Update input values
        fun3d.input.Alpha = inputs['Alpha']
        fun3d.input.Beta  = inputs['Beta']

        capsProblem.geometry.despmtr.area = inputs['area']
        capsProblem.geometry.despmtr.aspect = inputs['aspect']

        # Don't compute geometric sensitivities
        fun3d.input.Design_Sensitivity = False

        fun3d.preAnalysis()
        fun3d.system(f"mpirun -np {cores} nodet_mpi --design_run | tee nodet.out", 'Flow')
        fun3d.postAnalysis()

        # Grab functional and attach as an output
        outputs['Cd2']   = fun3d.dynout.Cd2
        outputs['Const'] = fun3d.dynout.Const

    def setup_partials(self):
    # declare the partial derivatives

        # Declare and attach partials to self
        self.declare_partials('Cd2','Alpha')
        self.declare_partials('Cd2','Beta')
        self.declare_partials('Cd2','area')
        self.declare_partials('Cd2','aspect')

        self.declare_partials('Const','Alpha')
        self.declare_partials('Const','Beta')
        self.declare_partials('Const','area')
        self.declare_partials('Const','aspect')

    def compute_partials(self,inputs,partials):
    # The actual value of the partial derivative is assigned here.

        cores = self.options['cores']
        capsProblem = self.options['capsProblem']
        fun3d = capsProblem.analysis['fun3d']

        # Request analysis and/or geometric design sensitivities for the Design_Functional
        fun3d.input.Design_Sensitivity = True

        fun3d.preAnalysis()
        fun3d.system(f"mpirun -np {cores} dual_mpi --getgrad --outer_loop_krylov | tee dual.out",'Adjoint')
        fun3d.postAnalysis()

        # Get derivatives and set partials
        partials['Cd2', 'Alpha']  = fun3d.dynout["Cd2"].deriv("Alpha")
        partials['Cd2', 'Beta']   = fun3d.dynout["Cd2"].deriv("Beta")
        partials['Cd2', 'area']   = fun3d.dynout["Cd2"].deriv("area")
        partials['Cd2', 'aspect'] = fun3d.dynout["Cd2"].deriv("aspect")

        partials['Const', 'Alpha']  = fun3d.dynout["Const"].deriv("Alpha")
        partials['Const', 'Beta']   = fun3d.dynout["Const"].deriv("Beta")
        partials['Const', 'area']   = fun3d.dynout["Const"].deriv("area")
        partials['Const', 'aspect'] = fun3d.dynout["Const"].deriv("aspect")

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Fun3D and OpenMDAO v3 PyTest Example',
                                 prog = 'fun3d_OpenMDAO_3_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
parser.add_argument("-cores", default = 1, type=int, help="Number of MPI cores")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "FUN3D_OpenMDAO_v3")

#######################################
# Create CAPS Problem                ##
#######################################

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

geometry = myProblem.geometry

#######################################
# Load AFLR4 aim                     ##
#######################################
aflr4 = myProblem.analysis.create(aim = "aflr4AIM", name = "aflr4")

# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Farfield growth factor
aflr4.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.5
aflr4.input.min_scale = 0.05

#######################################
# Load AFLR3 aim                     ##
#######################################
aflr3 = myProblem.analysis.create(aim = "aflr3AIM", name = "aflr3")

aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

# Set AIM verbosity
aflr3.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

#######################################
## Load FUN3D aim                    ##
#######################################

fun3d = myProblem.analysis.create(aim = "fun3dAIM", name = "fun3d")

fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])

# Set inputs
fun3d.input.Proj_Name = "fun3dOpenMDAO"
fun3d.input.Alpha = 1.0
fun3d.input.Beta = 1.0
fun3d.input.Viscous = "laminar"
fun3d.input.Mach = 0.5901
fun3d.input.Re = 10E5
fun3d.input.Equation_Type = "compressible"
fun3d.input.Num_Iter = 10
fun3d.input.CFL_Schedule = [0.5, 3.0]
fun3d.input.Restart_Read = "off"
fun3d.input.CFL_Schedule_Iter = [1, 40]
fun3d.input.Overwrite_NML = True

# Set boundary conditions
bc1 = {"bcType" : "Viscous", "wallTemperature" : 1}
bc2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
fun3d.input.Boundary_Condition = {"Wing1"   : bc1,
                                  "Wing2"   : bc2,
                                  "Farfield": "farfield"}

# Setup Cd**2 as the objective functional to minimize
# and Const = Cl + Cmz**2 for the constraint
fun3d.input.Design_Functional = {"Cd2": {"function":"Cd", "power":2},
                                 "Const": [{"function":"Cl"}, {"function":"Cmz","power":2}]}

# Declare design variables
fun3d.input.Design_Variable = {"Alpha":"",
                               "Beta":"",
                               "area":"",
                               "aspect":""}

#######################################
## Setup OpenMDAO                    ##
#######################################

# Setup the openmdao problem object
omProblem = om.Problem()

# Add the fluid analysis - all the options are added as arguments
fun3dSystem = Fun3dAnalysis(cores = args.cores,
                            capsProblem  = myProblem)

omProblem.model.add_subsystem('fun3dSystem', fun3dSystem)

# setup the optimization
omProblem.driver = om.ScipyOptimizeDriver()
omProblem.driver.options['optimizer'] = 'SLSQP'
#omProblem.driver.options['optimizer'] = "L-BFGS-B"
omProblem.driver.options['tol'] = 1.e-7
omProblem.driver.options['disp'] = True
omProblem.driver.options['maxiter'] = 3

# add design variables to model
area_limits   = geometry.despmtr["area"].limits
aspect_limits = geometry.despmtr["aspect"].limits

omProblem.model.add_design_var('fun3dSystem.Alpha' , lower=-25.0, upper=25.0)
omProblem.model.add_design_var('fun3dSystem.Beta'  , lower=-15.0, upper=15.0)
omProblem.model.add_design_var('fun3dSystem.area'  , lower=area_limits[0]  , upper=area_limits[1])
omProblem.model.add_design_var('fun3dSystem.aspect', lower=aspect_limits[0], upper=aspect_limits[1])

# minimize Cd**2 subject to (Cl + Cmz**2) == 0.1
omProblem.model.add_objective('fun3dSystem.Cd2')
omProblem.model.add_constraint('fun3dSystem.Const', equals = 0.1)

# Press go
omProblem.setup()
omProblem.run_driver()
omProblem.cleanup()
