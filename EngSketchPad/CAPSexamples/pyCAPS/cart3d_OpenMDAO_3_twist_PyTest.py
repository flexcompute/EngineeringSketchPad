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

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Cart3D and OpenMDAO v3 PyTest Example',
                                 prog = 'cart3d_OpenMDAO_3_twist_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
parser.add_argument("-plotFunctional", action='store_true', default = False, help="Plot the functional and gradient")
args = parser.parse_args()


class Cart3dAnalysis(om.ExplicitComponent):

    def initialize(self):
    # Setup variables so that arguments can be provided to the component.

        # CAPS Problem input
        self.options.declare('capsProblem', types=object)

    def setup(self):
    # Assign the actual values to the variables.

        capsProblem = self.options['capsProblem']
        geometry = capsProblem.geometry

        # attach free design parameters to OpenMDAO fluid object
        self.add_input('twist' , val=geometry.despmtr.twist)

        # Add outputs for both objective and constraints
        self.add_output('CL2')

    def compute(self, inputs, outputs):

        capsProblem = self.options['capsProblem']
        geometry = capsProblem.geometry
        cart3d = capsProblem.analysis['cart3d']

        # Update input values if changed
        if geometry.despmtr.twist != inputs['twist']:
            geometry.despmtr.twist = inputs['twist']

            # Don't compute geometric sensitivities
            cart3d.input.Design_Sensitivity = False

        print('-'*40)
        print('>>> twist', inputs['twist'])

        # Grab functional and attach as an output
        outputs['CL2'] = cart3d.dynout["CL2"].value
        print('>>> CL2', outputs['CL2'], " CL", cart3d.output.C_L)
        print('-'*40)

    def setup_partials(self):
    # declare the partial derivatives

        # Declare and attach partials to self
        self.declare_partials('CL2','twist')

    def compute_partials(self,inputs,partials):
    # The actual value of the partial derivative is assigned here.

        capsProblem = self.options['capsProblem']
        geometry = capsProblem.geometry
        cart3d = capsProblem.analysis['cart3d']

        # Update input values if changed
        if geometry.despmtr.twist != inputs['twist']:
            geometry.despmtr.twist = inputs['twist']

        # Request analysis and/or geometric design sensitivities for the Design_Functional
        cart3d.input.Design_Sensitivity = True

        # Get derivatives and set partials
        print('='*40)
        partials['CL2', 'twist']  = cart3d.dynout["CL2"].deriv("twist")
        print('>>> CL2_twist', partials['CL2', 'twist'])
        print('='*40)

# Working directory
workDir = os.path.join(str(args.workDir[0]), "Cart3D_OpenMDAO_v3_twist")

#######################################
# Create CAPS Problem                ##
#######################################

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","cfd_airfoilSection.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

geometry = myProblem.geometry

geometry.despmtr.taper = 1.0
geometry.despmtr.twist = 4.0
geometry.cfgpmtr.aerosystem = True

#######################################
## Load Cart3D aim                   ##
#######################################

cart3d = myProblem.analysis.create(aim = "cart3dAIM", name = "cart3d")

# Set inputs
cart3d.input.alpha = 0.0
cart3d.input.beta = 0.0
cart3d.input.Mach = 0.5901

cart3d.input.maxR = 7
cart3d.input.preSpec = True
cart3d.input.outer_box = 5.12345
cart3d.input.Tess_Params = [0.005, 0.00025, 10.00]

cart3d.input.y_is_spanwise = True
cart3d.input.symmY = True
cart3d.input.halfBody = True

cart3d.input.Model_X_axis = "-Xb"
cart3d.input.Model_Y_axis = "Yb"
cart3d.input.Model_Z_axis = "-Zb"

# Setup C_L**2 as the objective functional to minimize
cart3d.input.Design_Functional = {"CL2": {"function":"C_L", "power":2}}

# Declare design variables
cart3d.input.Design_Variable = {"twist":""}

cart3d.input.aerocsh = ["set it_fc = 150",
                        "set it_ad = 150",
                        #"set cubes_a = 3    # angle criterion (-a)",
                        #"set cubes_b = 5    # buffer layers   (-b)",
                        #"set cflmin = 0.2"
                        ]

# Plot the functional and gradient if requested
if args.plotFunctional:
    try:
        import matplotlib.pyplot as plt
        import numpy as np

        # Turn on sensitivity calculations
        cart3d.input.Design_Sensitivity = True

        twists = np.linspace(-5, 5, 21)
        CL2 = []
        CL2_twist = []
        for twist in twists:
            print("-"*20)
            print("twist", twist)
            print("-"*20)
            geometry.despmtr.twist = twist
        
            CL2.append(cart3d.dynout["CL2"].value)
            CL2_twist.append(cart3d.dynout["CL2"].deriv("twist"))
        
            # Create an archive of cart3d files
            #os.system("cd " + cart3d.analysisDir + " && c3d_objGrad.csh -clean")
            #os.rename(os.path.join(cart3d.analysisDir, "design"),  os.path.join(cart3d.analysisDir, "design_G"+str(twist)))

        # Compute central difference derivatives
        dCL2_twist = []
        for i in range(1,len(twists)-1):
            dCL2_twist.append( (CL2[i+1]-CL2[i-1])/(twists[i+1]-twists[i-1]) )

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        
        ax2.axhline(0., linewidth=1, linestyle=":", color='red')

        # Plot the functional
        lns1 = ax1.plot(twists, CL2, 's-', label = r"$C_L^2$", color="blue")
        
        # Plot plot the derivative
        lns2 = ax2.plot(twists, CL2_twist, 'o-', label = r"$\partial C_L^2/\partial \gamma$", color="red")

        # Plot plot the FD derivative
        lns3 = ax2.plot(twists[1:-1], dCL2_twist, '.--', label = r"$\Delta C_L^2/\Delta \gamma$", color='red', mec='k')

        ax1.set_xlabel(r"$\gamma$ $(^o)$")
        ax1.set_ylabel(r"$C_L^2$", color="blue")
        ax2.set_ylabel(r"$\partial C_L^2/\partial \gamma$", color="red")
        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        
        # Shrink current axis's height by 15% on the top
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width, box.height * 0.85])
        ax2.set_position([box.x0, box.y0, box.width, box.height * 0.85])
        ax1.set_title("(window must be closed to continue Python script)", y=1.0, pad=34)

        # add legend
        lns = lns1+lns2+lns3
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc='lower center', bbox_to_anchor=(0.5, 1.01), fancybox=True, ncol=3)

        #plt.savefig("cart3d_twist.png")
        plt.show()

        geometry.despmtr.twist = [1.0, 2.0]

    except ImportError:
        print ("Unable to import matplotlib.pyplot module. Plotting is not available")
        plt = None


#######################################
## Setup OpenMDAO                    ##
#######################################

# Setup the openmdao problem object
omProblem = om.Problem()

# Add the fluid analysis - all the options are added as arguments
cart3dSystem = Cart3dAnalysis(capsProblem  = myProblem)

omProblem.model.add_subsystem('cart3dSystem', cart3dSystem)

# setup the optimization
omProblem.driver = om.ScipyOptimizeDriver()
#omProblem.driver.options['optimizer'] = 'SLSQP'
omProblem.driver.options['optimizer'] = "L-BFGS-B"
omProblem.driver.options['tol'] = 1.e-6
omProblem.driver.options['disp'] = True
omProblem.driver.options['maxiter'] = 10

# add design variables to model
omProblem.model.add_design_var('cart3dSystem.twist', lower=-25.0, upper=25.0)

# minimize Cl**2
omProblem.model.add_objective('cart3dSystem.CL2')

# Press go
omProblem.setup()
omProblem.run_driver()
omProblem.cleanup()

opt_twist = omProblem.get_val("cart3dSystem.twist")
print("Optimized values:", opt_twist)
assert abs(opt_twist) < 0.2
