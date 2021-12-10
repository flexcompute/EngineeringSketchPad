# Import pyCAPS module
import pyCAPS

# Import other modules
import os
import argparse

# Import OpenMDAO v3 module
import openmdao
import openmdao.api as om

if openmdao.__version__[:openmdao.__version__.find(".")] != "3":
    print( "This example is designed for OpenMDAO version 3")
    print( "Installed OpenMDAO :", openmdao.__version__ )
    exit()

class avlAnalysis(om.ExplicitComponent):
    
    def initialize(self):
    # Setup variables so that arguments can be provided to the component.
        
        # CAPS Problem input
        self.options.declare('capsProblem', types=object)

    def setup(self):
    # Assign the actual values to the variables.
 
        capsProblem = self.options['capsProblem']
        avl = capsProblem.analysis['avl']

        # attach parameters to OpenMDAO fluid object
        self.add_input('Alpha' , val=avl.input.Alpha)
        self.add_input('Beta'  , val=avl.input.Beta)
        
        # Add output metric
        self.add_output('Cmtot2')

    def compute(self, inputs, outputs):

        capsProblem = self.options['capsProblem']
        avl = capsProblem.analysis['avl']

        # Update input values
        avl.input.Alpha = inputs['Alpha']
        avl.input.Beta  = inputs['Beta']
        
        # Grab objective and attach as an output
        outputs['Cmtot2']  = avl.output.Cmtot**2
        print("Alpha", inputs['Alpha'], "Beta", inputs['Beta'], "Cmtot2", outputs['Cmtot2'])

    def setup_partials(self):
    # declare the partial derivatives
        
        # Declare and attach partials to self
        self.declare_partials('Cmtot2','Alpha')
        self.declare_partials('Cmtot2','Beta')

    def compute_partials(self,inputs,partials):
    # The actual value of the partial derivative is assigned here. 
        # Assume input value are consisten with last call to compute
        
        capsProblem = self.options['capsProblem']
        avl = capsProblem.analysis['avl']

        # Get derivatives and set partials
        Cm  = avl.output["Cmtot"].value
        Cma = avl.output["Cmtot"].deriv("Alpha")
        Cmb = avl.output["Cmtot"].deriv("Beta")

        partials['Cmtot2', 'Alpha'] = 2*Cm*Cma
        partials['Cmtot2', 'Beta']  = 2*Cm*Cmb
        print("Cmtot2_Alpha", partials['Cmtot2', 'Alpha'])
        print("Cmtot2_Beta", partials['Cmtot2', 'Beta'])


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Avl and OpenMDAO v3 PyTest Example',
                                 prog = 'avl_OpenMDAO_3_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
parser.add_argument("-plotFunctional", action='store_true', default = False, help="Plot the functional and gradient")
args = parser.parse_args()
 
#######################################
# Create CAPS Problem                ##
#######################################

# Initiate problem and object
problemName = os.path.join(args.workDir[0], "avlOpenMDAO_v3")
capsProblem = pyCAPS.Problem(problemName="avlOpenMDAO_v3",
                             capsFile=os.path.join("..","csmData","avlWing.csm"),
                             outLevel=args.outLevel)

# Increase dihedral so Beta has more impact on moment
capsProblem.geometry.despmtr.dihedral = 15

# Setup AVL
avl = capsProblem.analysis.create(aim = "avlAIM", name="avl")

# Setup AVL surfaces
wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
        "numChord"          : 11,
        "spaceChord"        : 1.0,
        "numSpanPerSection" : 12,
        "spaceSpan"         : 1.0}

avl.input.AVL_Surface = {"Wing": wing}

# Set initial value angles
avl.input.Alpha = 2.
avl.input.Beta  = 3.

# Plot the functional and gradient if requested
if args.plotFunctional:
    try:
        import numpy as npy
        import matplotlib.pyplot as plt
        
        # Set Beta to 0
        avl.input.Beta  = 0.0

        Alphas = npy.linspace(-10, 10, 15)
        Cmtot2 = []
        Cmtot2_alpha = []
        for Alpha in Alphas:
            avl.input.Alpha = Alpha
        
            Cm  = avl.output["Cmtot"].value
            Cma = avl.output["Cmtot"].deriv("Alpha")

            Cmtot2.append(Cm**2)
            Cmtot2_alpha.append(2*Cm*Cma)

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()

        # Plot the functional
        lns1 = ax1.plot(Alphas, Cmtot2, 'o--', label = r"$C_M^2(\alpha)$", color='blue')
        
        # Plot plot the derivative
        lns2 = ax2.plot(Alphas, Cmtot2_alpha, 's-', label = r"$\partial C_M^2/\partial \alpha$", color='red')

        avl.input.Alpha  = -5.

        Betas = npy.linspace(-10, 10, 15)
        Cmtot2 = []
        Cmtot2_beta = []
        for Beta in Betas:

            avl.input.Beta = Beta
        
            Cm  = avl.output["Cmtot"].value
            Cmb = avl.output["Cmtot"].deriv("Beta")

            Cmtot2.append(Cm**2)
            Cmtot2_beta.append(2*Cm*Cmb)

        # Plot the functional
        lns3 = ax1.plot(Betas, Cmtot2, 'o-', label = r"$C_M^2(\beta)$", color='blue')
        
        # Plot plot the derivative
        lns4 = ax2.plot(Betas, Cmtot2_beta, 's-.', label = r"$\partial C_M^2/\partial \beta$", color='red')

        plt.title("(window must be closed to continue Python script)\n")
        ax1.set_xlabel(r"$\alpha$ and $\beta$ $(^o)$")
        ax1.set_ylabel(r"$C_M^2$", color='blue')
        ax2.set_ylabel(r"$\partial C_M^2$", color='red')

        # add legend
        lns = lns1+lns2+lns3+lns4
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc=0)

        fig.tight_layout()  # otherwise the right y-label is slightly clipped
        plt.show()

    except ImportError:
        print ("Unable to import matplotlib.pyplot module or numpy. Plotting is not available")
        plt = None

#######################################
## Setup OpenMDAO                    ##
#######################################

# Setup the openmdao problem object    
omProblem = om.Problem()

# Add the fluid analysis - all the options are added as arguments
avlSystem = avlAnalysis(capsProblem = capsProblem)

omProblem.model.add_subsystem('avlSystem', avlSystem)

# setup the optimization
omProblem.driver = om.ScipyOptimizeDriver()
#omProblem.driver.options['optimizer'] = 'SLSQP'
omProblem.driver.options['optimizer'] = "L-BFGS-B"
omProblem.driver.options['tol'] = 1.e-7
omProblem.driver.options['disp'] = True

# add design variables to model
omProblem.model.add_design_var('avlSystem.Alpha', lower=-25.0, upper=25.0)
omProblem.model.add_design_var('avlSystem.Beta' , lower=-15.0, upper=15.0)

# minimize Cm**2
omProblem.model.add_objective('avlSystem.Cmtot2')

# Press go
omProblem.setup()
omProblem.run_driver()
omProblem.cleanup()


print("Optimized values:", omProblem.get_val("avlSystem.Alpha"), omProblem.get_val("avlSystem.Beta"))

# Check the result
assert abs(omProblem.get_val("avlSystem.Alpha") - -5.505) < 1e-3