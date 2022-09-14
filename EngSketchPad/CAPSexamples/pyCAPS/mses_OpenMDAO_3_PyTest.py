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

class msesAnalysis(om.ExplicitComponent):
    
    def initialize(self):
    # Setup variables so that arguments can be provided to the component.
        
        # CAPS Problem input
        self.options.declare('capsProblem', types=object)

    def setup(self):
    # Assign the actual values to the variables.
 
        capsProblem = self.options['capsProblem']

        # attach parameters to OpenMDAO object
        self.add_input('camber', val=capsProblem.geometry.despmtr.camber)
        
        # Add output metric
        self.add_output('CD')

    def compute(self, inputs, outputs):

        capsProblem = self.options['capsProblem']
        mses = capsProblem.analysis['mses']

        # Update input values
        capsProblem.geometry.despmtr.camber = inputs['camber']
    
        # Grab objective and attach as an output
        outputs['CD']  = mses.output.CD
        print("camber", inputs['camber'], "CD", outputs['CD'])

    def setup_partials(self):
    # declare the partial derivatives
        
        # Declare and attach partials to self
        self.declare_partials('CD','camber')

    def compute_partials(self,inputs,partials):
    # The actual value of the partial derivative is assigned here. 
    # Assume input value are consisten with last call to compute
        
        capsProblem = self.options['capsProblem']
        mses = capsProblem.analysis['mses']

        # Get derivatives and set partials
        partials['CD', 'camber'] = mses.output["CD"].deriv("camber")
        print("CD_camber", partials['CD', 'camber'])


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'mses and OpenMDAO v3 PyTest Example',
                                 prog = 'mses_OpenMDAO_3_PyTest',
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
problemName = os.path.join(args.workDir[0], "msesOpenMDAO_v3")
capsProblem = pyCAPS.Problem(problemName="msesOpenMDAO_v3",
                             capsFile=os.path.join("..","csmData","airfoilSection.csm"),
                             outLevel=args.outLevel)

# Set a non-zero camber. The optimization should drive it back to zero
capsProblem.geometry.despmtr.camber = 0.1

# Setup MSES
mses = capsProblem.analysis.create(aim = "msesAIM", name="mses")

# Set flow condition
mses.input.Alpha = 0.0
mses.input.Mach  = 0.5
mses.input.Re = 5e6

# Trip the flow near the leading edge to get smooth gradient
mses.input.xTransition_Upper = 0.1
mses.input.xTransition_Lower = 0.1

# Use camber as the design variable
mses.input.Design_Variable = {"camber":{}}

# Plot the functional and gradient if requested
if args.plotFunctional:
    try:
        import numpy as npy
        import matplotlib.pyplot as plt

        cambers = npy.linspace(-0.1, 0.1, 15)
        CD = []
        CD_camber = []
        for camber in cambers:
            capsProblem.geometry.despmtr.camber = camber
        
            CD.append(       mses.output["CD"].value )
            CD_camber.append(mses.output["CD"].deriv("camber") )

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        
        # Compute central difference derivatives
        dCD_camber = []
        for i in range(1,len(cambers)-1):
            dCD_camber.append( (CD[i+1]-CD[i-1])/(cambers[i+1]-cambers[i-1]) )

        # Plot the functional
        lns1 = ax1.plot(cambers, CD, 'o--', label = r"$C_D(\mathrm{camber})$", color='blue')
        
        # Plot plot the derivative
        lns2 = ax2.plot(cambers, CD_camber, 's-', label = r"$\partial C_D/\partial \mathrm{camber}$", color='red')

        # Plot plot the FD derivative
        lns3 = ax2.plot(cambers[1:-1], dCD_camber, '.--', label = r"$\Delta C_D/\Delta \mathrm{camber}$", color='orange')

        ax2.axhline(0., linewidth=1, linestyle=":", color='r')
        #ax2.axvline(0., linewidth=1, linestyle=":", color='k')

        plt.title("(window must be closed to continue Python script)\n")
        ax1.set_xlabel(r"camber")
        ax1.set_ylabel(r"$C_D$", color='blue')
        ax2.set_ylabel(r"$\partial C_D$", color='red')

        # add legend
        lns = lns1+lns2+lns3
        labs = [l.get_label() for l in lns]
        ax1.legend(lns, labs, loc='best', facecolor='white', framealpha=1)

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
msesSystem = msesAnalysis(capsProblem = capsProblem)

omProblem.model.add_subsystem('msesSystem', msesSystem)

# setup the optimization
omProblem.driver = om.ScipyOptimizeDriver()
#omProblem.driver.options['optimizer'] = 'SLSQP'
omProblem.driver.options['optimizer'] = "L-BFGS-B"
omProblem.driver.options['tol'] = 1.e-7
omProblem.driver.options['disp'] = True

# add design variables to model
omProblem.model.add_design_var('msesSystem.camber', lower=-0.1, upper=0.1)

# minimize CD
omProblem.model.add_objective('msesSystem.CD')

# Press go
omProblem.setup()
omProblem.run_driver()
#omProblem.check_partials(step=0.0001)
omProblem.cleanup()


print("Optimized values:", omProblem.get_val("msesSystem.camber"))

# Check the result
assert abs(omProblem.get_val("msesSystem.camber")) < 1e-3