# Import pyCAPS module
import pyCAPS

# Import OpenMDAO v3 module
import openmdao
import openmdao.api as om

# Import other modules
import os
import argparse

if int(openmdao.__version__[:openmdao.__version__.find(".")]) < 3:
    print( "This example is designed for OpenMDAO version 3")
    print( "Installed OpenMDAO :", openmdao.__version__ )
    exit()

class avlAnalysis(om.ExplicitComponent):
    
    def initialize(self):
    # Setup variables so that arguments can be provided to the component.
        
        # CAPS Problem inputs
        self.options.declare('workDir' , types=str)
        self.options.declare('outLevel', types=int)
        
        # Free design parameters
        self.options.declare('Alpha', types=float)
        self.options.declare('Beta' , types=float)

    def setup(self):
    # Assign the actual values to the variables.
 
        # attach parameters to OpenMDAO fluid object
        self.add_input('Alpha', val=self.options['Alpha'])
        self.add_input('Beta' , val=self.options['Beta'])
        
        # Add output metric
        self.add_output('Cmtot2')

        # Initiate problem and object
        problemName = os.path.join(self.options['workDir'], "avlOpenMDAO_v3")
        self.capsProblem = pyCAPS.Problem(problemName="avlOpenMDAO_v3",
                                          capsFile=os.path.join("..","csmData","avlWing.csm"),
                                          outLevel=self.options['outLevel'])

        # Setup AVL
        self.avl = self.capsProblem.analysis.create(aim = "avlAIM")

        # Setup AVL surfaces
        wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
                "numChord"          : 11,
                "spaceChord"        : 1.0,
                "numSpanPerSection" : 12,
                "spaceSpan"         : 1.0}
        
        self.avl.input.AVL_Surface = {"Wing": wing}
        
        self.avl.input.Alpha = self.options['Alpha']
        self.avl.input.Beta  = self.options['Beta']

    def compute(self, inputs, outputs):
                
        # Update input values
        self.avl.input.Alpha = inputs['Alpha']
        self.avl.input.Beta  = inputs['Beta']
        
        # Grab objective and attach as an output
        outputs['Cmtot2']  = self.avl.output.Cmtot**2
        #print("Alpha", inputs['Alpha'], "Beta", inputs['Beta'], "Cmtot2", outputs['Cmtot2'])

    def setup_partials(self):
    # declare the partial derivatives
        
        # Declare and attach partials to self
        self.declare_partials('Cmtot2','Alpha')
        self.declare_partials('Cmtot2','Beta')

    def compute_partials(self,inputs,partials):
    # The actual value of the partial derivative is assigned here. 
        # Update input values
        self.avl.input.Alpha = inputs['Alpha']
        self.avl.input.Beta  = inputs['Beta']

        # Get derivatives and set partials
        Cm  = self.avl.output["Cmtot"].value
        Cma = self.avl.output["Cmtot"].deriv("Alpha")
        Cmb = self.avl.output["Cmtot"].deriv("Beta")

        partials['Cmtot2', 'Alpha'] = 2*Cm*Cma
        partials['Cmtot2', 'Beta']  = 2*Cm*Cmb
        #print("Cmtot2_Alpha", partials['Cmtot2', 'Alpha'])
        #print("Cmtot2_Beta", partials['Cmtot2', 'Beta'])


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Avl and OpenMDAO v3 PyTest Example',
                                 prog = 'avl_OpenMDAO_3_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()
 
# Setup the openmdao problem object    
omProblem = om.Problem()

# Add the fluid analysis - all the options are added as arguments
avlSystem = avlAnalysis(workDir  = str(args.workDir[0]),
                        outLevel = args.outLevel,
                        Alpha    = 0.,
                        Beta     = 3.)

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
