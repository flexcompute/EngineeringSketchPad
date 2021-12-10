# Import pyCAPS module
import pyCAPS

# Import OpenMDAO module - currently tested against version 2.8.0
from openmdao.api import Problem, Group, IndepVarComp, ExplicitComponent, ScipyOptimizeDriver

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Avl and OpenMDAO v2 PyTest Example',
                                 prog = 'avl_OpenMDAO_2_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "avlOpenMDAO_1")

# Load a *.csm file "../csmData/cfdMultiBody.csm" into our newly created problem.
print("Loading file into our Problem")
geometryScript = os.path.join("..","csmData","avlWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AVL aim
myAnalysis = myProblem.analysis.create(aim = "avlAIM",
                                       name = "avl",
                                       autoExec = False)

# Setup AVL surfaces
wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute
        "numChord"          : 8,
        "spaceChord"        : 1.0,
        "numSpanPerSection" : 12,
        "spaceSpan"         : 1.0}

myAnalysis.input.AVL_Surface = {"Wing": wing}

# Create OpenMDAOComponent - ExternalCode
avlComponent = myAnalysis.createOpenMDAOComponent(["Alpha", "Mach","AVL_Surface:Wing:numChord", # Analysis inputs parameters
                                                   "thick", "area"], # Geometry design variables
                                                  ["CDtot", "CLtot", "Cmtot"], # Output parameters
                                                  changeDir = True, # Change in the analysis directory during execution
                                                  executeCommand = ["avl", "caps"],
                                                  stdin     = "avlInput.txt", # Modify stdin and stdout
                                                  stdout    = "avlOutput.txt", # for avl execution
                                                  setSensitivity = {"type": "fd",  # Set sensitivity information for the
                                                                    "form" : "central", # component
                                                                    "step_size" : 1.0E-3})



# Setup and run OpenMDAO model
print("Setting up and running OpenMDAO model")

# Create a component that calculates the lift to drag ratio
class l2dComponent(ExplicitComponent):
    def __init__(self):
        super(l2dComponent,self).__init__()
        self.add_input("CL", val=0.0)
        self.add_input("CD", val= 0.0)

        self.add_output("L2D", shape=1)
        self.declare_partials(of="*", wrt="*", method="fd", form="central", step=1.E-3)

    def compute(self, params, unknowns):

        unknowns["L2D"] = -params["CL"]/params["CD"] # Minus since we want to maximize L/D in the objective

    def compute_partials(self, params, unknowns):
        J = {}
        J["L2D", "CL"] = -1/params["CD"]
        J["L2D", "CD"] = +params["CL"]/params["CD"]**2

        return J

# Create a Group class for our design problem
class MaxLtoD(Group):
    def __init__(self):
        super(MaxLtoD, self).__init__()

        # Add design variables
        self.add_subsystem('dvAlpha', IndepVarComp('Alpha', float(3.0)))
        self.add_subsystem('dvMach', IndepVarComp('Mach', float(0.1)))
        self.add_subsystem('dvAirfoil_Thickness', IndepVarComp('Thick', float(0.13)))
        self.add_subsystem('dvAirfoil_Area', IndepVarComp('Area', float(10)))

        # Add AVL component
        self.add_subsystem("AVL", avlComponent)

        self.add_subsystem("L2D", l2dComponent())

        # Make connections between design variables and inputs to AVL and geometry
        self.connect("dvAlpha.Alpha", "AVL.Alpha")
        self.connect("dvMach.Mach" , "AVL.Mach")

        self.connect("dvAirfoil_Thickness.Thick", "AVL.thick")
        self.connect("dvAirfoil_Area.Area" , "AVL.area")

        # Connect AVL outputs to L2D calculation inputs
        self.connect("AVL.CDtot", "L2D.CD")
        self.connect("AVL.CLtot", "L2D.CL")

# Iniatate the the OpenMDAO problem
openMDAOProblem = Problem()
openMDAOProblem.model = MaxLtoD()

# Set design driver parameters
openMDAOProblem.driver = ScipyOptimizeDriver()
openMDAOProblem.driver.options['optimizer'] = 'SLSQP'
openMDAOProblem.driver.options['tol'] = 1.0e-8
openMDAOProblem.driver.options['maxiter'] = 500 # Maximum number of solver iterations

# Set bounds on design variables
openMDAOProblem.model.add_design_var('dvAlpha.Alpha', lower=-5.0, upper=10) # Analysis design values
openMDAOProblem.model.add_design_var('dvMach.Mach'  , lower=0.05, upper=0.7)

openMDAOProblem.model.add_design_var('dvAirfoil_Thickness.Thick', lower=0.1, upper=0.5) # Geometry design values
openMDAOProblem.model.add_design_var('dvAirfoil_Area.Area', lower=5.0, upper=20.0)

# Set objective
openMDAOProblem.model.add_objective('L2D.L2D')

# Set Cm constaint
openMDAOProblem.model.add_constraint('AVL.Cmtot', equals= -0.2)

# Setup and run
openMDAOProblem.setup()
openMDAOProblem.run_driver()

# Print the output
print("Design Alpha = ", openMDAOProblem["dvAlpha.Alpha"])
print("Design Mach = " , openMDAOProblem["dvMach.Mach"])
print("Design Thickness = ", openMDAOProblem["dvAirfoil_Thickness.Thick"])
print("Design Area = " , openMDAOProblem["dvAirfoil_Area.Area"])

print("Lift to Drag Ratio = ", -openMDAOProblem["L2D.L2D"])

print("Lift Coeff  = " , myAnalysis.output.CLtot)
print("Drag Coeff  = " , myAnalysis.output.CDtot)
print("Moment (y) Coeff  = ", myAnalysis.output.Cmtot)
