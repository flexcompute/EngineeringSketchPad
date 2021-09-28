from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Import OpenMDAO module
from openmdao.api import Problem, Group, IndepVarComp

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. 
print("Loading file into our capsProblem")
myGeometry = myProblem.loadCAPS("./csmData/avlWing.csm")

# Load AVL aim 
myAnalysis = myProblem.loadAIM(aim = "avlAIM", 
                               altName = "avl", 
                               analysisDir = "OpenMDAOUnitTest",
                               capsIntent = "ALL")

# Create OpenMDAOComponet - external code 
myAnalysis.createOpenMDAOComponent(["Alpha", "Mach", "thick"], 
                                   ["CDtot"], 
                                   executeCommand = ["avl", myAnalysis.analysisDir + "/caps"],
                                   inputFile = myAnalysis.analysisDir + "/avlInput.txt", 
                                   stdin     = myAnalysis.analysisDir + "/avlInput.txt", # Modify stdin and stdout
                                   stdout    = myAnalysis.analysisDir + "/avlOutput.txt") # for avl execution 

# Setup avl surfaces 
wing = {"groupName"         : "Wing", # Notice Wing is the value for the capsGroup attribute 
        "numChord"          : 8, 
        "spaceChord"        : 1.0, 
        "numSpanPerSection" : 12, 
        "spaceSpan"         : 1.0}

myAnalysis.setAnalysisVal("AVL_Surface", [("Wing", wing)])

# Setup and run OpenMDAO model
print("Setting up and running OpenMDAO model")
openMDAOProblem = Problem()
openMDAOProblem.root = Group()

# Create and connect inputs
openMDAOProblem.root.add('param1', IndepVarComp('Alpha', 3.0))
openMDAOProblem.root.add('param2', IndepVarComp('Mach',  0.5))
openMDAOProblem.root.add('param3', IndepVarComp('thick', 0.15))
openMDAOProblem.root.add('myAnalysis', myAnalysis.openMDAOComponent)

openMDAOProblem.root.connect('param1.Alpha', 'myAnalysis.Alpha')
openMDAOProblem.root.connect('param2.Mach',  'myAnalysis.Mach')
openMDAOProblem.root.connect('param3.thick', 'myAnalysis.thick')

# Run the ExternalCode Component
openMDAOProblem.setup()
openMDAOProblem.run()

# Print the output
print("Unknown in OpenMDAO = ", openMDAOProblem.root.myAnalysis.unknowns['CDtot'])

# Close our CAPS problem
print("Closing our CAPS problem")
myProblem.closeCAPS()