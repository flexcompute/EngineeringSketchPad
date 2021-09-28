# Import pyCAPS module
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. The 
# project name "basicTest" may be optionally set here;  if no argument is provided
# the CAPS file provided is used as the project name. 
print("Loading file into our capsProblem")
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Load Tetgen aim 
tetgen = myProblem.loadAIM(aim = "tetgenAIM", 
                           altName = "tetgen", 
                           analysisDir = "TetGenAnalysisTest",
                           capsIntent = "CFD")

# Analysis information
state = tetgen.getAnalysisInfo()
print ("Cleanliness state", state) # Cleanliness status

# Get analysis information dictionary 
infoDict = tetgen.getAnalysisInfo(printinfo = False,  infoDict = True)

# Print dictionary 
print(infoDict)