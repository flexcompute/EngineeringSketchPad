# Use: Load a analysis module (AIM) into the problem

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. The 
# project name "basicTest" may be optionally set here;  if no argument is provided
# the CAPS file provided is used as the project name. 
print("Loading file into our capsProblem")
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Load FUN3D aim 
## [loadAIM]
myProblem.loadAIM(aim = "fun3dAIM", 
                  altName = "fun3d", 
                  analysisDir = "FUN3DAnalysisTest",
                  capsIntent = "CFD")
## [loadAIM]

## [referenceAIM]
print(myProblem.analysis["fun3d"].analysisDir)
## [referenceAIM]

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
