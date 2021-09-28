# Use: Duplicate an AIM 

from __future__ import print_function

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
fun3d = myProblem.loadAIM(aim = "fun3dAIM", 
                          altName = "fun3d", 
                          analysisDir = "FUN3DAnalysisTest",
                          capsIntent = "CFD")

# Print out memory
print(fun3d)

print("Duplicate")
# Duplicate the fun3d aim using the caps_dupAnalysis function
fun3d2 = myProblem.loadAIM(altName = "fun3d2", 
                           analysisDir = "FUN3DAnalysisTest2",
                           copyAIM = "fun3d")
# Print out memory 
print(fun3d2)
