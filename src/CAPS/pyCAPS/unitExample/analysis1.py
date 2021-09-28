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
                          analysisDir = "FUN3DAnalysisTest",
                          capsIntent = "CFD")

# Print out memory
print(fun3d)

print("Duplicate")
# Try to load the same AIM twice
fun3d = myProblem.loadAIM(aim = "fun3dAIM",
                          analysisDir = "FUN3DAnalysisTest2",
                          capsIntent = "CFD")
# Print out memory 
print(fun3d)

print(myProblem.analysis["fun3dAIM"])
print(myProblem.analysis.keys())

print("Overwrite")
# Overwrite aim to load the same AIM twice
fun3d = myProblem.loadAIM(aim = "tetgenAIM", 
                          altName = "tetgen", 
                          analysisDir = "FUN3DAnalysisTest",
                          capsIntent = "CFD")

print(fun3d)

print(myProblem.analysis["fun3dAIM"])
print(myProblem.analysis["fun3dAIM_1"])

print(myProblem.analysis["tetgen"])

print("Redundant directories")

# Try to load the same AIM twice with redundant directories expect an error. 
fun3d = myProblem.loadAIM(aim = "fun3dAIM",
                          analysisDir = "FUN3DAnalysisTest",
                          capsIntent = "CFD")
