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

##[analysisValueGetAll]
# Get all input values
valueDict = fun3d.getAnalysisVal()

print("Variable Dict = ", valueDict)
##[analysisValueGetAll]

##[analysisValueGetNamesAll]
# Get just the names of the input values
valueList = fun3d.getAnalysisVal(namesOnly = True)

print("Variable List =", valueList)
##[analysisValueGetNamesAll]


##[analysisOutValueGetNamesAll]
# Get just the names of the output values
valueList = fun3d.getAnalysisOutVal(namesOnly = True)

print("Variable List =", valueList)
##[analysisOutValueGetNamesAll]

##[analysisOutValueGetAll]
# Get all output values - pre- and post-analysis must have been run
valueDict = fun3d.getAnalysisOutVal()

print("Variable Dict = ", valueDict)
##[analysisOutValueGetAll]
