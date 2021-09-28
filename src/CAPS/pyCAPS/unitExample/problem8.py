# Use: Load a analysis module (AIM) into the problem, save the problem, then use that CAPS file
# to initiate a new problem

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
myAnalysis = myProblem.loadAIM(aim = "fun3dAIM", 
                               altName = "fun3d", 
                               analysisDir = "FUN3DAnalysisTest",
                               capsIntent = "CFD")

# Run pre and post to get the analysis in a "clean" state
myAnalysis.preAnalysis()
myAnalysis.postAnalysis() 

# Add a value/parameter to problem 
myValye = myProblem.createValue("Mach", 0.5)

# Save our problem - use default name
print("Saving out problem")
myProblem.saveCAPS()

# Close our problem
print("Closing our problem")
myProblem.closeCAPS()

#Create a new problem
myProblem = pyCAPS.capsProblem()

# Reload the saved problem
print("Reload the saved problem")
myGeometry = myProblem.loadCAPS("saveCAPS.caps")

# Get a copy of the analysis object 
myAnalysis = myProblem.analysis["fun3d"]

# Get analysis 
myAnalysis.getAnalysisInfo()

# Check other features of analysis object
print(myAnalysis.aimName)
print(myAnalysis.analysisDir)

# Get a copy of the value object 
myValue = myProblem.value["Mach"]
print("Value =", myValue.value)
