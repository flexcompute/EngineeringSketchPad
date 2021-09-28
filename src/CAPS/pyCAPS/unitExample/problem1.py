# Use: Create a multiple problems with geometry

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblems")
myProblem = pyCAPS.capsProblem()

# Create another problem 
myProblemNew = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. The 
# project name "basicTest" may be optionally set here;  if no argument is provided
# the CAPS file provided is used as the project name. 
print("Loading file into our capsProblems")

## [loadCAPS]
myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")
## [loadCAPS]

myProblemNew.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Close our problems
print("Closing our problems")
myProblem.closeCAPS()
myProblemNew.closeCAPS()
