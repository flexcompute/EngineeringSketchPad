# Use: Load geometry into the problem - compare geometry with returned geometry object

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. The 
# project name "basicTest" may be optionally set here;  if no argument is provided
# the CAPS file provided is used as the project name. 
print("Loading file into our capsProblem")

## [loadCAPS]
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")
## [loadCAPS]

print(myGeometry) 
print(myProblem.geometry)

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
