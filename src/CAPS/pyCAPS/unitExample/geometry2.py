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

##[geometryValueGetAll]
# Get all the geometery design parameters 
valueDict = myGeometry.getGeometryVal()
print("Design Parameter Dict =", valueDict)
##[geometryValueGetAll]

##[geometryValueGetNamesAll]
# Get just the names of the geometery design parameters 
valueList = myGeometry.getGeometryVal(namesOnly = True)
print("Design Parameter List =", valueList)
##[geometryValueGetNamesAll]


# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
