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

##[geometryValueSetandGet]
# Get current value of the leading edges sweep 
value = myGeometry.getGeometryVal("lesweep")
print("Current sweep Value =", value)
    
# Set a new value for the leading edges sweep 
myGeometry.setGeometryVal("lesweep", value + 20.0)

# Check to see if the value was set correctly
value = myGeometry.getGeometryVal("lesweep")
print("New Sweep Value =", value)
##[geometryValueSetandGet]

##[geometryView]
myGeometry.viewGeometry(filename = "GeomViewDemo_NewSweep",
                        title="DESPMTR: lesweep = " + str(value), 
                        showImage = True, 
                        ignoreBndBox = True, 
                        combineBodies = True, 
                        showTess = False,
                        showAxes = False,
                        directory = "GeometryImages",
                        viewType = "fourview")
##[geometryView]

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
