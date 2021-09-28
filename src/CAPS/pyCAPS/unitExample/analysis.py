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


##[analysisValueSetandGet]
# Get current value of the project name
value = fun3d.getAnalysisVal("Proj_Name")
print("Current project name =", value)
    
# Set a new value for the project name
fun3d.setAnalysisVal("Proj_Name", "pyCAPS_Demo")

# Check to see if the value was set correctly
value = fun3d.getAnalysisVal("Proj_Name")
print("New project name =", value)
##[analysisValueSetandGet]

# Set overwrite nml to true 
fun3d.setAnalysisVal("Overwrite_NML", True)

##[analysisPreAndPost]
# Execute pre-Analysis
fun3d.preAnalysis()

# Run AIM - os.system call, python interface, etc.....

# Execute post-Analysis   
fun3d.postAnalysis()
##[analysisPreAndPost]
