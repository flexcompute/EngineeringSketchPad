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

inviscid = {"bcType" : "Inviscid", "wallTemperature" : 1.1}
fun3d.setAnalysisVal("Boundary_Condition", [("Skin", inviscid),
                                            ("SymmPlane", "SymmetryY"),
                                            ("Farfield","farfield")])
    
# Create tree 
fun3d.createTree(internetAccess = False, embedJSON = True, analysisGeom = True, reverseMap = True)


# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
