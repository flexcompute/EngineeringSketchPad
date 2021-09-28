# Use: Retrieve attributes from the geometry loaded into an analysis object 

# Make print statement Python 2-3 agnostic
from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm")

# Load FUN3D aim 
fun3d = myProblem.loadAIM(aim = "fun3dAIM", 
                          altName = "fun3d", 
                          analysisDir = "FUN3DAnalysisTest",
                          capsIntent = "CFD")

# Retrieve "capsGroup" attributes for geometry loaded into fun3dAIM. 
# Only search down to the body level 
attributeList = fun3d.getAttributeVal("capsGroup", attrLevel = "Body") 
print(attributeList)

# Provide a wrong attribute level - issue warning and default to "Face"
attributeList = fun3d.getAttributeVal("capsGroup", attrLevel = 5) 
print(attributeList)

# Retrieve "capsGroup" attributes for geometry loaded into fun3dAIM. 
# Only search down to the node level 
attributeList = fun3d.getAttributeVal("capsGroup", attrLevel = "Node") 
print(attributeList)
