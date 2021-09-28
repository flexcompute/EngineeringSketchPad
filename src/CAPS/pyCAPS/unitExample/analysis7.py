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

# Load Tetgen aim 
tetgen = myProblem.loadAIM(aim = "tetgenAIM", 
                           altName = "tetgen", 
                           analysisDir = "TetGenAnalysisTest",
                           capsIntent = "CFD")

# Add attribute
tetgen.addAttribute("testAttr", [1, 2, 3])

# Add another attribute
tetgen.addAttribute("testAttr_2", "anotherAttribute")

myValue = tetgen.getAttribute("testAttr")
print("Value = ",  myValue)

myValue = tetgen.getAttribute("testAttr_2")
print("Value = ",  myValue)

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
