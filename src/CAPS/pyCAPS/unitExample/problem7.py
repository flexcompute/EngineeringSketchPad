# Use: Check setting and getting an Attribute on a problem object  

from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", verbosity = "debug")

# Add attribute
myProblem.addAttribute("testAttr", [1, 2, 3])

# Add another attribute
myProblem.addAttribute("testAttr_2", "anotherAttribute")

myValue = myProblem.getAttribute("testAttr")
print("Value = ",  myValue)

myValue = myProblem.getAttribute("testAttr_2")
print("Value = ",  myValue)

# Close our problems
myProblem.closeCAPS()