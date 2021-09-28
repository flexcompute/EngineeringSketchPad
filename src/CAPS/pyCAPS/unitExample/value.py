from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Create a value 
myValue = myProblem.createValue("Altitude", [0.0, 30000.0, 60000.0], units="ft", limits=[0.0, 70000])

print("Name = ",   myValue.name)
print("Value = ",  myValue.value)
print("Units = ",  myValue.units)
print("Limits = ", myValue.limits)
