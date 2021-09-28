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

# Change the value 
myValue.value = [0.0, 40000.0, 50000]
print("New Value = ",  myValue.value)

# Or change the value using the setVal() function - equivalent to the directly setting the value  
myValue.setVal([0.0, 45000.0, 55000])
print("New Value = ",  myValue.value)

# Get the value using the getVal() function - equivalent to directly calling the value
print("Current Value = ", myValue.getVal())

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()
