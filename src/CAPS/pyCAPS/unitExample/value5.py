from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Create a value - allow the shape and length to change 
myValue = myProblem.createValue("Mach", 0.0, fixedLength=False, fixedShape=False)

print("Value = ",  myValue.value)

myValue.value = [0.0, 1.0] # Try to change the length - should be ok  
print("Value = ",  myValue.value)

# Create a value - fix the shape but not the length 
myValue = myProblem.createValue("Altitude", [0.0, 30000.0, 60000.0], fixedLength=False)

print("Value = ",  myValue.value)

myValue.value = [100000.0, 200000.0] # Try to change the length - should be ok
print("Value = ",  myValue.value)

# Create an array of strings value 
myValue = myProblem.createValue("Strings", ["test1", "test2", "test3"])
print("Value = ",  myValue.value)

# Create a value - fix the shape and the length 
myValue = myProblem.createValue("Alpha", [0.0, 3.0, 6.0])

print("Value = ",  myValue.value)

myValue.value = 0.0 # Try to change the length/shape - should throw an error 
print("Value = ",  myValue.value)

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()