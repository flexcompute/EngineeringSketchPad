from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

# Create a value 
myValue = myProblem.createValue("Altitude", [0.0, 30000.0, 60000.0], units="ft")
print("Name = ",   myValue.name)
print("Value = ",  myValue.value)
print("Units = ",  myValue.units)

# Convert the units from feet to meters
convertValue = myValue.convertUnits("m")
print("Converted Value = ",  convertValue, "(m)")

# Create a new value 
myValue = myProblem.createValue("Freq", [[0.0, 1, 2], [.25, .5, .75]], units="Hz")
print("Name = ",   myValue.name)
print("Value = ",  myValue.value)
print("Units = ",  myValue.units)

# Convert the units from 1/s to 1/min
convertValue = myValue.convertUnits("1 per minute")
print("Converted Value = ",  convertValue, "(1/min)")

# Create a new value 
myValue = myProblem.createValue("Energy", 1.0, units="Btu")
print("Name = ",   myValue.name)
print("Value = ",  myValue.value)
print("Units = ",  myValue.units)

# Convert the units from Btu to Joules
convertValue = myValue.convertUnits("J")
print("Converted Value = ",  convertValue, "(J)")
               
# Close our problems
print("Closing our problem")
myProblem.closeCAPS()