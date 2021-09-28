# Use: Test autolinking of capsValue objects 

from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm", "basicTest")

myAnalysis = myProblem.loadAIM(aim = "fun3dAIM", 
                               analysisDir = "FUN3DAnalysisTest",
                               capsIntent = "CFD")

# Create a value  
myValue = myProblem.createValue("Alpha", -10.0, units="degree")
print("Alpha Value = ",  myValue.value)

# Create another value  
myValue = myProblem.createValue("Mach", 0.25)
print("Mach Value = ",  myValue.value)

# Autolink all capsValues with all AIMs loaded into the problem
myProblem.autoLinkValue()

# Values are updated in a lazy manner - the value shouldn't reflect the linked value yet 
print("Mach = ", myAnalysis.getAnalysisVal("Mach"))
print("Alpha = ", myAnalysis.getAnalysisVal("Alpha"))

# Create a new value  
myValue = myProblem.createValue("Beta", 1.3, units="degree")
print("Beta Value = ",  myValue.value)

# Autolink just the new value with all AIMs loaded into the problem
myProblem.autoLinkValue(myValue) # or could have used - myProblem.autoLinkValue("Beta")

# Run preAnalysis to force the value update
myAnalysis.preAnalysis()

# Re-print the values - the values should reflect the linked values 
print("Again... Mach = ", myAnalysis.getAnalysisVal("Mach"))
print("Again... Alpha = ", myAnalysis.getAnalysisVal("Alpha"))
print("Again... Beta = ", myAnalysis.getAnalysisVal("Beta"))

# Close our problems
print("Closing our problem")
myProblem.closeCAPS()