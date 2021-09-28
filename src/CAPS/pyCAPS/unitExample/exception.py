# Use: Test autolinking of capsValue objects 

from __future__ import print_function

# Import pyCAPS module (Linux and OSx = pyCAPS.so file; Windows = pyCAPS.pyd file) 
import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Initiating capsProblem")
myProblem = pyCAPS.capsProblem()

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. 
print("Loading file into our capsProblem")
myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm")

# Try to load another geometry into the problem - this is forbidden 
try: 
    myGeometry = myProblem.loadCAPS("./csmData/cfdMultiBody.csm")

except pyCAPS.CAPSError as e:
    print("We caught the error!")
    
    print("Error code =", e.errorCode)
    print("Error name =", e.errorName)
    
    # Do something based on error captured.
    
except:
    print ("We did not catch the error!")
