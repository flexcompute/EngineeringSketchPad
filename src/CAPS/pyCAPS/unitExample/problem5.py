#Use case set verbosity of the problem
import pyCAPS

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem. The 
# project name "basicTest" may be optionally set here;  if no argument is provided
# the CAPS file provided is used as the project name. 
print("Loading file into our Problem")
myProblem = myProblem.Porblem(problemName = "outLevelExample", 
                              capsFile="csmData/cfdMultiBody.csm", 
                              outLevel="debug")


# Change verbosity to minimal - 0 (integer value) 
myProblem.setOutLevel("minimal")

# Change verbosity to standard - 1 (integer value)
myProblem.setOutLevel("standard")

# Change verbosity to back to minimal using integer value - 0
myProblem.setOutLevel(0)

# Change verbosity to back to debug using integer value - 2
myProblem.setOutLevel(2)

# Give wrong value (raises and Error)
myProblem.setOutLevel(10)