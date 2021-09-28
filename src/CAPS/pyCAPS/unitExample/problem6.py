# Use: Check creating a tree on the problem
import pyCAPS

# Load a *.csm file "./csmData/cfdMultiBody.csm" into our newly created problem.
myProblem = myProblem.Porblem(problemName = "outLevelExample", 
                              capsFile="csmData/cfdMultiBody.csm", 
                              outLevel="debug")

# Create problem tree 
myProblem.createTree()
