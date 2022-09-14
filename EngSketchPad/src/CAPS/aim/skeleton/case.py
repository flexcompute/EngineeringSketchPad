import pyCAPS

# Define length unit
m = pyCAPS.Unit("m")

# Instantiate our CAPS problem "myProblem"
print("Loading file into our Problem")
myProblem = pyCAPS.Problem(problemName="skeletonExample",
                           capsFile="case.csm")

# Load our skeletal aim
skel = myProblem.analysis.create(aim = "skeletonAIM")

# Get current value of our first input
value = skel.input.num
print("Default num =", value)
skel.input.num = 16.0
value = skel.input.num
print("Current num =", value)

# AIM autoExecutes

# Get an output
value = skel.output.sqrtNum
print("Computed sqrtNum =", value)
