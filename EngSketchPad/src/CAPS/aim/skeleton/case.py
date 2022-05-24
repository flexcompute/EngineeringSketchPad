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
value = skel.input.SkeletonAIMin
print("Default SkeletonAIMin =", value)
skel.input.SkeletonAIMin = 6.0*m
value = skel.input.SkeletonAIMin
print("Current SkeletonAIMin =", value)

# AIM autoExecutes

# Get an output
value = skel.output.SkeletonAIMout
print("Computed SkeletonAIMout =", value)
