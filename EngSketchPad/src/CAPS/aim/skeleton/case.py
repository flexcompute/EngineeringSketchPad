import pyCAPS

# Instantiate our CAPS problem "myProblem" 
print("Loading file into our Problem")
myProblem = pyCAPS.Problem(problemName="skeletonExample",
                           capsFile="case.csm")

# Load our skeletal aim 
skel = myProblem.analysis.create(aim = "skeletonAIM")

##[analysisValueSetandGet]
# Get current value of our first input
value = skel.input.skeletonAIMin
print("Default skeletonAIMin =", value)
skel.input.skeletonAIMin = 6.0
value = skel.input.skeletonAIMin
print("Current skeletonAIMin =", value)

##[analysisPreAndPost]
# Execute pre-Analysis
skel.preAnalysis()

# Run AIM - os.system call, python interface, etc.....

# Execute post-Analysis (notify CAPS that analysis is done)
skel.postAnalysis()

# Get an output
value = skel.output.skeletonAIMout
print("Computed skeletonAIMout =", value)
##[analysisPreAndPost]
