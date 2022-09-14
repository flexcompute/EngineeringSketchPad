###################################################################
#                                                                 #
# conserveUpper.py - make/view conservative Upper Bound           #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model5.csm",
                           outLevel    = 1)

# load the aims
if ("skeletonAIM_1" not in myProblem.analysis):
    print("creating skeletonAIM_1")
    myProblem.analysis.create(aim        = "skeletonAIM",
                              capsIntent = "Body_1",
                              name       = "skeletonAIM_1")

if ("skeletonAIM_2" not in myProblem.analysis):
    print("creating skeletonAIM_2")
    myProblem.analysis.create(aim        = "skeletonAIM",
                              capsIntent = "Body_2",
                              name       = "skeletonAIM_2")

# Conservative Bound on upper surface
if ("upperWing" not in myProblem.bound):
    print("creating upperWing")
    boundUpper = myProblem.bound.create("upperWing")
else:
    print("reusing upperWing")
    boundUpper = myProblem.bound["upperWing"]

vset1 = boundUpper.vertexSet.create(myProblem.analysis["skeletonAIM_1"])
vset2 = boundUpper.vertexSet.create(myProblem.analysis["skeletonAIM_2"])

dset1x = vset1.dataSet.create("x",   pyCAPS.fType.FieldOut)
dset2x = vset2.dataSet.create("in1", pyCAPS.fType.FieldIn)
dset2x.link(dset1x, "Conserve")

dset1y = vset1.dataSet.create("y",   pyCAPS.fType.FieldOut)
dset2y = vset2.dataSet.create("in2", pyCAPS.fType.FieldIn)
dset2y.link(dset1y, "Conserve")

dset1z = vset1.dataSet.create("z",   pyCAPS.fType.FieldOut)
dset2z = vset2.dataSet.create("in3", pyCAPS.fType.FieldIn)
dset2z.link(dset1z, "Conserve")

dset1pi = vset1.dataSet.create("pi",  pyCAPS.fType.FieldOut)
dset2pi = vset2.dataSet.create("in4", pyCAPS.fType.FieldIn)
dset2pi.link(dset1pi, "Conserve")

boundUpper.close()

# Get an output (needed to make sure AIM runs)
myProblem.analysis["skeletonAIM_1"].input.num = 16
value = myProblem.analysis["skeletonAIM_1"].output.sqrtNum
print("Computed sqrt =", value)
#myProblem.analysis["skeletonAIM_1"].execute()

# show the flowchart of the AIMS and Bounds
esp.TimLoad("flowchart", esp.GetEsp("pyscript"), "")
esp.TimMesg("flowchart", "show");
esp.TimQuit("flowchart");

# load the viewer
esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")

# view the upperWing Bound
print("==> Viewing x in skeletonAIM_1")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_1|x")

print("==> Viewing in1 in skeletonAIM_2 (conserve)")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_2|in1")

print("==> Viewing y in skeletonAIM_1")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_1|y")

print("==> Viewing in2 in skeletonAIM_2 (conserve)")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_2|in2")

print("==> Viewing z in skeletonAIM_1")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_1|z")

print("==> Viewing in3 in skeletonAIM_2 (conserve)")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_2|in3")

print("==> Viewing pi in skeletonAIM_1")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_1|pi")

print("==> Viewing in4 in skeletonAIM_2 (conserve)")
esp.TimMesg("viewer", "BOUND|upperWing|skeletonAIM_2|in4")

esp.TimQuit("viewer")
