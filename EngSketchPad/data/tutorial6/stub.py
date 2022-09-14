###################################################################
#                                                                 #
# stub.py - stub pyCAPS script to hook up with ESP                #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model1.csm",
                           outLevel    = 1)

################   user script goes here  ########################

#--- myProblem.parameter["badValue"].markForDelete()

myProblem.analysis["skeletonAIM_1"].markForDelete()
myProblem.analysis["skeletonAIM_2"].markForDelete()
myProblem.bound["upperWing"].markForDelete()

##################################################################
