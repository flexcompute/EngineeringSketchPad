###################################################################
#                                                                 #
# bendWing.py - find wing tip deflection for given weight         #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

import math

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model3.csm",
                           outLevel    = 1)

# constants
Eyoung   = 77.5e6         # lb/ft2

# design variables
sparrad  = 0.0104         # ft
W        = 0.36           # lb
b        = myProblem.geometry.outpmtr["wing:span"].value

# compute performance
Inertia  = math.pi / 4 * sparrad**4
deflect  = (W/2) * (b/2)**3 / (3 * Eyoung * Inertia)

# report the results
print("deflect=", deflect)

# put the result into the CAPS Problem
myProblem.geometry.despmtr["wing:ztip"].value = deflect
myProblem.geometry.build()
