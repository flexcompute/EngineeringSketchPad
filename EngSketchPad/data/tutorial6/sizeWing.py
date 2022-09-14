###################################################################
#                                                                 #
# sizeWing.py - pick (S, AR, V) to maximize L/D                   #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

# other imports
import math

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model1.csm",
                           outLevel    = 1)

# constants
rho      = 0.002377       # slug/ft3
Wball    = 0.10           # lb/ball
Wfixed   = 0.02           # lb
Wwing    = 0.03           # lb/ft3
CD0      = 0.04
oswald   = 0.90
CL_max   = 1.00

# get the number of balls from the CAPS Problem
if ("nball" in myProblem.parameter):
    nball = myProblem.parameter["nball"].value
else:
    nball = 2

print("nball", nball)

# initialize _best
S_best   = 0
AR_best  = 0
V_best   = 0
W_best   = 0
CL_best  = 0
LoD_best = 0

# loops
for S in [0.25, 0.50, 0.75, 1.00, 1.25, 1.50, 1.75, 2.00, 2.25, 2.50, 2.75, 3.00]:
    for AR in [3, 4, 5, 6, 7, 8]:
        for V in [2.5, 5.0, 7.5, 10.0, 12.5, 15.0, 17.5, 20.0]:

            # compute performance
            b   = math.sqrt(S * AR)
            q   = 0.5 * rho * V**2
            W   = Wfixed + nball * Wball + S * b * Wwing
            CL  = W / (q * S)
            CD  = CD0 + CL**2 / (math.pi * AR * oswald)
            LoD = CL / CD

            # save if best result so far
            if (LoD > LoD_best and CL <= CL_max):
                S_best   = S
                AR_best  = AR
                W_best   = W
                V_best   = V
                CL_best  = CL
                LoD_best = LoD

# report the _best results
print("S_best  =", S_best  )
print("AR_best =", AR_best )
print("W_best  =", W_best  )
print("V_best  =", V_best  )
print("CL_best =", CL_best )
print("LoD_best=", LoD_best)

# put the _best results into the CAPS Problem
myProblem.geometry.despmtr["wing:area"  ].value = S_best
myProblem.geometry.despmtr["wing:aspect"].value = AR_best

if "CLcruise" in myProblem.parameter:
    myProblem.parameter["CLcruise"].value = CL_best
else:
    myProblem.parameter.create("CLcruise", CL_best)
