###################################################################
#                                                                 #
# optTaper.py - pick wing:taper to maximize oswald                #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# get access to pyCAPS
import pyCAPS
from   pyOCSM import esp

import time

# if running in serveESP, the following load is ignored
myProblem = pyCAPS.Problem(problemName = "RunWithoutESP",
                           capsFile    = "../data/tutorial6/model2.csm",
                           outLevel    = 1)

# create the AVL aim
if ("avl" in myProblem.analysis):
    avl = myProblem.analysis["avl"]
else:
    avl = myProblem.analysis.create(aim  = "avlAIM",
                                    name = "avl")

# set the non-geometric AVL inputs
if ("CLcriuse" in myProblem.parameter):
    avl.input["CL"].value = myProblem.parameter["CLcruise"].value
else:
    avl.input["CL"].value = 0.95
    
avl.input["Mach"].value = 0
avl.input["Beta"].value = 0
avl.input.AVL_Surface   = {"wing" : {"numChord"     : 4,
                                     "numSpanTotal" : 24}}

# initialize _best
taper_best  = 0
oswald_best = 0

# initialize _data
taper_data  = ""
oswald_data = ""

# loop over the tapers
for taper in [0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    print("executing taper =", taper)

    # compute performance (getting output automagically runs AVL)
    myProblem.geometry.despmtr["wing:taper"].value = taper
    oswald = avl.output["e"].value

    # keep track of plotter data
    taper_data  += str(taper)  + ";"
    oswald_data += str(oswald) + ";"

    # save if best result so far
    if (oswald > oswald_best):
        taper_best  = taper
        oswald_best = oswald

    # sleep for 2 seconds so that we can suspend
    time.sleep(2)

# plot the _data
esp.TimLoad("plotter", esp.GetEsp("pyscript"), "")
esp.TimMesg("plotter", "new|Taper optimization|taper ratio|oswald|")
esp.TimMesg("plotter", "add|"+taper_data+"|"+oswald_data+"|k:+|")
esp.TimMesg("plotter", "add|"+str(taper_best)+"|"+str(oswald_best)+"|ro|")
esp.TimMesg("plotter", "show")
esp.TimQuit("plotter")
            
# report the results
print("taper_best =", taper_best )
print("oswald_best=", oswald_best)

# put the _best result into the CAPS Problem
print("re-setting optimal taper")
myProblem.geometry.despmtr["wing:taper"].value = taper_best
myProblem.geometry.build()

if ("oswald" in myProblem.parameter):
    myProblem.parameter["oswald"].value = oswald_best
else:
    myProblem.parameter.create("oswald", oswald_best)
