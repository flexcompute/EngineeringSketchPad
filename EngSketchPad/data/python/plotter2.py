###################################################################
#                                                                 #
# plotter2.py -- sample line plotting (animation)                 #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

import math
import time

esp.TimLoad("plotter", esp.GetEsp("pyscript"), "")

# moving sine wave
nstep = 32
tdata = ""
sdata = ""
for i in range(nstep+1):
    t = i/nstep * 2 * math.pi
    s = math.sin(t)

    tdata += str(t) + ";"
    sdata += str(s) + ";"

for i in range(nstep+1):
    t = i/nstep * 2 * math.pi
    s = math.sin(t)

    esp.TimMesg("plotter", "new|Moving dot (t="+("%.2f"%t)+")|t|sin(t)|")
    esp.TimMesg("plotter", "add|"+tdata+"|"+sdata+"|k:|")
    esp.TimMesg("plotter", "add|"+str(t)+"|"+str(s)+"|ro|")
    esp.TimMesg("plotter", "show")
    time.sleep(0.5)

esp.TimHold("pyscript", "plotter")
esp.TimQuit("plotter")
