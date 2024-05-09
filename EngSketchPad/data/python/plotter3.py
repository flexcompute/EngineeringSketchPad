###################################################################
#                                                                 #
# plotter3.py -- sample line plotting                             #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

import math

esp.TimLoad("plotter", esp.GetEsp("pyscript"), "")

# data
nstep = 32
xdata = ""
sdata = ""
cdata = ""
for i in range(nstep+1):
    x = i/nstep * 180
    s = math.sin(x * math.pi / 180)
    c = math.cos(x * math.pi / 180)

    xdata += str(x) + ";"
    sdata += str(s) + ";"
    cdata += str(c) + ";"
    
# sine and cosine waves (left axis only)
esp.TimMesg("plotter", "new|Left axis only|theta (deg)|sin and cos|")
esp.TimMesg("plotter", "add|"+xdata+"|"+sdata+"|b-s|")
esp.TimMesg("plotter", "add|"+xdata+"|"+cdata+"|g:+|")
esp.TimMesg("plotter", "show")

esp.TimMesg("plotter", "new|Left (sine) and right (cosine) axes|theta (deg)|sin|cos|")
esp.TimMesg("plotter", "add|"+xdata+"|"+sdata+"|b-s|")
esp.TimMesg("plotter", "add|"+xdata+"|"+cdata+"|g:+2|")
esp.TimMesg("plotter", "show")

esp.TimMesg("plotter", "new|Left axis only|theta (deg)|sin and cos|")
esp.TimMesg("plotter", "add|"+xdata+"|"+sdata+"|b-s|")
esp.TimMesg("plotter", "add|"+xdata+"|"+cdata+"|g:+|")
esp.TimMesg("plotter", "show")

esp.TimQuit("plotter")
