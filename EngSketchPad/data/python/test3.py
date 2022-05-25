###################################################################
#                                                                 #
# test3.py --- start with: serveESP ../data/python/test           #
#                  Tool->Pyscript   ../data/python/test1.py       #
#                  Tool->Pyscript   ../data/python/test2.py       #
#                  Tool->Pyscript   ../data/python/test3.py       #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyOCSM import ocsm
from pyOCSM import esp

# espModl is inherited from serveESP
espModl = ocsm.Ocsm(esp.GetModl(esp.GetEsp("pyscript")))

# get the current volume
ipmtr = espModl.FindPmtr("myVol", 0, 0, 0)
myVol = espModl.GetValu(ipmtr, 1, 1)
print("From previous build")
print("    myVol:", myVol);   assert (abs(myVol[0]-12) < 1e-12)

# set front face to magenta
espModl.SetAttr(4, "_color", "$magenta")
ocsm.SetOutLevel(0)
espModl.Build(0, 0)
ocsm.SetOutLevel(1)

# write newTest.csm
print("Creating \"newTest.csm\"")
espModl.Save("newTest.csm")

# generate a purposeful error
purposefulError()
