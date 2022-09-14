###################################################################
#                                                                 #
# test2.py --- start with: serveESP ../data/python/test           #
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
print("    myVol  :", myVol);     assert (abs(myVol[0]-24) < 1e-12)

# change L and rebuild
ipmtr = espModl.FindPmtr("L", 0, 0, 0)
value = espModl.SetValuD(ipmtr, 1, 1, 2.0)

(builtTo, nbody, bodys) = espModl.Build(0, 0)
print("Rebuild with L=2")
print("    builtTo:", builtTo);   assert (builtTo == 4   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

# make sure we get the correct volume now
ipmtr = espModl.FindPmtr("myVol", 0, 0, 0)
myVol = espModl.GetValu(ipmtr, 1, 1)
print("    myVol  :", myVol);     assert (abs(myVol[0]-12) < 1e-12)

# add a sleep so that .py script runs more slowly
print("sleeping for 5 seconds...")
import time
time.sleep(5)
