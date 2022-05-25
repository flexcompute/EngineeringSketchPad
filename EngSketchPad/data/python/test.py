###################################################################
#                                                                 #
# test.py --- start with: serveESP ../data/python/test.py         #
#                 Tool->Pyscript   ../data/python/test.py         #
#              or:        python   ../data/python/test.py         #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

# set tolerance for assertions
TOL = 1e-3

# callback functions
def pyMesgCB(text):
    print(" ")
    print("======= in pyMesgCB =======")
    print("   ", text.decode())
    print("===========================")
    return

def pySizeCB(modl, ipmtr, nrow, ncol):
    print(" ")
    print("======= in pySizeCB =======")
    print("    ipmtr:", ipmtr)
    print("    nrow :", nrow )
    print("    ncol :", ncol )
    print("===========================")
    return

print("\ncalling ocsm.Version()")
(imajor, iminor) = ocsm.Version()
print("    imajor :", imajor)
print("    iminor :", iminor)

# make a new Modl
print("\nmaking modl(../data/python/test.csm)")
modl = ocsm.Ocsm("../data/python/test.csm")

esp.SetModl(modl, esp.GetEsp("pyscript"))

print("\ncalling modl.RegMesgCB(pyMesgCB)")
modl.RegMesgCB(pyMesgCB)

print("\ncalling modl.Check()")
modl.Check()

print("\ncalling modl.GetFilelist()")
filelist = modl.GetFilelist()
print("    filelist:", filelist)

print("\ncalling modl.Build(0, 0)")
(builtTo, nbody, bodys) = modl.Build(0, 0)
print("")
print("    builtTo:", builtTo);  assert (builtTo == 4   )
print("    nbody  :", nbody  );  assert (nbody   == 0   )
print("    bodys  :", bodys  );  assert (bodys   == None)

print("\ncalling modl.Info()")
(nbrch, npmtr, nbody) = modl.Info()
print("    nbrch  :", nbrch);    assert (nbrch == 4 )
print("    npmtr  :", npmtr);    assert (npmtr == 48)
print("    nbody  :", nbody);    assert (nbody == 1 )

# view the MODL
esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")
esp.TimMesg("viewer", "MODL")
esp.TimQuit("viewer")

# add a sleep so that .py script runs more slowly
print("sleeping for 5 seconds...")
import time
time.sleep(5)

# make sure we get the correct volume
ipmtr = modl.FindPmtr("myVol", 0, 0, 0)
myVol = modl.GetValu(ipmtr, 1, 1)
print("    myVol  :", myVol);     assert (abs(myVol[0]-24) < 1e-12)

# visualize the model in ESP
#esp.ViewModl(modl)

# change L and rebuild
ipmtr = modl.FindPmtr("L", 0, 0, 0)
value = modl.SetValuD(ipmtr, 1, 1, 2.0)

print("Rebuild with L=2")
(builtTo, nbody, bodys) = modl.Build(0, 0)
print("")
print("    builtTo:", builtTo);   assert (builtTo == 4   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

# make sure we get the correct volume now
ipmtr = modl.FindPmtr("myVol", 0, 0, 0)
myVol = modl.GetValu(ipmtr, 1, 1)
print("    myVol  :", myVol);     assert (abs(myVol[0]-12) < 1e-12)

# visualize the new (smaller) modl
#esp.ViewModl(modl)

# return modl to ESP
#esp.SetModl(modl, esp.GetEsp("pyscript"))
