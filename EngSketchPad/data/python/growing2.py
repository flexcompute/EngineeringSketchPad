###################################################################
#                                                                 #
# growing2.py --- start with: serveESP ../data/python/growing2.py #
#                     Tool->Pyscript   ../data/python/growing2.py #
#                 or:         python   ../data/python/growing2.py #
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
print("\nmaking modl(../data/python/growing.csm)")
modl = ocsm.Ocsm("../data/python/growing.csm")

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
print("    builtTo:", builtTo);
print("    nbody  :", nbody  );
print("    bodys  :", bodys  );

print("\ncalling modl.Info()")
(nbrch, npmtr, nbody) = modl.Info()
print("    nbrch  :", nbrch);
print("    npmtr  :", npmtr);
print("    nbody  :", nbody);

# inform ESP of modl
esp.SetModl(modl, esp.GetEsp("pyscript"))

# loop through successively larger cylinders
esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")
for xsize in [0.5, 1.0, 1.5, 2.001, 2.5, 3.0]:
    esp.TimMesg("viewer", "MODL")
    esp.TimHold("pyscript", "viewer")

    ipmtr = modl.FindPmtr("xsize", 0, 0, 0)
    modl.SetValuD(ipmtr, 1, 1, xsize)

    modl.Build(0, 0)

# show final configuration
esp.TimMesg("viewer", "MODL")
esp.TimQuit("viewer")

