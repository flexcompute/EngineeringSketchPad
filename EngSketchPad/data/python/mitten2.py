###################################################################
#                                                                 #
# mitten2.py --- start with: serveESP ../data/python/mitten2.py   #
#                or:         python   ../data/python/mitten2.py   #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

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
print("\nmaking modl(../data/basic/mitten1.csm)")
modl = ocsm.Ocsm("../data/basic/mitten1.csm")
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
print("    builtTo:", builtTo);
print("    nbody  :", nbody  );
print("    bodys  :", bodys  );

print("\ncalling modl.Info()")
(nbrch, npmtr, nbody) = modl.Info()
print("    nbrch  :", nbrch);
print("    npmtr  :", npmtr);
print("    nbody  :", nbody);

# load the "mitten" TIM
esp.TimLoad("mitten", esp.GetEsp("pyscript"), "A")

# move block to left
esp.TimMesg("mitten", "xcent|-1|")
esp.TimMesg("mitten", "xcent|-1|")
esp.TimMesg("mitten", "xcent|-1|")
esp.TimMesg("mitten", "xcent|-1|")

# wait 5 seconds
esp.TimMesg("mitten", "countdown|5|")

# save the "mitten" results in the mitten1.csm file
esp.TimSave("mitten")

print("Please choose \"File->Open mitten1.csm\"")
print("Then   choose \"File->Edit mitten1.csm\" to make your changes permanent")

# return modl to ESP
esp.SetModl(modl, esp.GetEsp("pyscript"))
