###################################################################
#                                                                 #
# shapes.py --- start with: serveESP ../data/python/box           #
#                   Tool->Pyscript   ../data/python/shapes.py     #
#                   Tool->Pyscript   ../data/python/shapes.py     #
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

# boxModl is inherited from serveESP
print("\ninheriting boxModl")
boxModl = ocsm.Ocsm(esp.GetModl(esp.GetEsp("pyscript")))
(builtTo, nbody, bodys) = boxModl.Build(0, 0)
print("")
print("    builtTo:", builtTo);
print("    nbody  :", nbody  );
print("    bodys  :", bodys  );

# cylModl is from a file
print("\nloading cylModl")
cylModl = ocsm.Ocsm("../data/python/cyl.csm")
(nbrch, npmtr, nbody) = cylModl.Info()
cylModl.NewBrch(nbrch, cylModl.GetCode("select"), "<none>", 0,
                "face", "",  "",  "",  "",  "",  "",  "",  "")
cylModl.SetAttr(nbrch+1, "_makeQuads", "1")
(builtTo, nbody, bodys) = cylModl.Build(0, 0)
print("")
print("    builtTo:", builtTo);
print("    nbody  :", nbody  );
print("    bodys  :", bodys  );

# sphModl is made from scratch
print("\nmaking sphModl")
sphModl = ocsm.Ocsm("")
i = sphModl.FindPmtr("radius", ocsm.DESPMTR, 1, 1)
sphModl.SetValuD(i, 1, 1, 2)
sphModl.NewBrch(0, sphModl.GetCode("sphere"), "<none>", 0,
                "0", "0", "0", "radius", "", "", "", "", "")
sphModl.SetAttr(1, "_makeQuads", "1")
(builtTo, nbody, bodys) = sphModl.Build(0, 0)
print("")
print("    builtTo:", builtTo);
print("    nbody  :", nbody  );
print("    bodys  :", bodys  );

# conModl comes from pyEgads
print("\nloading conModl from EGADS model")
context = egads.Context()
model   = context.loadModel("../data/python/cone.egads")
conModl = ocsm.Ocsm(model)

# loop through the modls with the viewer
esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")
for i in range(3):
    print("showing box...")
    esp.SetModl(boxModl, esp.GetEsp("pyscript"))
    esp.TimMesg("viewer", "MODL")
    esp.TimHold("pyscript", "viewer")

    print("showing cylinder...")
    esp.SetModl(cylModl, esp.GetEsp("pyscript"))
    esp.TimMesg("viewer", "MODL")
    esp.TimHold("pyscript", "viewer")

    print("showing sphere...")
    esp.SetModl(sphModl, esp.GetEsp("pyscript"))
    esp.TimMesg("viewer", "MODL")
    esp.TimHold("pyscript", "viewer")
    
    print("showing cone...")
    esp.SetModl(conModl, esp.GetEsp("pyscript"))
    esp.TimMesg("viewer", "MODL")
    esp.TimHold("pyscript", "viewer")

esp.TimQuit("viewer")

esp.SetModl(cylModl, esp.GetEsp("pyscript"))
print("cylinder should be active MODL in ESP")

# boxModl and sphModl should get cleaned up after pyscript exits
