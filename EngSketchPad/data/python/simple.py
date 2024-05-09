###################################################################
#                                                                 #
# simple.py --- example python script                             #
#                                                                 #
# this can be executed with either:                               #
#    python simple.py                                             #
# or:                                                             #
#    serveESP simple.py                                           #
# or:                                                             #
#    serveESP                                                     #
#        Tools->Pyscript  simple.py                               #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

# if running from "serveESP", then this will allow the user to view the modl once the script completes
modl = ocsm.Ocsm(esp.GetModl(esp.GetEsp("pyscript")))

# if running from "python simple.py", create an empty modl
try:
    ipmtr = modl.FindPmtr("@vol", 0, 0, 0)
except:
    modl = ocsm.Ocsm("")

# add SPHERE and DUMP branches
modl.NewBrch(0, modl.GetCode("sphere"), "<none>", 0,
             "0", "0", "0", "1", "0", "0", "0", "0", "0")

modl.NewBrch(1, modl.GetCode("dump"), "<none>", 0,
             "simple.step", "0", "0", "0", "0", "0", "0", "0", "0");

# build the modl
(builtTo, nbody, bodys) = modl.Build(0, 0)


