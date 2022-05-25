###################################################################
#                                                                 #
# mitten1.py --- start with: serveESP ../data/basic/mitten1       #
#                    Tool->Pyscript   ../data/python/mitten1.py   #
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

# load the "mitten" TIM
esp.TimLoad("mitten", esp.GetEsp("pyscript"), "A")

# move block to left
esp.TimMesg("mitten", "xcent|-1|")
esp.TimMesg("mitten", "xcent|-1|")

# wait 5 seconds
esp.TimMesg("mitten", "countdown|5|")

# save the "mitten" results in the mitten1.csm file
esp.TimSave("mitten")

# load the "mitten" TIM
esp.TimLoad("mitten", esp.GetEsp("pyscript"), "B")

# move block to left
esp.TimMesg("mitten", "xcent|+1|")
esp.TimMesg("mitten", "xcent|+1|")

# wait 5 seconds
esp.TimMesg("mitten", "countdown|5|")

# save the "mitten" results in the mitten1.csm file
esp.TimSave("mitten")

print("Please choose \"File->Edit mitten1.csm\" to make your changes permanent")
