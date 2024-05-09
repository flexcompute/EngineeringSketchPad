###################################################################
#                                                                 #
# areaRule.py -- adjust fuse:radius to satisfy Sears-Haack        #
#                                                                 #
#              can be executed in either of two ways:             #
#                  serveESP filename   (such as areaRule1.csm)    #
#                      Tool->Pyscript areaRule.py                 #
#              or                                                 #
#                  python areaRule.py                             #
#                      filename        (such as areaRule1.csm)    #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

from pyEGADS import egads
from pyOCSM  import ocsm
from pyOCSM  import esp

import math
import os

#------------------------------------------------------------------------------#

# callback function
def pyMesgCB(text):
    print(" ")
    print("======= in pyMesgCB =======")
    print("   ", text.decode())
    print("===========================")
    return

# make a semi-colon-separated string from a list
def makeString(array):
    out = ""
    for i in array:
        out += str(i) + ";"
    return out

#------------------------------------------------------------------------------#

# run quietly
ocsm.SetOutLevel(0)

# if we are running via serveESP, link to that MODL
try:
    modl = ocsm.Ocsm(esp.GetModl(esp.GetEsp("pyscript")))
    modl.RegMesgCB(pyMesgCB)
    print("==> getting MODL from ESP")

# an error means that we are probably running from the python prompt,
#    so get the filename from the user to create a new MODL
except ocsm.OcsmError:
    filename = ""
    while (".csm" not in filename):
        filename = input("Enter name of .csm file: ")
        if (not os.path.exists(filename)):
            print("\""+filename+"\" does not exist")
            filename = ""
    modl = ocsm.Ocsm(filename)
    print("==> making new MODL from \""+filename+"\"")

# check and build original MODL
modl.Check()
modl.Build(0, 0)

# get the pmtr indicies
ixloc   = modl.FindPmtr("fuse:xsect",    0, 0, 0)
iradius = modl.FindPmtr("fuse:radius",   0, 0, 0)
iarea   = modl.FindPmtr("aircraft:area", 0, 0, 0)

# get values from the MODL
nsect  = int(modl.GetValu(modl.FindPmtr("fuse:nsect",  0, 0, 0), 1, 1)[0])
length =     modl.GetValu(modl.FindPmtr("fuse:length", 0, 0, 0), 1, 1)[0]

xloc        = []
radius      = []
radius_lbnd = []
area        = []
for i in range(nsect):
    xloc.append(       modl.GetValu(ixloc,   i+1, 0)[0])
    radius.append(     modl.GetValu(iradius, i+1, 0)[0])
    radius_lbnd.append(modl.GetBnds(iradius, i+1, 0)[0])
    area.append(       modl.GetValu(iarea,   i+1, 0)[0])

    # find the maximum cross-sectional area
area_max   = 0
for i in range(nsect):
    if (area[i] > area_max):
        area_max = area[i]

# compute the sears-haack distribution
sears = []
for i in range(nsect):
    sears.append(area_max * math.pow(4 * xloc[i]/length * (1 - xloc[i]/length), 1.5))

# big iteration loop
niter = 10
for iter in range(niter+1):

    # show the current configuration
    if (iter > 0):
        esp.TimLoad("viewer", esp.GetEsp("pyscript"), "")
        esp.TimMesg("viewer", "MODL")
        esp.TimQuit("viewer")

    area   = []
    for i in range(nsect):
        area.append(modl.GetValu(iarea, i+1, 0)[0])
        
    # compute the rms error between the area and sears-haack
    rms = 0
    for i in range(nsect):
        rms += (area[i] - sears[i]) * (area[i] - sears[i])
    rms = round(math.sqrt(rms / nsect), 4)
    txt = "--> Iteration "+str(iter)+"  rms="+str(rms)
    print(txt)

    # show the area distribution
    esp.TimLoad("plotter", esp.GetEsp("pyscript"), "")
    esp.TimMesg("plotter", "new|"+txt+"|xloc|area|")
    esp.TimMesg("plotter", "add|"+makeString(xloc)+"|"+makeString(area)+"|k-+|")
    esp.TimMesg("plotter", "add|"+makeString(xloc)+"|"+makeString(sears)+"|g:|")
    esp.TimMesg("plotter", "show")
    esp.TimQuit("plotter")

    # check convergence
    if (rms < 0.01):
        print("--> converged")
        break
    elif (iter >= niter):
        print("--> out of iterations")
        break
                
    # compute new fuselage radii
    for i in range(nsect):
        radius[i] *= math.sqrt(sears[i] / area[i])
        if (radius[i] < radius_lbnd[i]):
            radius[i] = radius_lbnd[i]
        modl.SetValuD(iradius, i+1, 0, radius[i])

    # rebuild with the new radii
    modl.Build(0, 0)

print("==> saving final DESPMTRs in \"./areaRule.despmtrs\"")
modl.SaveDespmtrs("./areaRule.despmtrs")
