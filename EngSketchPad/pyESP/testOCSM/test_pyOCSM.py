###################################################################
#                                                                 #
# test_pyOCSM --- test program to exercise pyOCSM                 #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

import os
import sys

from   pyEGADS import egads
from   pyOCSM  import ocsm

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

# get help on ocsm
skipHelp = 0;
if (len(sys.argv) > 1):
    if (sys.argv[1] == "-skipHelp"):
        skipHelp = 1
if (skipHelp == 0):
    help(ocsm)
else:
    print("skipping help")

print("\ntest 001: calling ocsm.Version()")
(imajor, iminor) = ocsm.Version()
print("    imajor :", imajor)
print("    iminor :", iminor)

# first MODL
print("\ntest 002: making modl1(sample.csm)")
modl1 = ocsm.Ocsm("sample.csm")

print("\ntest 003: calling modl1.LoadDict(sample.dict)")
modl1.LoadDict("sample.dict")

print("\ntest 004: calling modl1.RegMesgCB(pyMesgCB)")
modl1.RegMesgCB(pyMesgCB)

print("\ntest 005: calling modl1.Check()")
modl1.Check()

print("\ntest 006: calling modl1.GetFilelist()")
filelist = modl1.GetFilelist()
print("    filelist:", filelist)

print("\ntest 007: calling modl1.Build(0, 0)")
(builtTo, nbody, bodys) = modl1.Build(0, 0)
print("    builtTo:", builtTo);   assert (builtTo == 49  )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 008: calling modl1.Info()")
(nbrch, npmtr, nbody) = modl1.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 3 )

print("\ntest 009: calling modl1.PrintBodys()")
modl1.PrintBodys("")

print("\ntest 010: calling modl2=modl1.Copy()")
modl2 = modl1.Copy()
print("    modl2  :", modl2)

if (os.path.exists("sample.out")):
    os.remove("sample.out")

print("\ntest 011: calling modl2.PrintPmtrs(sample.out)")
f = open("sample.out", "w")
f.write("\ncalling modl2.PrintPmtrs\n")
f.close()
modl2.PrintPmtrs("sample.out")

print("\ntest 012: calling modl2.PrintAttrs(sample.out)")
f = open("sample.out", "a")
f.write("\ncalling modl2.PrintAttrs\n")
f.close()
modl2.PrintAttrs("sample.out")

print("\ntest 013: calling modl2.PrintBrchs(sample.out)")
f = open("sample.out", "a")
f.write("\ncalling modl2.PrintBrchs\n")
f.close()
modl2.PrintBrchs("sample.out")

print("\ntest 014: calling modl2.FindPmtr(rad)")
irad = modl2.FindPmtr("rad", 0, 0, 0)
print("    irad  :", irad);   assert (irad == 8)

print("\ntest 015: calling modl2.SetValuD(irad, 1, 1, 0.5)")
modl2.SetValuD(irad, 1, 1, 0.5)

print("\ntest 016: calling modl2.Build(0, 20)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl2.Build(0, 20)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (builtTo  == 49)
print("    nbody  :", nbody  );   assert (nbody    == 2 )
print("    bodys  :", bodys  );   assert (bodys[0] == 4 );   assert (bodys[1] == 3)

print("\ntest 017: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 4 )

print("\ntest 018: calling modl2.SaveDespmtrs(sample.despmtrs)")
modl2.SaveDespmtrs("sample.despmtrs")

print("\ntest 019: calling modl2.Free() --- should NOT remove EGADS context")
modl2.Free()

# this should cause an error
print("\ntest 020: calling modl2.Info() --- after modl2.Free")
try:
    (nbrch, npmtr, nbody) = modl2.Info()
except ocsm.OcsmError as error:
    print("    OcsmError rasied (as expected)")
else:
    raise Exception

print("\ntest 021: calling modl1.Info()")
(nbrch, npmtr, nbody) = modl1.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 3 )

print("\ntest 022: calling modl1.UpdateDespmtrs(sample.despmtrs)")
modl1.UpdateDespmtrs("sample.despmtrs")

print("\ntest 023: calling modl1.Build(0, 0)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl1.Build(0, 0)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (nbrch == 49)
print("    nbody  :", nbody  );   assert (npmtr == 62)
print("    bodys  :", bodys  );   assert (nbody == 0 )

print("\ntest 024: calling modl1.Info()")
(nbrch, npmtr, nbody) = modl1.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 4 )

print("\ntest 025: calling modl1.FindPmtr(oper)")
ioper = modl1.FindPmtr("oper", 0, 0, 0)
print("    ioper  :", ioper);   assert (ioper == 1)

print("\ntest 026: calling modl1.SetValu(ioper, 1, 1, 5-4)")
modl1.SetValu(ioper, 1, 1, "5-4")

print("\ntest 027: calling modl1.Build(0, 20)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl1.Build(0, 20)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (builtTo  == 49)
print("    nbody  :", nbody  );   assert (nbody    == 1 )
print("    bodys  :", bodys  );   assert (bodys[0] == 5 )

print("\ntest 028: calling modl1.Info()")
(nbrch, npmtr, nbody) = modl1.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 029: calling modl1.EvalExpr(dx*dy*dz)")
(value, dot, str) = modl1.EvalExpr("dx*dy*dz")
print("    value  :", value);   assert (value == 24)
print("    dot    :", dot  );   assert (dot   == 0 )
print("    str    :", str  );   assert (str   == "")

print("\ntest 030: calling modl1.EvalExpr(val2str(3.1415926,4))")
(value, dot, str) = modl1.EvalExpr("val2str(3.1415926,4)")
print("    value  :", value);   assert (value == 0 )
print("    dot    :", dot  );   assert (dot   == 0 )
print("    str    :", str  );   assert (str == "3.1416")

print("\ntest 031: calling modl1.Save(sample2.cpc)")
modl1.Save("sample2.cpc")

print("\ntest 032: calling PyOcsm.Load(sample2.cpc)")
modl2 = ocsm.Ocsm("sample2.cpc")

print("\ntest 033: calling modl2.LoadDict(sample.dict)")
modl2.LoadDict("sample.dict")

print("\ntest 034: calling modl2.Check()")
modl2.Check()

print("\ntest 035: calling modl2.Build(0, 0)")
(builtTo, nbody, bodys) = modl2.Build(0, 0)
print("    builtTo:", builtTo);   assert (nbrch == 49)
print("    nbody  :", nbody  );   assert (npmtr == 62)
print("    bodys  :", bodys  );   assert (nbody == 0 )
assert (nbody == 0)

print("\ntest 036: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 037: calling modl2.EvalExpr(@stack.size)")
(value, dot, str) = modl2.EvalExpr("@stack[1]")
print("    value  :", value);   assert (value == 5.0)
print("    dot    :", dot  );   assert (dot   == 0.0)
print("    str    :", str  );   assert (str   == "")

print("\ntest 038: calling modl2.PrintBrep(int(value), )")
modl2.PrintBrep(int(value), "")

print("\ntest 039: calling modl2.GetBrch(1)")
(type, bclass, actv, ichld, ileft, irite, narg, nattr) = modl2.GetBrch(1)
print("    type   :", type  );   assert (type   == 112)
print("    bclass :", bclass);   assert (bclass == 201)
print("    actv   :", actv  );   assert (actv   == 300)
print("    ichld  :", ichld );   assert (ichld  == 3  )
print("    ileft  :", ileft );   assert (ileft  == -1 )
print("    irite  :", irite );   assert (irite  == -1 )
print("    narg   :", narg  );   assert (narg   == 6  )
print("    nattr  :", nattr );   assert (nattr  == 0  )

print("\ntest 040: converting cylinder to sphere")
for ibrch in range(1, nbrch+1):
    (type, bclass, actv, ichld, ileft, irite, narg, nattr) = modl2.GetBrch(ibrch)

    if (modl2.GetText(type) == "cylinder"):
        print("    cylinder is ibrch", ibrch)

        print("    calling modl2.NewBrch(", ibrch, ", sphere x0+dx/2 y0+dy/2 z0+dz")
        modl2.NewBrch(ibrch, modl2.GetCode("sphere"), "<none>", 0,
                      "x0+dx/2", "y0+dy/2", "z0+dz", "rad", "", "", "", "", "")

        print("    calling modl2.SetName(", ibrch+1, ", sphere)")
        modl2.SetName(ibrch+1, "sphere")

        print("    calling modl2.SetAttr(", ibrch+1, ", _color $red)")
        modl2.SetAttr(ibrch+1, "_color", "$red")

        print("    calling modl2.SetAttr(", ibrch+1, ", _name $the_sphere)")
        modl2.SetAttr(ibrch+1, "_name",  "$the_sphere")

        print("    call modl2.GetAttr(", ibrch+1, ", _color)")
        avalue = modl2.GetAttr(ibrch+1, "_color")
        print("        avalue:", avalue);   assert (avalue == "$red")

        print("    getting all Attributes for ibrch=", ibrch+1)
        for iattr in range(1, 3):
            (aname, avalue) = modl2.RetAttr(ibrch+1, iattr)
            print("        ", aname, ":", avalue)

        print("    calling modl2.DelBrch(", ibrch, ")")
        modl2.DelBrch(ibrch)

print("\ntest 041: calling modl2.PrintBrchs()")
modl2.PrintBrchs("")

print("\ntest 042: calling modl2.Build(0, 0)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl2.Build(0, 0)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (builtTo == 49  )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 043: calling modl2.NewPmtr(Rad, DESPMTR, 1, 1)")
modl2.NewPmtr("Rad", ocsm.DESPMTR, 1, 1)
iRad = modl2.FindPmtr("Rad", 0, 0, 0)
print("    iRad   :", iRad);   assert (iRad == 19)

print("\ntest 044: calling modl2.SetValuD(iRad, 1, 1, 0.75)")
modl2.SetValuD(iRad, 1, 1, 0.75)

print("\ntest 045: calling modl2.SetBnds(iRad, 1, 1, -1, +2)")
modl2.SetBnds(iRad, 1, 1, -1.0, +2.0)

print("\ntest 046: calling modl2.SetValu(iRad, 1, 1, 2/3)")
modl2.SetValu(iRad, 1, 1, "2/3")

print("\ntest 047: calling modl2.GetBnds(iRad, 1, 1)")
(lbound, ubound) = modl2.GetBnds(iRad, 1, 1)
print("    lbound :", lbound);   assert (lbound == -1.0)
print("    ubound :", ubound);   assert (ubound ==  2.0)

print("\ntest 048: calling modl2.PrintPmtrs()")
modl2.PrintPmtrs("")

print("\ntest 049: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 19)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 050: converting 4th arg of sphere to Rad")
for ibrch in range(1, nbrch+1):
    if (modl2.GetName(ibrch) == "sphere"):

        print("    calling modl2.GetArg(", ibrch, ", 4)")
        (defn, value, dot) = modl2.GetArg(ibrch, 4)
        print("        defn :", defn );   assert (defn == "rad")
        print("        value:", value);   assert (abs(value-0.5) < TOL)
        print("        dot  :", dot  );   assert (dot == 0)

        print("    calling modl2.SetArg(", ibrch, ", 4, Rad")
        modl2.SetArg(ibrch, 4, "Rad")

        print("    calling modl2.GetArg(", ibrch, ", 4)")
        (defn, value, dot) = modl2.GetArg(ibrch, 4)
        print("        defn :", defn );   assert (defn == "Rad")
        print("        value:", value);   assert (abs(value-0.6666667) < TOL)
        print("        dot  :", dot  );   assert (dot == 0)

print("\ntest 051: calling modl2.NewBrch(nbrch, point 5  5  5")
modl2.NewBrch(nbrch, modl2.GetCode("point"), "<none>", 0,
              "5", "5", "5", "", "", "", "", "", "")

print("\ntest 052: calling modl2.Build(0, 0)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl2.Build(0, 20)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (builtTo  == 50)
print("    nbody  :", nbody  );   assert (nbody    == 2 )
print("    bodys  :", bodys  );   assert (bodys[0] == 6 );   assert (bodys[1] == 5)

print("\ntest 053: calling modl2.GetBody(bodys[0])")
(type, ichld, ileft, irite, vals, nnode, nedge, nface) = modl2.GetBody(bodys[0])
print("    type   :", type );   assert (type  == 111)
print("    ichld  :", ichld);   assert (ichld == 0  )
print("    ileft  :", ileft);   assert (ileft == -1 )
print("    irite  :", irite);   assert (irite == -1 )
print("    vals   :", vals )
print("    nnode  :", nnode);   assert (nnode == 1  )
print("    nedge  :", nedge);   assert (nedge == 0  )
print("    nface  :", nface);   assert (nface == 0  )

print("\ntest 054: calling modl2.SetBrch(", nbrch+1, ", SUPPRESSED)")
modl2.SetBrch(nbrch+1, ocsm.SUPPRESSED)

print("\ntest 055: calling modl2.Build(0, 0)")
ocsm.SetOutLevel(0)
(builtTo, nbody, bodys) = modl2.Build(0, 20)
ocsm.SetOutLevel(1)
print("    builtTo:", builtTo);   assert (builtTo  == 49)
print("    nbody  :", nbody  );   assert (nbody    == 1 )
print("    bodys  :", bodys  );   assert (bodys[0] == 5 )

print("\ntest 056: calling modl2.GetBody(bodys[0])")
(type, ichld, ileft, irite, vals, nnode, nedge, nface) = modl2.GetBody(bodys[0])
print("    type   :", type );   assert (type  == 143)
print("    ichld  :", ichld);   assert (ichld == 0  )
print("    ileft  :", ileft);   assert (ileft == 3  )
print("    irite  :", irite);   assert (irite == 4  )
print("    vals   :", vals )
print("    nnode  :", nnode);   assert (nnode == 11 )
print("    nedge  :", nedge);   assert (nedge == 18 )
print("    nface  :", nface);   assert (nface == 8  )

print("\ntest 057: calling modl2.DelBrch(", nbrch+1, ")")
modl2.DelBrch(nbrch+1)

print("\ntest 058: calling modl2.FindPmtr(foo, DESPMTR, 1, 1)")
modl2.FindPmtr("foo", ocsm.DESPMTR, 1, 1)

print("\ntest 059: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 20)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 060: alling modl2.DelPmtr(npmtr)")
modl2.DelPmtr(npmtr)

print("\ntest 061: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 19)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 062: making list of all DESPMTRs")
for ipmtr in range(1, npmtr+1):
    (type, nrow, ncol, name) = modl2.GetPmtr(ipmtr)

    if (type == ocsm.DESPMTR):
        print("   ", name, " is a", nrow, "*", ncol, "DESPMTR")

print("\ntest 063: making list of all OUTPMTRs")
for ipmtr in range(1, npmtr+1):
    (type, nrow, ncol, name) = modl2.GetPmtr(ipmtr)

    if (type == ocsm.OUTPMTR):
        print("  ", name, " is a", nrow, "*", ncol, "OUTPMTR")

print("\ntest 064: calling modl2.GetPmtr(iMyCG)")
iMyCG = modl2.FindPmtr("myCG", 0, 0, 0)
(type, nrow, ncol, name) = modl2.GetPmtr(iMyCG)

for irow in range(1, nrow+1):
    for icol in range(1, ncol+1):
        (value, dot) = modl2.GetValu(iMyCG, irow, icol)
        print("   ", name, "[", irow, ",", icol, "]: ", value, dot)

print("\ntest 065: calling modl2.GetValuS(ititle)")
ititle = modl2.FindPmtr("title", 0, 0, 0)
title  = modl2.GetValuS(ititle)
print("    title  :", title);   assert (title == "This_is_the_title")

print("\ntest 066: calling modl2.SetCsys(1, boxCsys, csysName)")
modl2.SetCsys(1, "boxCsys", "csysName")

print("\ntest 067: calling modl2.GetCsys(1, boxCsys)")
cvalue = modl2.GetCsys(1, "boxCsys")
print("    cvalue :", cvalue);   assert (cvalue == "csysName")

print("\ntest 068: calling modl2.RetCsys(1, 1)")
(cname, cvalue) = modl2.RetCsys(1, 1)
print("    cname  :", cname );   assert (cname  == "boxCsys" )
print("    cvalue :", cvalue);   assert (cvalue == "csysName")

print("\ntest 069: calling modl2.SetCsys(1, boxCsys, )")
modl2.SetCsys(1, "boxCsys", "")

print("\ntest 070: calling modl2.FindEnt(nbody, ocsm.FACE, [1,6,2])")
ient = modl2.FindEnt(nbody, ocsm.FACE, [1,6,1])
print("    ient   :", ient);   assert (ient == 3)

print("\ntest 071: getting coordinates on Face ient")
(xx, dot, str) = modl2.EvalExpr("x0+0.1")
print("    xx     :", xx);   assert (abs(xx-1.1) < TOL)
(yy, dot, str) = modl2.EvalExpr("y0+dy/2")
print("    yy     :", yy);   assert (abs(yy-3.5) < TOL)
(zz, dot, str) = modl2.EvalExpr("z0+dz")
print("    zz     :", zz);   assert (abs(zz-5.0) < TOL)

print("\ntest 072: calling modl2.GetUV(nbody, ocsm.FACE, ient, 1, [xx,yy,zz])")
(uu, vv) = modl2.GetUV(nbody, ocsm.FACE, ient, 1, [xx, yy, zz])
print("    uu     :", uu);   assert (abs(uu-0.1) < TOL)
print("    vv     :", vv);   assert (abs(vv-1.5) < TOL)

print("\ntest 073: calling modl2.GetXYZ(nbody, ocsm.FACE, ient, 2, [uu,vv,uu,vv])")
(xx, yy, zz, xxx, yyy, zzz) = modl2.GetXYZ(nbody, ocsm.FACE, ient, 2, [uu, vv, uu, vv])
print("    xxx    :", xxx);   assert (abs(xx-xxx) < TOL)
print("    yyy    :", yyy);   assert (abs(yy-yyy) < TOL)
print("    zzz    :", zzz);   assert (abs(zz-zzz) < TOL)

print("\ntest 074: calling modl2.GetNorm(nbody, ient, 1, [uu,vv])")
(normx, normy, normz) = modl2.GetNorm(nbody, ient, 1, [uu, vv])
print("    normx  :", normx);   assert (abs(normx    ) < TOL)
print("    normy  :", normy);   assert (abs(normy    ) < TOL)
print("    normz  :", normz);   assert (abs(normz-1.0) < TOL)

print("\ntest 075: calling modl2.FindPmtr(dx, 0, 0, 0)")
iz0 = modl2.FindPmtr("z0", 0, 0, 0)
print("    iz0    :", iz0);   assert (iz0 == 4)

print("\ntest 076: calling modl2.SetVel(iz0, 1, 1, 1/2)")
modl2.SetVel(iz0, 1, 1, "1/2");

print("\ntest 077: calling modl2.Build(0, 0)")
#ocsm.SetOutLevel(0)
modl2.Build(0, 0)
#ocsm.SetOutLevel(1)

print("\ntest 078: calling modl2.FindPmtr(myBbox, 0, 0, 0)")
iMyBbox = modl2.FindPmtr("myBbox", 0, 0, 0)
print("    iMyBbox:", iMyBbox);   assert (iMyBbox == 14)

values = [1.0, 2.0, 3.0, 5.0, 5.0, 5.0]
dots   = [0.0, 0.0, 0.5, 0.0, 0.0, 0.5]
indx   = 0
(type, nrow, ncol, name) = modl2.GetPmtr(iMyBbox)
for irow in range(1, nrow+1):
    for icol in range(1, ncol+1):
        (value, dot) = modl2.GetValu(iMyBbox, irow, icol)
        print("   ", name, "[", irow, ",", icol, "]=", value, dot)
        assert (abs(value-values[indx]) < TOL)
        assert (abs(dot  -dots[  indx]) < TOL)
        indx += 1

print("\ntest 079: calling modl2.SetVelD(0, 0, 0, 0)")
modl2.SetVelD(0, 0, 0, 0)

print("\ntest 080: calling modl2.SetVel(iz0, 1, 1, 1/2)")
modl2.SetVelD(iz0, 1, 1, 1)

print("\ntest 081: calling modl2.Build(0, 0)")
ocsm.SetOutLevel(0)
modl2.Build(0, 0)
ocsm.SetOutLevel(1)

print("\ntest 082: calling modl2.FindPmtr(myBbox, 0, 0, 0)")
iMyBbox = modl2.FindPmtr("myBbox", 0, 0, 0)
print("    iMyBbox:", iMyBbox);   assert (iMyBbox == 14)

values = [1.0, 2.0, 3.0, 5.0, 5.0, 5.0]
dots   = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
indx   = 0
(type, nrow, ncol, name) = modl2.GetPmtr(iMyBbox)
for irow in range(1, nrow+1):
    for icol in range(1, ncol+1):
        (value, dot) = modl2.GetValu(iMyBbox, irow, icol)
        print("   ", name, "[", irow, ",", icol, "]=", value, dot)
        assert (abs(value-values[indx]) < TOL)
        assert (abs(dot  -dots[  indx]) < TOL)
        indx += 1

print("\ntest 083: calling modl1.Info()")
(nbrch, npmtr, nbody) = modl1.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 62)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 084: calling modl2.GetBody(nbody)")
(type, ichld, ileft, irite, vals, nnode, nedge, nface) = modl2.GetBody(nbody)
print("    type   :", type );   assert (type  == 143)
print("    ichld  :", ichld);   assert (ichld == 0  )
print("    ileft  :", ileft);   assert (ileft == 3  )
print("    irite  :", irite);   assert (irite == 4  )
print("    vals   :", vals )
print("    nnode  :", nnode);   assert (nnode == 11 )
print("    nedge  :", nedge);   assert (nedge == 18 )
print("    nface  :", nface);   assert (nface == 8  )

print("\ntest 085: getting modl2.GetVel(nbody, ocsm.NODE, inode, 1, None)")
for inode in range(1, nnode+1):
    dxyz = modl2.GetVel(nbody, ocsm.NODE, inode, 1, None)
    print("    dxyz(inode=", inode, "): ", dxyz[0], dxyz[1], dxyz[2])
    assert (abs(dxyz[0]  ) < TOL)
    assert (abs(dxyz[1]  ) < TOL)
    assert (abs(dxyz[2]-1) < TOL)

print("\ntest 086: calling modl2.SetDtime(0.1)")
modl2.SetDtime(0.1)

print("\ntest 087: calling modl2.Build(0, 0)")
ocsm.SetOutLevel(0)
modl2.Build(0, 0)
ocsm.SetOutLevel(1)

print("\ntest 088: calling modl2.SetDtime(0.)")
modl2.SetDtime(0.)

print("\ntest 089: calling modl2.FindPmtr(myBbox, 0, 0, 0)")
iMyBbox = modl2.FindPmtr("myBbox", 0, 0, 0)
print("    iMyBbox:", iMyBbox);   assert (iMyBbox == 14)

values = [1.0, 2.0, 3.0, 5.0, 5.0, 5.0]
dots   = [0.0, 0.0, 1.0, 0.0, 0.0, 1.0]
indx   = 0
(type, nrow, ncol, name) = modl2.GetPmtr(iMyBbox)
for irow in range(1, nrow+1):
    for icol in range(1, ncol+1):
        (value, dot) = modl2.GetValu(iMyBbox, irow, icol)
        print("   ", name, "[", irow, ",", icol, "]=", value, dot)
        assert (abs(value-values[indx]) < TOL)
        assert (abs(dot  -dots[  indx]) < TOL)
        indx += 1

print("\ntest 090: calling modl2.SetVelD(0, 0, 0, 0)")
modl2.SetVelD(0, 0, 0, 0)

ix0 = modl2.FindPmtr("x0", 0, 0, 0)
iy0 = modl2.FindPmtr("y0", 0, 0, 0)
iz0 = modl2.FindPmtr("z0", 0, 0, 0)

print("\ntest 091: calling modl2.Perturb(3, [ix0, iy0, iz0], [1, 1, 1], [1, 1, 1], [1, 2, 3])")
modl2.Perturb(3, [ix0, iy0, iz0], [1, 1, 1], [1, 1, 1], [11, 22, 33])

print("\ntest 092; calling modl2.Perturb(0, 0, 0, 0, 0)")
modl2.Perturb(0, 0, 0, 0, 0)

print("\ntest 093: calling modl2.Info()")
(nbrch, npmtr, nbody) = modl2.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 49)
print("    npmtr  :", npmtr);   assert (npmtr == 63)
print("    nbody  :", nbody);   assert (nbody == 5 )

print("\ntest 094: calling modl2.GetBody(nbody)")
(type, ichld, ileft, irite, vals, nnode, nedge, nface) = modl2.GetBody(nbody)
print("    type   :", type );   assert (type  == 143)
print("    ichld  :", ichld);   assert (ichld == 0  )
print("    ileft  :", ileft);   assert (ileft == 3  )
print("    irite  :", irite);   assert (irite == 4  )
print("    vals   :", vals )
print("    nnode  :", nnode);   assert (nnode == 11 )
print("    nedge  :", nedge);   assert (nedge == 18 )
print("    nface  :", nface);   assert (nface == 8  )

print("\ntest 095: getting npnt from tessellations")
for iedge in range(1, nedge+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.EDGE, iedge)
    print("    edge", iedge, "has", npnt, "tessellation points")
for iface in range(1, nface+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.FACE, iface)
    print("    face", iface, "has", npnt, "tessellation points")

print("\ntest 096: calling modl2.NewBrch(nbrch, dump, <none>, 0, $sample.tess, 0, 0, 0, ...)")
modl2.NewBrch(nbrch, modl2.GetCode("dump"), "<none>", 0,
              "$sample.tess", "0", "0", "0", "", "", "", "", "")

print("\ntest 097: calling modl2.Build(0, 0)")
modl2.Build(0, 0)

print("\ntest 098: calling modl2.DelBrch(mbrch+1)")
modl2.DelBrch(nbrch+1)

print("\ntest 099: calling modl2.FindPmtr(tfact, 0, 0, 0)")
itfact = modl2.FindPmtr("tfact", 0, 0, 0)
print("    itfact :", itfact);   assert (itfact == 9)

print("\ntest 100: calling modl2.SetValuD(itfact, 1, 1, 1)")
modl2.SetValuD(itfact, 1, 1, 1.0)

print("\ntest 101: calling modl2.Build(-1, 0)   to ensure rebuild")
modl2.Build(-1, 0)

print("\ntest 102: getting npnt from tessellations")
for iedge in range(1, nedge+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.EDGE, iedge)
    print("    edge", iedge, "has", npnt, "tessellation points")
for iface in range(1, nface+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.FACE, iface)
    print("    face", iface, "has", npnt, "tessellation points")

print("\ntest 103: calling modl2.Updatetess(nbody, sample.tess)")
modl2.UpdateTess(nbody, "sample.tess")

print("\ntest 104: getting npnt from tessellations")
for iedge in range(1, nedge+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.EDGE, iedge)
    print("    edge", iedge, "has", npnt, "tessellation points")
for iface in range(1, nface+1):
    npnt = modl2.GetTessNpnt(nbody, ocsm.FACE, iface)
    print("    face", iface, "has", npnt, "tessellation points")

# print info about the Body and its parts
print("\nBody EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.BODY, 0)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.BODY     )
print("mtype :", mtype ); assert (mtype  == egads.SOLIDBODY)
ocsm.PrintEgo(theEgo)

print("\nTess EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.BODY, 1)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.TESSELLATION)
ocsm.PrintEgo(theEgo)

print("\nContext EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.BODY, 2)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.CONTXT)
ocsm.PrintEgo(theEgo)

print("\nNode[1] EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.NODE, 1)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.NODE   )
print("mtype :", mtype ); assert (mtype  == egads.NOMTYPE)
ocsm.PrintEgo(theEgo)

print("\nEdge[1] EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.EDGE, 1)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.EDGE   )
print("mtype :", mtype ); assert (mtype  == egads.TWONODE)
ocsm.PrintEgo(theEgo)

print("\nFace[1] EGO:")
theEgo = modl2.GetEgo(nbody, ocsm.FACE, 1)
(oclass, mtype, topRef, eprev, enext) = theEgo.getInfo()
print("oclass:", oclass); assert (oclass == egads.FACE    )
print("mtype :", mtype ); assert (mtype  == egads.SREVERSE)
ocsm.PrintEgo(theEgo)
del theEgo
del topRef
del eprev
del enext

print("\ntest 105: calling modl2.Free() --- should remove EGADS context")
modl2.Free()

# test multiple MODLs with inline files
print("\ntest 106: creating modl3 from inline.csm")
modl3 = ocsm.Ocsm("inline.csm")

print("\ntest 107: creating modl4 from inline.csm")
modl4 = ocsm.Ocsm("inline.csm")

print("\ntest 108: calling modl3.Build(0, 0)")
modl3.Build(0, 0)

print("\ntest 109: setting dx=5 for modl4")
modl4.SetValuD(modl4.FindPmtr("dx", 0, 0, 0), 1, 1, 5)

print("\ntest 110: calling modl4.Build(0, 0)")
modl4.Build(0, 0)

(xleft1, dot)  = modl3.GetValu(modl3.FindPmtr("xleft", 0, 0, 0), 1, 1)
print("xleft1 :", xleft1);   assert (xleft1 == 0)
(xrite1, dot) = modl3.GetValu(modl3.FindPmtr("xrite", 0, 0, 0), 1, 1)
print("xrite1 :", xrite1);   assert (xrite1 == 4)

(xleft2, dot) = modl4.GetValu(modl4.FindPmtr("xleft", 0, 0, 0), 1, 1)
print("xleft2 :", xleft2);   assert (xleft2 == 0)
(xrite2, dot) = modl4.GetValu(modl4.FindPmtr("xrite", 0, 0, 0), 1, 1)
print("xrite2 :", xrite2);   assert (xrite2 == 5)

print("\ntest 111: calling modl3.Free() --- should remove EGADS context")
modl3.Free()     # should get rid of tmp_OpenCSM_00 directory

print("\ntest 112: setting x0=1 for modl4")
modl4.SetValuD(modl4.FindPmtr("x0", 0, 0, 0), 1, 1, 1)

print("\ntest 113: calling modl4.Build(0, 0)")
modl4.Build(0, 0)

(xleft2, dot) = modl4.GetValu(modl4.FindPmtr("xleft", 0, 0, 0), 1, 1)
print("xleft2 :", xleft2);   assert (xleft2 == 1)
(xrite2, dot) = modl4.GetValu(modl4.FindPmtr("xrite", 0, 0, 0), 1, 1)
print("xrite2 :", xrite2);   assert (xrite2 == 6)

print("\ntest 114: calling modl5=modl4.Copy()")
modl5 = modl4.Copy();

print("\ntest 115: calling modl4.GetArg(8, 3)")
(defn, value, dot) = modl4.GetArg(8, 3)
print("defn  :", defn);   assert (defn == "$<<inline/0>>")

print("\ntest 116: calling modl4.PrintBrchs()")
modl4.PrintBrchs("")

print("\ntest 117: calling mod5.GetArg(8, 3)")
(defn, value, dot) = modl5.GetArg(8, 3)
print("defn  :", defn);   assert (defn == "$<<inline/0>>")

print("\ntest 118: calling modl5.PrintBrchs()")
modl5.PrintBrchs("")

# make a new modl
print("\ntest 119: making modl6()")
modl6 = ocsm.Ocsm("")

print("\ntest 120: calling modl6.RegMesgCB(pyMesgCB)")
modl6.RegMesgCB(pyMesgCB)

print("\ntest 121: calling modl6.RegSizeCB(pySizeCB)")
modl6.RegSizeCB(pySizeCB)

print("\ntest 122: calling modl6.FindPmtr(n)")
iN = modl6.FindPmtr("n", ocsm.CFGPMTR, 1, 1)
print("    iN    :", iN);   assert (iN == 1)

print("\ntest 123: calling SetValuD(n, 1, 1, 2)")
modl6.SetValuD(iN, 1, 1, 2)

print("\ntest 124: calling modl6.NewBrch(DIMENSION)")
modl6.NewBrch(0, modl6.GetCode("dimension"), "<none>", 0,
              "$dxyz", "1", "n", "", "", "", "", "", "")

print("\ntest 125: calling modl6.FindPmtr(dxyz)")
iDxyz = modl6.FindPmtr("dxyz", ocsm.DESPMTR, 1, 1)
print("    iDxyz :", iDxyz);    assert (iDxyz == 2)

print("\ntest 126: calling modl6.SetValuD(iDxyz, 1, 1, 1)")
print("(should trigger pySizeCB)")
modl6.SetValuD(iDxyz, 1, 1, 1)

print("\ntest 127: calling modl6.SetValuD(iDxyz, 1, 2, 2)")
modl6.SetValuD(iDxyz, 1, 2, 2)

print("\ntest 128: calling modl6.FindPmtr(myVol)")
iMyVol = modl6.FindPmtr("myVol", ocsm.OUTPMTR, 1, 1)
print("    iMyVol:", iMyVol);    assert (iMyVol == 3)

modl6.PrintPmtrs("")

print("\ntest 129: calling modl6.SetValuD(n, 1, 1, 3)")
print("(should trigger pySizeCB)")
modl6.SetValuD(iN, 1, 1, 3)

modl6.PrintPmtrs("")

print("\ntest 130: calling modl6.NewBrch(BOX)")
modl6.NewBrch(1, modl6.GetCode("box"), "<none>", 0,
              "0", "0", "0", "dxyz[1]", "dxyz[2]", "dxyz[3]", "", "", "")

print("\ntest 131: calling modl6,NewBrch(SET)")
modl6.NewBrch(2, modl6.GetCode("set"), "<none>", 0,
              "$myVol", "@volume", "", "", "", "", "", "", "")

print("\ntest 132: calling modl6.Build(0, 0)")
modl6.Build(0, 0)

print("\ntest 133: calling modl6.GetValu(myVol)")
(myVol, dot) = modl6.GetValu(iMyVol, 1, 1)
print("    myVol:", myVol);   assert (abs(myVol-4) < TOL)

print("\ntest 134: calling modl6.SetValuD(dxyz, 1, 3, 3)")
modl6.SetValuD(iDxyz, 1, 3, 3)

modl6.PrintPmtrs("")

print("\ntest 135: calling modl6.Build(0, 0)")
modl6.Build(0, 0)

print("\ntest 136: calling modl6.GetValu(myVol)")
(myVol, dot) = modl6.GetValu(iMyVol, 1, 1)
print("    myVol:", myVol);   assert (abs(myVol-6) < TOL)

print("\ntest 137: calling modl6.Free() --- should remove EGADS context")
modl6.Free()

print("\ntest 138: calling modl4.Info()")
(nbrch, npmtr, nbody) = modl4.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 13)
print("    npmtr  :", npmtr);   assert (npmtr == 53)
print("    nbody  :", nbody);   assert (nbody == 3 )

print("\ntest 139: adding \"UDPRIM applyTparams\" after BOX Branch")
for ibrch in range(1, nbrch+1):
    (type, bclass, actv, ichld, ileft, irite, narg, nattr) = modl4.GetBrch(ibrch)
    if (modl4.GetText(type) == "box"):
        modl4.NewBrch(ibrch, modl4.GetCode("udprim"), "<none>", 0,
                      "$$$/applyTparams", "$factor", "0.2", "", "", "", "", "", "")
modl4.PrintBrchs("")

print("\ntest 140: calling modl4.Build(0, 0)")
modl4.Build(0, 0)

print("\ntest 141: calling modl4.Info()")
(nbrch, npmtr, nbody) = modl4.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 25)
print("    npmtr  :", npmtr);   assert (npmtr == 53)
print("    nbody  :", nbody);   assert (nbody == 3 )

modl4.PrintBrchs("")

print("\ntest 142: calling modl4.Save(sample3.cpc)")
modl4.Save("sample3.cpc")

print("\ntest 143: making modl6(sample3.cpc)")
modl6 = ocsm.Ocsm("sample3.cpc")

modl6.PrintBrchs("")

print("\ntest 144: calling modl6.Check()")
modl6.Check()

print("\ntest 145: calling modl6.Build(0, 0)")
(builtTo, nbody, bodys) = modl6.Build(0, 0)
print("    builtTo:", builtTo);   assert (nbrch == 25  )
print("    nbody  :", nbody  );   assert (nbody == 0   )
print("    bodys  :", bodys  );   assert (bodys == None)
assert (nbody == 0)

print("\ntest 146: calling modl6.Info()")
(nbrch, npmtr, nbody) = modl6.Info()
print("    nbrch  :", nbrch);   assert (nbrch == 25)
print("    npmtr  :", npmtr);   assert (npmtr == 53)
print("    nbody  :", nbody);   assert (nbody == 3 )

print("\ntest 147: calling modl6.Free() --- should remove EGADS context")
modl6.Free()

print("\ntest 148: removing temporary sample* files")
os.remove("sample2.cpc")
os.remove("sample3.cpc")
os.remove("sample.despmtrs")
os.remove("sample.out")
os.remove("sample.tess")

print("\ntest 149: making modl7()")
modl7 = ocsm.Ocsm("")

print("\ntest 150: adding DESPMTRs")
modl7.NewPmtr("xbeg", ocsm.DESPMTR, 1, 1)
modl7.NewPmtr("ybeg", ocsm.DESPMTR, 1, 1)
modl7.NewPmtr("zbeg", ocsm.DESPMTR, 1, 1)
modl7.NewPmtr("xend", ocsm.DESPMTR, 1, 1)
modl7.NewPmtr("yend", ocsm.DESPMTR, 1, 1)
modl7.NewPmtr("zend", ocsm.DESPMTR, 1, 1)

modl7.NewPmtr("myLen", ocsm.OUTPMTR, 1, 1)

print("\ntest 151: add Branches")
print("    calling modl7.NewBrch(", 0, ") skbeg   xbeg  ybeg  zbeg  0")
modl7.NewBrch(0, modl7.GetCode("skbeg"), "<none>", 0,
              "xbeg", "ybeg", "zbeg", "0", "", "", "", "", "")
print("    calling modl7.NewBrch(", 1, ") linseg  xend  yend  zend")
modl7.NewBrch(1, modl7.GetCode("linseg"), "<none>", 0,
              "xend", "yend", "zend", "", "", "", "", "", "")
print("    calling modl7.NewBrch(", 2, ") skbeg   0")
modl7.NewBrch(2, modl7.GetCode("skend"),  "<none>", 0,
              "0", "", "", "", "", "", "", "", "")
print("    calling modl7.NewBrch(", 3, ") set     myLen  @length")
modl7.NewBrch(3, modl7.GetCode("set"), "<none>", 0,
              "$myLen", "@length", "", "", "", "", "", "", "")

print("\ntest 152: setting DESPMTR values")
modl7.SetValuD(modl7.FindPmtr("xbeg", 0, 0, 0), 1, 1, -1)
modl7.SetValuD(modl7.FindPmtr("ybeg", 0, 0, 0), 1, 1, -2)
modl7.SetValuD(modl7.FindPmtr("zbeg", 0, 0, 0), 1, 1, -3)
modl7.SetValuD(modl7.FindPmtr("xend", 0, 0, 0), 1, 1, +1)
modl7.SetValuD(modl7.FindPmtr("yend", 0, 0, 0), 1, 1, +2)
modl7.SetValuD(modl7.FindPmtr("zend", 0, 0, 0), 1, 1, +3)
modl7.PrintPmtrs("")

print("\ntest 153: calling modl7.Build(0, 0)")
(builtTo, nbody, bodys) = modl7.Build(0, 0)
print("    builtTo:", builtTo);   assert (builtTo == 4   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 154: calling modl7.GetTessVel(3, NODE, 1)")
vels = modl7.GetTessVel(3, ocsm.NODE, 1)
print("    node1:xvel :", vels[0]);   assert (abs(vels[0]-0) < TOL)
print("    node1:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node1:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 155: calling modl7.GetTessVel(3, NODE, 2)")
vels = modl7.GetTessVel(3, ocsm.NODE, 2)
print("    node2:xvel :", vels[0]);   assert (abs(vels[0]-0) < TOL)
print("    node2:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node2:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 156: calling modl7.GetTessVel(3, EDGE, 1)")
vels = modl7.GetTessVel(3, ocsm.EDGE, 1)
print("    beg.xvel :", vels[ 0]);   assert (abs(vels[ 0]-0) < TOL)
print("    beg.yvel :", vels[ 1]);   assert (abs(vels[ 1]-0) < TOL)
print("    beg.zvel :", vels[ 2]);   assert (abs(vels[ 2]-0) < TOL)
print("    end.xvel :", vels[-3]);   assert (abs(vels[-3]-0) < TOL)
print("    end.yvel :", vels[-2]);   assert (abs(vels[-2]-0) < TOL)
print("    end.zvel :", vels[-1]);   assert (abs(vels[-1]-0) < TOL)

print("\ntest 157: call modl.GetValu(myLen)")
(value, dot) = modl7.GetValu(modl7.FindPmtr("myLen", 0, 0, 0), 1, 1)
print("    value:", value);   assert (abs(value-7.48332) < TOL)
print("    dot  :", dot  );   assert (abs(dot  -0      ) < TOL)

print("\ntest 158: calling modl7.SetVelD(xbeg, +1)")
modl7.SetVelD(modl7.FindPmtr("xbeg", 0, 0, 0), 1, 1, +1)

print("\ntest 159: calling modl7.SetVelD(zend, -1)")
modl7.SetVelD(modl7.FindPmtr("zend", 0, 0, 0), 1, 1, -1)

print("\ntest 160: calling modl7.Build(0, 0)")
(builtTo, nbody, bodys) = modl7.Build(0, 0)
print("    builtTo:", builtTo);   assert (builtTo == 4   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 161: call modl.GetValu(myLen)")
(value, dot) = modl7.GetValu(modl7.FindPmtr("myLen", 0, 0, 0), 1, 1)
print("    value:", value);   assert (abs(value-7.48332) < TOL)
print("    dot  :", dot  );   assert (abs(dot  +1.06904) < TOL)

print("\ntest 162: calling modl7.GetTessVel(3, NODE, 1)")
vels = modl7.GetTessVel(3, ocsm.NODE, 1)
print("    node1:xvel :", vels[0]);   assert (abs(vels[0]-1) < TOL)
print("    node1:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node1:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 163: calling modl7.GetTessVel(3, NODE, 2)")
vels = modl7.GetTessVel(3, ocsm.NODE, 2)
print("    node2:xvel :", vels[0]);   assert (abs(vels[0]-0) < TOL)
print("    node2:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node2:zvel :", vels[2]);   assert (abs(vels[2]+1) < TOL)

print("\ntest 164: calling modl7.GetTessVel(3, EDGE, 1)")
vels = modl7.GetTessVel(3, ocsm.EDGE, 1)
print("    beg.xvel :", vels[ 0]);   assert (abs(vels[ 0]-1) < TOL)
print("    beg.yvel :", vels[ 1]);   assert (abs(vels[ 1]-0) < TOL)
print("    beg.zvel :", vels[ 2]);   assert (abs(vels[ 2]-0) < TOL)
print("    end.xvel :", vels[-3]);   assert (abs(vels[-3]-0) < TOL)
print("    end.yvel :", vels[-2]);   assert (abs(vels[-2]-0) < TOL)
print("    end.zvel :", vels[-1]);   assert (abs(vels[-1]+1) < TOL)

print("\ntest 165: calling modl7.Free() --- should remove EGADS context")
modl7.Free()

print("\ntest 166: making modl8()")
modl8 = ocsm.Ocsm("")

print("\ntest 167: adding DESPMTRs")
modl8.NewPmtr("xbeg", ocsm.DESPMTR, 1, 1)
modl8.NewPmtr("ybeg", ocsm.DESPMTR, 1, 1)
modl8.NewPmtr("zbeg", ocsm.DESPMTR, 1, 1)
modl8.NewPmtr("xlen", ocsm.DESPMTR, 1, 1)

modl8.NewPmtr("myLen", ocsm.OUTPMTR, 1, 1)

print("\ntest 168: add Branches (backwards)")
print("    calling modl8.NewBrch(", 0, ") set     myLen  @length")
modl8.NewBrch(0, modl8.GetCode("set"), "<none>", 0,
              "$myLen", "@length", "", "", "", "", "", "", "")
print("    calling modl8.NewBrch(", 0, ") translate xbeg+xlen/2  ybeg  zbeg")
modl8.NewBrch(0, modl8.GetCode("translate"), "<none>", 0,
              "xbeg+xlen/2", "ybeg", "zbeg", "", "", "", "", "", "")
print("    calling modl8.NewBrch(", 0, ") udprim    box          dx    xlen")
modl8.NewBrch(0, modl8.GetCode("udprim"), "<none>", 0,
              "$box", "$dx", "xlen", "", "", "", "", "", "")

print("\ntest 169: setting DESPMTR values")
modl8.SetValuD(modl8.FindPmtr("xbeg", 0, 0, 0), 1, 1, -1)
modl8.SetValuD(modl8.FindPmtr("ybeg", 0, 0, 0), 1, 1, -2)
modl8.SetValuD(modl8.FindPmtr("zbeg", 0, 0, 0), 1, 1, -3)
modl8.SetValuD(modl8.FindPmtr("xlen", 0, 0, 0), 1, 1,  4)
modl8.PrintPmtrs("")

print("\ntest 170: calling modl8.Build(0, 0)")
(builtTo, nbody, bodys) = modl8.Build(0, 0)
print("    builtTo:", builtTo);   assert (builtTo == 3   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 171: calling modl8.GetTessVel(2, NODE, 1)")
vels = modl8.GetTessVel(2, ocsm.NODE, 1)
print("    node1:xvel :", vels[0]);   assert (abs(vels[0]-0) < TOL)
print("    node1:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node1:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 172: calling modl8.GetTessVel(2, NODE, 2)")
vels = modl8.GetTessVel(2, ocsm.NODE, 2)
print("    node2:xvel :", vels[0]);   assert (abs(vels[0]-0) < TOL)
print("    node2:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node2:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 173: calling modl8.GetTessVel(2, EDGE, 1)")
vels = modl8.GetTessVel(2, ocsm.EDGE, 1)
print("    beg.xvel :", vels[ 0]);   assert (abs(vels[ 0]-0) < TOL)
print("    beg.yvel :", vels[ 1]);   assert (abs(vels[ 1]-0) < TOL)
print("    beg.zvel :", vels[ 2]);   assert (abs(vels[ 2]-0) < TOL)
print("    end.xvel :", vels[-3]);   assert (abs(vels[-3]-0) < TOL)
print("    end.yvel :", vels[-2]);   assert (abs(vels[-2]-0) < TOL)
print("    end.zvel :", vels[-1]);   assert (abs(vels[-1]-0) < TOL)

print("\ntest 174: call modl.GetValu(myLen)")
(value, dot) = modl8.GetValu(modl8.FindPmtr("myLen", 0, 0, 0), 1, 1)
print("    value:", value);   assert (abs(value-4) < TOL)
print("    dot  :", dot  );   assert (abs(dot  -0) < TOL)

print("\ntest 175: calling modl8.SetVelD(xbeg, +1)")
modl8.SetVelD(modl8.FindPmtr("xbeg", 0, 0, 0), 1, 1, +1)

print("\ntest 176: calling modl8.SetVelD(xlen, +1)")
modl8.SetVelD(modl8.FindPmtr("xlen", 0, 0, 0), 1, 1, +1)

print("\ntest 177: calling modl8.Build(0, 0)")
(builtTo, nbody, bodys) = modl8.Build(0, 0)
print("    builtTo:", builtTo);   assert (builtTo == 3   )
print("    nbody  :", nbody  );   assert (nbody   == 0   )
print("    bodys  :", bodys  );   assert (bodys   == None)

print("\ntest 178: call modl.GetValu(myLen)")
(value, dot) = modl8.GetValu(modl8.FindPmtr("myLen", 0, 0, 0), 1, 1)
print("    value:", value);   assert (abs(value-4) < TOL)
print("    dot  :", dot  );   assert (abs(dot  -1) < TOL)

print("\ntest 179: calling modl8.GetTessVel(3, NODE, 1)")
vels = modl8.GetTessVel(2, ocsm.NODE, 1)
print("    node1:xvel :", vels[0]);   assert (abs(vels[0]-1) < TOL)
print("    node1:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node1:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 180: calling modl8.GetTessVel(3, NODE, 2)")
vels = modl8.GetTessVel(2, ocsm.NODE, 2)
print("    node2:xvel :", vels[0]);   assert (abs(vels[0]-2) < TOL)
print("    node2:yvel :", vels[1]);   assert (abs(vels[1]-0) < TOL)
print("    node2:zvel :", vels[2]);   assert (abs(vels[2]-0) < TOL)

print("\ntest 181: calling modl8.GetTessVel(3, EDGE, 1)")
vels = modl8.GetTessVel(2, ocsm.EDGE, 1)
print("    beg.xvel :", vels[ 0]);   assert (abs(vels[ 0]-1) < TOL)
print("    beg.yvel :", vels[ 1]);   assert (abs(vels[ 1]-0) < TOL)
print("    beg.zvel :", vels[ 2]);   assert (abs(vels[ 2]-0) < TOL)
print("    end.xvel :", vels[-3]);   assert (abs(vels[-3]-2) < TOL)
print("    end.yvel :", vels[-2]);   assert (abs(vels[-2]-0) < TOL)
print("    end.zvel :", vels[-1]);   assert (abs(vels[-1]-0) < TOL)

print("\ntest 182: get a reference to a Body in modl8")
(nbrch, npmtr, nbody) = modl8.Info()
newBody = modl8.GetEgo(nbody, ocsm.BODY, 0)

print("\ntest 183: calling newBody.getInfo()")
(oclass, mtype, topRed, eprev, enext) = newBody.getInfo()
print("oclass:", oclass);   assert (oclass == egads.BODY    )
print("mtype :", mtype );   assert (mtype  == egads.WIREBODY)

print("\ntest 184: calling modl8.Free() --- should remove EGADS context")
modl8.Free()

print("\ntest 185: calling modl5.Free() --- should NOT remove EGADS context")
modl5.Free()

print("\ntest 186: making modl3()")
modl3 = ocsm.Ocsm("")

print("\ntest 187: calling modl3.NewBrch(udprim,naca)")
modl3.NewBrch(0, modl3.GetCode("udprim"), "<none>", 0,
              "naca", "$thickness", "0.12", "$camber", "0.04", "", "", "", "")

print("\ntest 188: calling modl3.NewBrch(extrude)")
modl3.NewBrch(1, modl3.GetCode("extrude"), "<none>", 0,
              "0", "0", "3", "", "", "", "", "", "")

print("\ntest 189: calling modl3.NewBrch(udprim,printBbox)")
modl3.NewBrch(2, modl3.GetCode("udprim"), "<none>", 0,
              "$printBbox", "", "", "", "", "", "", "", "")

print("\ntest 190: calling modl3.NewBrch(message)")
modl3.NewBrch(3, modl3.GetCode("message"), "<none>", 0,
              "$area=+@area", "$_", "", "", "", "", "", "", "")

print("\ntest 191: calling modl3.PrintBrchs()")
modl3.PrintBrchs("")

print("\ntest 192: calling modl3.Build(0,0)")
modl3.Build(0, 0)

print("\ntest 193: calling modl3,GetTessNpnt(nbody, FACE, 1)")
(nbrch, npmtr, nbody) = modl3.Info()
npnt = modl3.GetTessNpnt(nbody, ocsm.FACE, 1)
print("    npnt     :", npnt);   assert (npnt == 2517)

print("\ntest 194: calling modl3.GetEgo(BODY)")
theBody = modl3.GetEgo(nbody, ocsm.BODY, 0)

print("\ntest 195: calling MakeTessBody")
theTess = theBody.makeTessBody((0.01, 0.01, 5))

print("\ntest 196: calling modl3.SetEgo(TESS)")
modl3.SetEgo(nbody, 1, theTess)

print("\ntest 197: calling modl3.GetTessNpnt(nbody, FACE, 1)")
npnt = modl3.GetTessNpnt(nbody, ocsm.FACE, 1)
print("    npnt     :", npnt);   assert (npnt == 73359)

print("\ntest 198: calling modl3.Free() --- should remove EGADS context")
modl3.Free()

print("\nNOTE:")
print("    GetSketch    is not tested")
print("    SolveSketch  is not tested")
print("    UpdateSketch is not tested")
print("    SetEgg       is not tested")

print("\nAt end:")
print("    EGADS context associated with modl1    should be removed")
print("    EGADS context associated with modl4/5  should be removed")

print("\ntest_ocsm finished successfully\n")
