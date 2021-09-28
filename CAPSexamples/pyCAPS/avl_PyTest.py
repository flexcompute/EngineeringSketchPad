# Import pyCAPS and os module
## [import]
import pyCAPS

import os
## [import]

import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Pytest Example',
                                 prog = 'avl_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# -----------------------------------------------------------------
# Load CSM file and Change a design parameter - area in the geometry
# Any despmtr from the avlWing.csm file are available inside the pyCAPS script
# They are: thick, camber, area, aspect, taper, sweep, washout, dihedral
# -----------------------------------------------------------------

# Create working directory variable
## [localVariable]
workDir = "AVLAnalysisTest"
## [localVariable]
workDir = os.path.join(str(args.workDir[0]), workDir)

## [geometry]
geometryScript = os.path.join("..","csmData","avlWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

myProblem.geometry.despmtr.area = 10.0
## [geometry]

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
## [loadAIM]
myAnalysis = myProblem.analysis.create(aim = "avlAIM", name = "avl")
## [loadAIM]
# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

## [setInputs]
myAnalysis.input.Mach  = 0.5
myAnalysis.input.Alpha = 1.0
myAnalysis.input.Beta  = 0.0

wing = {"groupName"    : "Wing", # Notice Wing is the value for the capsGroup attribute
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24,
        "spaceSpan"    : 1.0}

myAnalysis.input.AVL_Surface = {"Wing": wing}
## [setInputs]

# -----------------------------------------------------------------
# Get Output Data from AVL
# These calls access aimOutput data
# -----------------------------------------------------------------
## [output]
print ("CXtot  ", myAnalysis.output["CXtot" ].value)
print ("CYtot  ", myAnalysis.output["CYtot" ].value)
print ("CZtot  ", myAnalysis.output["CZtot" ].value)
print ("Cltot  ", myAnalysis.output["Cltot" ].value)
print ("Cmtot  ", myAnalysis.output["Cmtot" ].value)
print ("Cntot  ", myAnalysis.output["Cntot" ].value)
print ("Cl'tot ", myAnalysis.output["Cl'tot"].value)
print ("Cn'tot ", myAnalysis.output["Cn'tot"].value)
print ("CLtot  ", myAnalysis.output["CLtot" ].value)
print ("CDtot  ", myAnalysis.output["CDtot" ].value)
print ("CDvis  ", myAnalysis.output["CDvis" ].value)
print ("CLff   ", myAnalysis.output["CLff"  ].value)
print ("CYff   ", myAnalysis.output["CYff"  ].value)
print ("CDind  ", myAnalysis.output["CDind" ].value)
print ("CDff   ", myAnalysis.output["CDff"  ].value)
print ("e      ", myAnalysis.output["e"     ].value)
## [output]

## [strip forces]
StripForces = myAnalysis.output.StripForces
print ("StripForces ")
for name in StripForces.keys():
    print (name, " cl = ", StripForces[name]["cl"])
    print (name, " cd = ", StripForces[name]["cd"])
## [strip forces]

## [sensitivity]
#sensitivity = myAnalysis.output["CLtot"].deriv("Alpha")
## [sensitivity]

Cl = myAnalysis.output.CLtot
Cd = myAnalysis.output.CDtot

# Check assertation
assert abs(Cl-0.30126) <= 1E-4 and abs(Cd-0.00465) <= 1E-4
