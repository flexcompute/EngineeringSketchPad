# Import pyCAPS and os module
import pyCAPS

import os

import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Pytest Example',
                                 prog = 'avl_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = "AVLAutoSpanAnalysisTest"
workDir = os.path.join(str(args.workDir[0]), workDir)

# -----------------------------------------------------------------
# Load CSM file and Change a design parameter - area in the geometry
# Any despmtr from the avlWing.csm file are available inside the pyCAPS script
# They are: thick, camber, area, aspect, taper, sweep, washout, dihedral
# -----------------------------------------------------------------

geometryScript = os.path.join("..","csmData","avlWings.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
myAnalysis = myProblem.analysis.create(aim = "avlAIM")

# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

myAnalysis.input.Mach  = 0.5
myAnalysis.input.Alpha = 1.0
myAnalysis.input.Beta  = 0.0

AVL_Surface = {}
for i in range(1,4):
    wing = {"groupName"    : "Wing"+str(i), # Notice Wing is the value for the capsGroup attribute
            "numChord"     : 12,
            "spaceChord"   : 1.0,
            "numSpanTotal" : 40,
            "spaceSpan"    : 1.0}

    AVL_Surface["Wing"+str(i)] = wing

myAnalysis.input.AVL_Surface = AVL_Surface

# -----------------------------------------------------------------
# Get Output Data from AVL
# These calls access aimOutput data
# -----------------------------------------------------------------
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

Cl = myAnalysis.output["CLtot"].value
Cd = myAnalysis.output["CDtot"].value
