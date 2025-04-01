# Import pyCAPS module
import pyCAPS

import os

import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Cruciform Fuselage Pytest Example',
                                 prog = 'avl_CruciformFuselage_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = [""], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# -----------------------------------------------------------------
# Define units
# -----------------------------------------------------------------
m    = pyCAPS.Unit("meter")
kg   = pyCAPS.Unit("kg")
s    = pyCAPS.Unit("s")
K    = pyCAPS.Unit("Kelvin")
deg  = pyCAPS.Unit("degree")
ft   = pyCAPS.Unit("ft")
slug = pyCAPS.Unit("slug")

# -----------------------------------------------------------------
# Initialize Problem object
# -----------------------------------------------------------------
problemName = str(args.workDir[0]) + "AVLFuselageTest"
geometryScript = os.path.join("..","csmData","avlGlider.csm")
capsProblem = pyCAPS.Problem(problemName, capsFile=geometryScript, outLevel=args.outLevel)

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
## [loadAIM]
avl = capsProblem.analysis.create(aim = "avlAIM",
                                  name = "avl",
                                  unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

avl.input.Mach  = 0.5
avl.input.Alpha = 1.0*deg
avl.input.Beta  = 0.0*deg

# -----------------------------------------------------------------
# Set lifitng surface/slender body discretizations
# -----------------------------------------------------------------
fuseH = {"numChord"          : 24,
         "numSpanPerSection" : 1,
         "component"         : 1}

fuseW = {"numChord"          : 24,
         "numSpanPerSection" : 1,
         "component"         : 1}

wing = {"yMirror"      : True,
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24,
        "component"    : 1}

htail = {"yMirror"      : True,
         "numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 16,
         "component"    : 1}

vtail = {"numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 12,
         "component"    : 1}

avl.input.AVL_Surface = {"FuseH":fuseH, 
                         "FuseW":fuseW, 
                         "Wing": wing, 
                         "Htail": htail, 
                         "Vtail": vtail
                         }


# -----------------------------------------------------------------
# Set control surface parameters
# -----------------------------------------------------------------
Aileron = {"deflectionAngle" : 0*deg,
           "leOrTe" : 1,
           "controlGain" : -1.0,
           "hingeLine" : [0, 0, 0],
           "deflectionDup"  : -1.0}

Rudder = {"deflectionAngle" : 0*deg,
           "leOrTe" : 1,
           "controlGain" : 1.0,
           "hingeLine" : [0, 0, 0],
           "deflectionDup"  : 1.0}

Stabilizer = {"deflectionAngle" : 0*deg,
              "leOrTe" : 0,
              "controlGain" : 1.0,
              "hingeLine" : [0, 1, 0],
              "deflectionDup"  : 1.0}

avl.input.AVL_Control = {"Aileron": Aileron, "Stabilizer":Stabilizer, "Rudder":Rudder}


# -----------------------------------------------------------------
# Get Output Data from AVL
# These calls automatically run avl and access aimOutput data
# -----------------------------------------------------------------

## [output]
print ("CXtot  ", avl.output["CXtot" ].value)
print ("CYtot  ", avl.output["CYtot" ].value)
print ("CZtot  ", avl.output["CZtot" ].value)
print ("Cltot  ", avl.output["Cltot" ].value)
print ("Cmtot  ", avl.output["Cmtot" ].value)
print ("Cntot  ", avl.output["Cntot" ].value)
print ("Cl'tot ", avl.output["Cl'tot"].value)
print ("Cn'tot ", avl.output["Cn'tot"].value)
print ("CLtot  ", avl.output["CLtot" ].value)
print ("CDtot  ", avl.output["CDtot" ].value)
print ("CDvis  ", avl.output["CDvis" ].value)
print ("CLff   ", avl.output["CLff"  ].value)
print ("CYff   ", avl.output["CYff"  ].value)
print ("CDind  ", avl.output["CDind" ].value)
print ("CDff   ", avl.output["CDff"  ].value)
print ("e      ", avl.output["e"     ].value)
## [output]

# Check assertation
assert abs(avl.output["CLtot"].value-0.4868865261792446) <= 1E-4 
assert abs(avl.output["CDtot"].value-0.0088170924475462) <= 1E-4
assert abs(avl.output["e"    ].value-0.7753342232808564) <= 1E-4
