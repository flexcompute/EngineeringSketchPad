# Import pyCAPS module
import pyCAPS

import os

import argparse
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Eigen Value Pytest Example',
                                 prog = 'avl_EigenValue_PyTest.py',
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
problemName = str(args.workDir[0]) + "AVLEigenTest"
geometryScript = os.path.join("..","csmData","avlPlaneVanilla.csm")
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
fuse = {"groupName"    : "Fuselage",
        "numChord"     : 36}

wing = {"groupName"    : "Wing",
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24}

htail = {"numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 16}

vtail = {"numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 12}

avl.input.AVL_Surface = {"Fuse":fuse, "Wing": wing, "Htail": htail, "Vtail": vtail}


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
# Set mass properties
# -----------------------------------------------------------------
mass = 0.1773
x    =  0.02463
y    = 0.
z    = 0.2239
Ixx  = 1.350
Iyy  = 0.7509
Izz  = 2.095

avl.input.MassProp = {"Aircraft":{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}
avl.input.Gravity  = 32.18 * ft/s**2
avl.input.Density  = 0.002378 * slug/ft**3
avl.input.Velocity = 64.5396 * ft/s


# -----------------------------------------------------------------
# Set trim conditions (Eigen analysis should be done at a trimmed state)
# -----------------------------------------------------------------
Operation = {"Alpha":{"CL":0.2},
             "Stabilizer":{"Cm":0.0}
             }
avl.input.AVL_Operation = Operation

# -----------------------------------------------------------------
# Get Output Data from AVL
# These calls automatically run avl and access aimOutput data
# -----------------------------------------------------------------

EigenValues = avl.output.EigenValues
print ("EigenValues ", EigenValues)

# Print the CG location based on mass properties
print("Xcg", avl.output.Xcg)
print("Ycg", avl.output.Ycg)
print("Zcg", avl.output.Zcg)
