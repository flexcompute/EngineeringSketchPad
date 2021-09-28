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
geometryScript = os.path.join("..","csmData","avlWing.csm")
myProblem = pyCAPS.Problem(problemName, capsFile=geometryScript, outLevel=args.outLevel)

# -----------------------------------------------------------------
# Change a design parameter - area in the geometry
# Any despmtr from the avlWing.csm file are available inside the pyCAPS script
# They are: thick, camber, area, aspect, taper, sweep, washout, dihedral
# -----------------------------------------------------------------

myProblem.geometry.despmtr.area = 10.0

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
## [loadAIM]
myAnalysis = myProblem.analysis.create(aim = "avlAIM",
                                       name = "avl",
                                       unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

myAnalysis.input.Mach  = 0.5
myAnalysis.input.Alpha = 1.0*deg
myAnalysis.input.Beta  = 0.0*deg

wing = {"groupName"    : "Wing", # Notice Wing is the value for the capsGroup attribute
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24,
        "spaceSpan"    : 1.0}

myAnalysis.input.AVL_Surface = {"Wing": wing}

mass = 0.1773
x    =  0.02463
y    = 0.
z    = 0.2239
Ixx  = 1.350
Iyy  = 0.7509
Izz  = 2.095

myAnalysis.input.MassProp = {"Aircraft":{"mass":mass * kg, "CG":[x,y,z] * m, "massInertia":[Ixx, Iyy, Izz] * kg*m**2}}
myAnalysis.input.Gravity  = 32.18 * ft/s**2
myAnalysis.input.Density  = 0.002378 * slug/ft**3
myAnalysis.input.Velocity = 64.5396 * ft/s

# -----------------------------------------------------------------
# Get Output Data from AVL
# These calls automatically run avl and access aimOutput data
# -----------------------------------------------------------------

EigenValues = myAnalysis.output.EigenValues
print ("EigenValues ", EigenValues)
