# Import pyCAPS and os module
## [import]
import pyCAPS
import os
import argparse
## [import]

## [argparse]
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Cart3D Pytest Example',
                                 prog = 'cart_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "cart3dTest")
## [argparse]

# -----------------------------------------------------------------
# Load CSM file and Change a design parameter - area in the geometry
# Any despmtr from the avlWing.csm file are available inside the pyCAPS script
# They are: thick, camber, area, aspect, taper, sweep, washout, dihedral
# -----------------------------------------------------------------

## [geometry]
geometryScript = os.path.join("..","csmData","cfd_airfoilSection.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [geometry]

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
## [loadAIM]
myAnalysis = myProblem.analysis.create(aim = "cart3dAIM")
## [loadAIM]
# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

## [setInputs]
myAnalysis.input.Mach      = 0.95
myAnalysis.input.alpha     = 2.0
myAnalysis.input.maxCycles = 10
myAnalysis.input.nDiv      = 6
myAnalysis.input.maxR      = 12
## [setInputs]

# -----------------------------------------------------------------
# Cart3D auto-executes
# -----------------------------------------------------------------

## [output]
print ("C_A  " + str(myAnalysis.output.C_A))
print ("C_Y  " + str(myAnalysis.output.C_Y))
print ("C_N  " + str(myAnalysis.output.C_N))
print ("C_D  " + str(myAnalysis.output.C_D))
print ("C_S  " + str(myAnalysis.output.C_S))
print ("C_L  " + str(myAnalysis.output.C_L))
print ("C_l  " + str(myAnalysis.output.C_l))
print ("C_m  " + str(myAnalysis.output.C_m))
print ("C_n  " + str(myAnalysis.output.C_n))
print ("C_M_x  " + str(myAnalysis.output.C_M_x))
print ("C_M_y  " + str(myAnalysis.output.C_M_y))
print ("C_M_z  " + str(myAnalysis.output.C_M_z))
## [output]
