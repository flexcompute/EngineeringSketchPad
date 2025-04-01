# Import pyCAPS and os module
## [import]
import pyCAPS
import os
import argparse
## [import]

## [argparse]
# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Cart3D Pytest Example',
                                 prog = 'cart3d_PyTest',
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
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)
## [geometry]

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
## [loadAIM]
cart3d = capsProblem.analysis.create(aim = "cart3dAIM")
## [loadAIM]
# -----------------------------------------------------------------
# Also available are all aimInput values
# Set new Mach/Alt parameters
# -----------------------------------------------------------------

## [setInputs]
cart3d.input.Mach      = 0.5
cart3d.input.alpha     = 2.0
cart3d.input.maxCycles = 10
cart3d.input.nDiv      = 2
cart3d.input.maxR      = 11
cart3d.input.y_is_spanwise = True
## [setInputs]

# -----------------------------------------------------------------
# Cart3D auto-executes
# -----------------------------------------------------------------

## [output]
# List entire outputs
print ("entire C_A  " + str(cart3d.output.C_A))
print ("entire C_Y  " + str(cart3d.output.C_Y))
print ("entire C_N  " + str(cart3d.output.C_N))
print ("entire C_D  " + str(cart3d.output.C_D))
print ("entire C_S  " + str(cart3d.output.C_S))
print ("entire C_L  " + str(cart3d.output.C_L))
print ("entire C_l  " + str(cart3d.output.C_l))
print ("entire C_m  " + str(cart3d.output.C_m))
print ("entire C_n  " + str(cart3d.output.C_n))
print ("entire C_M_x  " + str(cart3d.output.C_M_x))
print ("entire C_M_y  " + str(cart3d.output.C_M_y))
print ("entire C_M_z  " + str(cart3d.output.C_M_z))

# List Wing component outputs
Wing = cart3d.dynout.Wing
print ("Wing C_A  " + str(Wing["C_A"]))
print ("Wing C_Y  " + str(Wing["C_Y"]))
print ("Wing C_N  " + str(Wing["C_N"]))
print ("Wing C_D  " + str(Wing["C_D"]))
print ("Wing C_S  " + str(Wing["C_S"]))
print ("Wing C_L  " + str(Wing["C_L"]))
print ("Wing C_l  " + str(Wing["C_l"]))
print ("Wing C_m  " + str(Wing["C_m"]))
print ("Wing C_n  " + str(Wing["C_n"]))
print ("Wing C_M_x  " + str(Wing["C_M_x"]))
print ("Wing C_M_y  " + str(Wing["C_M_y"]))
print ("Wing C_M_z  " + str(Wing["C_M_z"]))
## [output]
