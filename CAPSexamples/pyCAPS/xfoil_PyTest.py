## [importModules]
# Import pyCAPS module
import pyCAPS

import os
import argparse
## [importModules]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'xFoil Pytest Example',
                                 prog = 'xfoil_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

## [localVariable]
# Create working directory variable
workDir = "xFoilAnalysisTest"
workDir = os.path.join(str(args.workDir[0]), workDir)
## [localVariable]

# Load CSM file
## [loadGeom]
geometryScript = os.path.join("..","csmData","airfoilSection.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [loadGeom]

## [setGeom]
# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.camber = 0.1
## [setGeom]

## [loadXFOIL]
# Load xfoil aim
xfoil = myProblem.analysis.create(aim = "xfoilAIM")
## [loadXFOIL]

##[setXFOIL]
# Set Mach number, Reynolds number
xfoil.input.Mach = 0.5
xfoil.input.Re   = 1.0e6

# Set custom AoA
xfoil.input.Alpha = [0.0, 3.0, 5.0, 7.0, 9.0, 11, 13, 14, 15.0]

# Set AoA seq
xfoil.input.Alpha_Increment = [1.0, 2.0, 0.10]

# Set custom Cl
xfoil.input.CL = 0.1

# Set Cl seq
xfoil.input.CL_Increment = [0.8, 3, .25]

# Append the polar file if it already exists - otherwise the AIM will delete the file
xfoil.input.Append_PolarFile = True
##[setXFOIL]

## [results]
# Retrieve Cl and Cd
Cl = xfoil.output.CL
print("Cl = ", Cl)

Cd = xfoil.output.CD
print("Cd = ", Cd)

# Angle of attack
Alpha = xfoil.output.Alpha
print("Alpha = ", Alpha)

# Transition location
TranX = xfoil.output.Transition_Top
print("Transition location = ", TranX)
## [results]

if args.noPlotData == False:
    # Import pyplot module
    try:
        from matplotlib import pyplot as plt

        # Plot first sweep angle
        plt.plot(Alpha,
                 Cl,
                 'o',
                  label = "Cl")

        plt.plot(Alpha,
                 Cd,
                 'x',
                 label="Cd")

        plt.legend(loc = "upper left")
        plt.title("Lift and Drag Coefficients for Various Angles of Attack\n (window must be closed to terminate Python script)")
        plt.xlabel("$Angle\ of\ Attack\ (^o)$")
        plt.ylabel("$C_L,\ C_D$")
        plt.show()
        
    except:
        print ("Unable to import matplotlib.pyplot module. Drag polar will not be plotted")

# Check assertation
assert abs(Cl[0]-1.2455) <= 1E-4 and abs(Cd[0]-0.01565) <= 1E-4
