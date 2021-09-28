## [import]
# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse
## [import]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FRICTION Pytest Example',
                                 prog = 'friction_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

## [localVariable]
# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "FrictionAnalysisTest")
## [localVariable]

# -----------------------------------------------------------------
# Load CSM file
# -----------------------------------------------------------------
## [geometry]
geometryScript = os.path.join("..","csmData","frictionWingTailFuselage.csm")
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
myAnalysis = myProblem.analysis.create( aim = "frictionAIM" )
## [loadAIM]

# -----------------------------------------------------------------
# Set new Mach/Alt parameters
# -----------------------------------------------------------------
print ("Setting Mach & Altitude Values")
## [setInputs]

myAnalysis.input.Mach = [0.5, 1.5]

# Note: friction wants kft (defined in the AIM) - Automatic unit conversion to kft
myAnalysis.input.Altitude = [9000, 18200.0]*pyCAPS.Unit("m")
## [setInputs]

# -----------------------------------------------------------------
# Get Output Data from Friction (execution is automatic just-in-time)
# -----------------------------------------------------------------
## [output]
Cdtotal = myAnalysis.output.CDtotal
CdForm  = myAnalysis.output.CDform
CdFric  = myAnalysis.output.CDfric
## [output]

print("Total drag =", Cdtotal )
print("Form drag =", CdForm)
print("Friction drag =", CdFric)
