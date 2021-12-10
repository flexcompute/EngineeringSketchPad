## [import]
# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse
## [import]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Awave Pytest Example',
                                 prog = 'awave_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ".", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

## [localVariable]
# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AwaveAnalysisTest")
## [localVariable]

# Load CSM file
## [geometry]
geometryScript = os.path.join("..","csmData","awaveWingTailFuselage.csm")

myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

myProblem.geometry.despmtr.area = 10.0
## [geometry]

# Load desired aim
## [loadAIM]
myAnalysis = myProblem.analysis.create( aim = "awaveAIM" )
## [loadAIM]

# Set new Mach and Angle of Attack parameters
## [setInputs]
myAnalysis.input.Mach  = [ 1.2, 1.5]
myAnalysis.input.Alpha = [ 0.0, 2.0]
## [setInputs]

## [output]
CdWave = myAnalysis.output.CDwave
## [output]
print("CdWave = ", CdWave)

# Check assertation
assert abs(CdWave[0]-1.133332) <= 1E-3
assert abs(CdWave[1]-0.216380) <= 1E-4
