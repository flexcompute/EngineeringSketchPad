
# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Interference PyTest Example',
                                 prog = 'interference_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "InterferenceTest")

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","interference.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)

# Load Interference aim
interference = capsProblem.analysis.create(aim = "interferenceAIM")

# Set new EGADS body tessellation parameters
interference.input.Tess_Params = [.1, 0.1, 20.0]

# Run AIM
interference.runAnalysis()

print(interference.output.Names)
print(interference.output.Distances)
print(interference.output.Volumes)
print(interference.output.Areas)
print(interference.output.CGs)
print(interference.output.Inertias)
