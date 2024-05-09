# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 Generic Missile PyTest Example',
                                 prog = 'aflr4_Generic_Missile_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot surface meshes")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLR4GenericMissileAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","generic_missile.csm")
capsProblem = pyCAPS.Problem(problemName = workDir,
                             capsFile=geometryScript, 
                             outLevel=args.outLevel)

# Load AFLR4 aim
aflr4 = capsProblem.analysis.create(aim = "aflr4AIM")

# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Optional: Explicitly write mesh files
aflr4.input.Mesh_Format = ["Tecplot", "ugrid"]

# Farfield growth factor
aflr4.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.1
aflr4.input.min_scale = 0.01

aflr4.input.Mesh_Length_Factor = 0.25

# Dissable curvature refinement when AFLR4 cannot generate a mesh
# aflr4.input.Mesh_Gen_Input_String = "auto_mode=0"

# Run AIM
aflr4.runAnalysis()

