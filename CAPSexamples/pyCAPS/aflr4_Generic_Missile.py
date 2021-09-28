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
myProblem = pyCAPS.Problem(problemName = workDir,
                                capsFile=geometryScript, 
                                outLevel=args.outLevel)

# Load AFLR4 aim
myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

# Set AIM verbosity
myAnalysis.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set output grid format since a project name is being supplied - Tecplot  file
myAnalysis.input.Mesh_Format = "Tecplot"

# Farfield growth factor
myAnalysis.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
myAnalysis.input.max_scale = 0.1
myAnalysis.input.min_scale = 0.01

myAnalysis.input.Mesh_Length_Factor = 0.25

# Set project name so a mesh file is generated for each configuration
myAnalysis.input.Proj_Name = "pyCAPS_AFLR4_Missile_Test"

# Dissable curvature refinement when AFLR4 cannot generate a mesh
# myAnalysis.input.Mesh_Gen_Input_String = "auto_mode=0"

# Run AIM
myAnalysis.runAnalysis()

