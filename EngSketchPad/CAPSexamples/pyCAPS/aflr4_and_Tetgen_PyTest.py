# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 and Tetgen Pytest Example',
                                 prog = 'aflr4_and_Tetgen_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Set analysis directory
workDir = os.path.join(str(args.workDir[0]), "AFLR4TetgenAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName = workDir,
                           capsFile = geometryScript, 
                           outLevel = args.outLevel)

# Load AFLR4 aim
mySurfMesh = myProblem.analysis.create(aim = "aflr4AIM")

# Set project name - so a mesh file is generated
mySurfMesh.input.Proj_Name = "pyCAPS_AFLR4_Tetgen"

# Set AIM verbosity
mySurfMesh.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set output grid format since a project name is being supplied - Tecplot  file
mySurfMesh.input.Mesh_Format = "Tecplot"

# Farfield growth factor
mySurfMesh.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
mySurfMesh.input.max_scale = 0.2
mySurfMesh.input.min_scale = 0.01

######################################
## AFRL4 executes automatically     ##
######################################

#######################################
## Build volume mesh off of surface  ##
##  mesh(es) using TetGen            ##
#######################################

# Load TetGen aim
myVolMesh = myProblem.analysis.create(aim = "tetgenAIM")

# Link the surface mesh
myVolMesh.input["Surface_Mesh"].link(mySurfMesh.output["Surface_Mesh"])

# Set AIM verbosity
myVolMesh.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set surface mesh preservation
myVolMesh.input.Preserve_Surf_Mesh = True

# Set project name - so a mesh file is generated
myVolMesh.input.Proj_Name = "pyCAPS_AFLR4_Tetgen_VolMesh"

# Set output grid format since a project name is being supplied - Tecplot tetrahedral file
myVolMesh.input.Mesh_Format = "Tecplot"

# Don't extract the update the local-to-global node connectivity
myVolMesh.input.Ignore_Surface_Mesh_Extraction = True

# Run AIM
myVolMesh.runAnalysis()
