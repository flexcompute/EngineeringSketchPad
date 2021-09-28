# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 Symmetry PyTest Example',
                                 prog = 'aflr4_Symmetry_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLR4SymmetryAnalysisTest")

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","cfdSymmetry.csm")
myProblem = pyCAPS.Problem(problemName = workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Load AFLR4 aim
myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

# Mesing boundary conditions
myAnalysis.input.Mesh_Sizing = {"Farfield": {"bcType":"Farfield"},
                                "Symmetry": {"bcType":"Symmetry"}}

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = "pyCAPS_AFLR4_Test"

# Set AIM verbosity
myAnalysis.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set output grid format since a project name is being supplied - Tecplot  file
myAnalysis.input.Mesh_Format = "Tecplot"

# Farfield growth factor
myAnalysis.input.ff_cdfr = 1.4

# Generate quads and tris
# myAnalysis.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"

# Set maximum and minimum edge lengths relative to capsMeshLength
myAnalysis.input.max_scale = 0.2
myAnalysis.input.min_scale = 0.01

# Run AIM
myAnalysis.runAnalysis()