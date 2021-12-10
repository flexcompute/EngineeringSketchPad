# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Tetgen Pytest Example',
                                 prog = 'tetgen_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "TetgenAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load EGADS Tess aim
myProblem.analysis.create(aim = "egadsTessAIM", name = "egadsTess")

# Set project name so a mesh file is generated
myProblem.analysis["egadsTess"].input.Proj_Name = "egadsTessMesh"

# Set new EGADS body tessellation parameters
myProblem.analysis["egadsTess"].input.Tess_Params = [0.5, 0.1, 20.0]

# Set output grid format since a project name is being supplied - Tecplot file
myProblem.analysis["egadsTess"].input.Mesh_Format = "Tecplot"


##########################################
## egadsTess runs automatically         ##
##########################################


# Load TetGen aim
myProblem.analysis.create(aim = "tetgenAIM", name = "tetgen")

myProblem.analysis["tetgen"].input["Surface_Mesh"].link(myProblem.analysis["egadsTess"].output["Surface_Mesh"])

# Set Tetgen: maximum radius-edge ratio
myProblem.analysis["tetgen"].input.Quality_Rad_Edge = 1.5

# Set surface mesh preservation
myProblem.analysis["tetgen"].input.Preserve_Surf_Mesh = True

# Set project name
myProblem.analysis["tetgen"].input.Proj_Name = "pyCAPS_Tetgen_Test"

# Set output grid format since a project name is being supplied - Tecplot tetrahedral file
myProblem.analysis["tetgen"].input.Mesh_Format = "Tecplot"
myProblem.analysis["tetgen"].input.Mesh_ASCII_Flag = False


# Run AIM
myProblem.analysis["tetgen"].runAnalysis()
