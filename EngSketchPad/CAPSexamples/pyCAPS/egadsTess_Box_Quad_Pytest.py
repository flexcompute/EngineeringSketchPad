# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Quading Example',
                                 prog = 'egadsTess_Box_Quad_Pytest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "egadsTessBoxQuad")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaBoxes.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)

# Load egadsTess aim
egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM", name = "egads")

# Optional: Explicitly write mesh files
egadsTess.input.Mesh_Format = ["Tecplot", "stl", "ugrid"]

# Global tessellation paramters
egadsTess.input.Tess_Params = [0.1, 0.001, 15]

# Minimum number of points on an edge for quadding
egadsTess.input.Edge_Point_Min = 10

# Maximum number of points on an edge for quadding
egadsTess.input.Edge_Point_Max = 15

# Generate quad meshes
egadsTess.input.Mesh_Elements = "Quad"

# Run AIM
egadsTess.runAnalysis()

#egadsTess.geometry.view()
