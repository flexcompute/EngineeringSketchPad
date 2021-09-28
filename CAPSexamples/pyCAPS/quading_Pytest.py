# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Quading Example',
                                 prog = 'quading_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AstrosQuadding")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaBoxes.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
astros = myProblem.analysis.create(aim = "astrosAIM", name = "astros", autoExec = False)

# Set project name so a mesh file is generated
projectName = "quadding_test"
astros.input.Proj_Name = projectName

# Global tessellation paramters
astros.input.Tess_Params = [0.1, 0.001, 15]

# Minimum number of points on an edge for quadding
astros.input.Edge_Point_Min = 10

# Maximum number of points on an edge for quadding
astros.input.Edge_Point_Max = 15

# Generate quad meshes
astros.input.Quad_Mesh = True


# Set analysis
eigen = { "extractionMethod"     : "MGIV",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormaliztion"    : "MASS"}
astros.input.Analysis = {"EigenAnalysis": eigen}


# Run AIM pre-analysis to genearte bdf file that demonstrates quadding
astros.preAnalysis()
