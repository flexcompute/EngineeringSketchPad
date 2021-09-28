
# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'EGADS Tess Shheres Quad PyTest Example',
                                 prog = 'egadsTess_Nose_Quad_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "egadsTessNoseAnalysisTest")

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","nosecone1.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Load AFLR4 aim
myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = "egadsTessMesh"

# Set new EGADS body tessellation parameters
myAnalysis.input.Tess_Params = [.1, 0.05, 20.0]

# Set output grid format since a project name is being supplied - Tecplot file
myAnalysis.input.Mesh_Format = "Tecplot"

myAnalysis.input.Mesh_Elements = "Quad"

#Mesh_Sizing = {"nose":  {"tessParams" : [0.5, .1, 30]}}

#myAnalysis.input.Mesh_Sizing = Mesh_Sizing

# Run AIM
myAnalysis.runAnalysis()

#myAnalysis.geometry.view()
