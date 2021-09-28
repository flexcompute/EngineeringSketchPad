
# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'EGADS Tess PyTest Example',
                                 prog = 'egadsTess_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "egadsTessAnalysisTest")

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","feaWingBEMAero.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Load EGADS Tess aim
myAnalysis = myProblem.analysis.create(aim = "egadsTessAIM")

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = "egadsTessMesh"

# Set new EGADS body tessellation parameters
myAnalysis.input.Tess_Params = [.1, 0.1, 20.0]

# Set output grid format since a project name is being supplied - Tecplot file
myAnalysis.input.Mesh_Format = "Tecplot"

myAnalysis.input.Mesh_Elements = "Quad"

Mesh_Sizing = {"LeadingEdge":  {"tessParams" : [0.5, .1, 30]}}
#Mesh_Sizing = {"LeadingEdge":  {"numEdgePoints" : 3}}

myAnalysis.input.Mesh_Sizing = Mesh_Sizing

# Run AIM
myAnalysis.runAnalysis()

#myAnalysis.geometry.view()
