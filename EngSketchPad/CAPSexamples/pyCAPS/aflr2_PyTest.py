# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR2 Pytest Example',
                                 prog = 'aflr2_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLR2AnalysisTest")

# Project name
projectName = "pyCAPS_AFLR2_Test"

# Load CSM file
geometryScript = os.path.join("..","csmData","cfd2D.csm")
myProblem = pyCAPS.Problem(problemName=workDir, 
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Load aflr2 aim
myMesh = myProblem.analysis.create(aim = "aflr2AIM")

# Set project name
myMesh.input.Proj_Name = projectName

# Set AIM verbosity
myMesh.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set airfoil edge parameters
airfoil = {"numEdgePoints" : 100, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.001, 0.001]}

# Set meshing parameters
myMesh.input.Mesh_Sizing = {"Airfoil"   : airfoil,
                            "TunnelWall": {"numEdgePoints" : 50},
                            "InFlow"    : {"numEdgePoints" : 25},
                            "OutFlow"   : {"numEdgePoints" : 25},
                            "2DSlice"   : {"tessParams" : [0.50, .01, 45]}}

# Make quad/tri instead of just quads
myMesh.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"

# Set output grid format since a project name is being supplied - Tecplot  file
myMesh.input.Mesh_Format = "Tecplot"

# Run AIM
myMesh.runAnalysis()
