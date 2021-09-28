# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Delaundo Pytest Example',
                                 prog = 'delaundo_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = str(args.workDir[0]) + "/DelaundoAnalysisTest"

# Load CSM file
geometryScript = os.path.join("..","csmData","cfd2D.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Set the sharp trailing edge geometry design parameter
myProblem.geometry.cfgpmtr.sharpTE = 0

# Load delaundo aim
myMesh = myProblem.analysis.create(aim = "delaundoAIM")

# Set airfoil edge parameters
airfoil = {"numEdgePoints" : 100, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.001, 0.001]}

# Set meshing parameters
myMesh.input.Mesh_Sizing = {"Airfoil"   : airfoil,
                            "AirfoilTE" : {"numEdgePoints" : 4},
                            "TunnelWall": {"numEdgePoints" : 5},
                            "InFlow"    : {"numEdgePoints" : 5},
                            "OutFlow"   : {"numEdgePoints" : 5}}

# Thickness of stretched region
myMesh.input.Delta_Thickness = .1

# Maximum aspect ratio
myMesh.input.Max_Aspect = 90.0

# Set project name and output mesh type
projectName = "delaundoMesh"
myMesh.input.Proj_Name = projectName

myMesh.input.Mesh_Format = "Tecplot"

# Run AIM
myMesh.runAnalysis()
