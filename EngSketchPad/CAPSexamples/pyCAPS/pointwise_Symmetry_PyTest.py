
# Import pyCAPS class file
import pyCAPS

# Import time to sleep in-between attempts to call pointwise
import time

# Import argparse module
import os
import argparse
import platform

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Pointwise Pytest Example',
                                 prog = 'pointwise_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "PointwiseSymmetryAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdSymmetry.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AFLR4 aim
pointwise = myProblem.analysis.create(aim = "pointwiseAIM")

# Set project name so a mesh file is generated
pointwise.input.Proj_Name = "pyCAPS_Pointwise_Test"
pointwise.input.Mesh_Format = "Tecplot"

# Block level for viscous mesh
pointwise.input.Block_Full_Layers = 1
pointwise.input.Block_Max_Layers  = 100

# Set mesh sizing parmeters, only Wing2 is viscous
viscousBC  = {"boundaryLayerSpacing" : 0.001}
pointwise.input.Mesh_Sizing = {"Wing2": viscousBC}

# Run AIM pre-analysis
pointwise.preAnalysis()

# Pointwise batch execution notes (PW_HOME referes to the Pointwise installation directory):
# Windows (DOS prompt) -
#     PW_HOME\win64\bin\tclsh %CAPS_GLYPH%\GeomToMesh.glf caps.egads capsUserDefaults.glf
# Linux/macOS (terminal) -
#     PW_HOME/pointwise -b $CAPS_GLYPH/GeomToMesh.glf caps.egads capsUserDefaults.glf

# Attempt to run pointwise repeatedly in case a licensce is not readily available
CAPS_GLYPH = os.environ["CAPS_GLYPH"]
for i in range(30):
    try:
        if platform.system() == "Windows":
            PW_HOME = os.environ["PW_HOME"]
            pointwise.system(PW_HOME + "\win64\bin\tclsh.exe " + CAPS_GLYPH + "\\GeomToMesh.glf caps.egads capsUserDefaults.glf")
        elif "CYGWIN" in platform.system():
            PW_HOME = os.environ["PW_HOME"]
            pointwise.system(PW_HOME + "/win64/bin/tclsh.exe " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
        else:
            pointwise.system("pointwise -b " + CAPS_GLYPH + "/GeomToMesh.glf caps.egads capsUserDefaults.glf")
    except pyCAPS.CAPSError:
        time.sleep(10) # wait and try again
        continue

    time.sleep(1) # let the harddrive breathe
    if os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.gma')) and \
      (os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.ugrid')) or \
       os.path.isfile(os.path.join(pointwise.analysisDir,'caps.GeomToMesh.lb8.ugrid'))): break
    time.sleep(10) # wait and try again

# Run AIM post-analysis
pointwise.postAnalysis()

#pointwise.geometry.view()
