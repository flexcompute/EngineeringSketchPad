# Import pyCAPS module
import pyCAPS

# Import modules
import os
import argparse
import platform
import time

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 and Pointwise Pytest Example',
                                 prog = 'su2_and_pointwise_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "SU2PointwiseAnalysisTest")

projectName = "pyCAPS_SU2_Pointwise"

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load pointwise aim
pointwise = myProblem.analysis.create(aim = "pointwiseAIM",
                                      name = "pointwise")

# Global Min/Max number of points on edges
pointwise.input.Connector_Initial_Dim = 10
pointwise.input.Connector_Min_Dim     = 3
pointwise.input.Connector_Max_Dim     = 15

# Run AIM pre-analysis
pointwise.preAnalysis()

####### Run pointwise ####################
CAPS_GLYPH = os.environ["CAPS_GLYPH"]
for i in range(30):
    try:
        if platform.system() == "Windows":
            PW_HOME = os.environ["PW_HOME"]
            pointwise.system(PW_HOME + "\\win64\\bin\\tclsh.exe " + CAPS_GLYPH + "\\GeomToMesh.glf caps.egads capsUserDefaults.glf")
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
##########################################

# Run AIM post-analysis
pointwise.postAnalysis()

# Load SU2 aim
su2 = myProblem.analysis.create(aim = "su2AIM",
                                name = "su2")

su2.input["Mesh"].link(myProblem.analysis["pointwise"].output["Volume_Mesh"])

# Set SU2 Version
su2.input.SU2_Version = "Blackbird"

# Set project name
su2.input.Proj_Name = projectName

# Set AoA number
su2.input.Alpha = 1.0

# Set Mach number
su2.input.Mach = 0.5901

# Set equation type
su2.input.Equation_Type = "Compressible"

# Set the output file formats
su2.input.Output_Format = "Tecplot"

# Set number of iterations
su2.input.Num_Iter = 5

# Set boundary conditions
inviscidBC1 = {"bcType" : "Inviscid"}
inviscidBC2 = {"bcType" : "Inviscid"}
su2.input.Boundary_Condition = {"Wing1"   : inviscidBC1,
                                "Wing2"   : inviscidBC2,
                                "Farfield":"farfield"}

# Specifcy the boundares used to compute forces
myProblem.analysis["su2"].input.Surface_Monitor = ["Wing1", "Wing2"]

# Run AIM pre-analysis
su2.preAnalysis()

####### Run SU2 ######################
print ("\n\nRunning SU2......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(su2.analysisDir) # Move into test directory

su2Run(projectName + ".cfg", args.numberProc) # Run SU2

os.chdir(currentDirectory) # Move back to top directory
#######################################

# Run AIM post-analysis
su2.postAnalysis()

# Get force results
print ("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print ("Cl = " , su2.output.CLtot,
       "Cd = " , su2.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx = " , su2.output.CMXtot,
       "Cmy = " , su2.output.CMYtot,
       "Cmz = " , su2.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print ("Cx = " , su2.output.CXtot,
       "Cy = " , su2.output.CYtot,
       "Cz = " , su2.output.CZtot)
