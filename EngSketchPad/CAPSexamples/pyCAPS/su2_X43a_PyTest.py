# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 X43a Pytest Example',
                                 prog = 'su2_X43a_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "SU2X43aAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdX43a.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.tailLength = 150


# Load EGADS Tess aim
mySurfMesh = myProblem.analysis.create(aim = "egadsTessAIM",
                                       name = "tess")

# Set project name so a mesh file is generated
mySurfMesh.input.Proj_Name = "egadsTessMesh"

# Set new EGADS body tessellation parameters
mySurfMesh.input.Tess_Params = [0.5, 0.1, 20.0]

# Set output grid format since a project name is being supplied - Tecplot file
mySurfMesh.input.Mesh_Format = "Tecplot"

##########################################
## egadsTess is executed automatically  ##
##########################################


# Load Tetgen aim
myMesh = myProblem.analysis.create(aim = "tetgenAIM",
                                   name = "myMesh")

myMesh.input["Surface_Mesh"].link(mySurfMesh.output["Surface_Mesh"])

# Preserve surface mesh while meshing
myMesh.input.Preserve_Surf_Mesh = True

##########################################
## TetGen is executed automatically     ##
##########################################


# Load SU2 aim - child of Tetgen AIM
myAnalysis = myProblem.analysis.create(aim = "su2AIM", name = "su2")

myAnalysis.input["Mesh"].link(myMesh.output["Volume_Mesh"])

projectName = "x43a_Test"

# Set SU2 Version
myAnalysis.input.SU2_Version = "Blackbird"

# Set project name
myAnalysis.input.Proj_Name = projectName

# Set AoA number
myAnalysis.input.Alpha = 3.0

# Set Mach number
myAnalysis.input.Mach = 3.5

# Set equation type
myAnalysis.input.Equation_Type = "Compressible"

# Set number of iterations
myAnalysis.input.Num_Iter = 5

# Set overwrite cfg file
myProblem.analysis["su2"].input.Overwrite_CFG = True

# Set boundary conditions
inviscidBC = {"bcType" : "Inviscid"}
myProblem.analysis["su2"].input.Boundary_Condition = {"x43A": inviscidBC,
                                                      "Farfield":"farfield"}

# Specifcy the boundares used to compute forces
myProblem.analysis["su2"].input.Surface_Monitor = ["x43A"]

# Run AIM pre-analysis
myAnalysis.preAnalysis()

####### Run SU2 ####################
print ("\n\nRunning SU2......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myAnalysis.analysisDir) # Move into test directory

su2Run(myAnalysis.input.Proj_Name + ".cfg", args.numberProc) # Run SU2

os.chdir(currentDirectory) # Move back to top directory
#######################################

# Run AIM post-analysis
myAnalysis.postAnalysis()

# Get force results
print("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print("Cl = " , myAnalysis.output.CLtot,
      "Cd = " , myAnalysis.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print("Cmx = " , myAnalysis.output.CMXtot,
      "Cmy = " , myAnalysis.output.CMYtot,
      "Cmz = " , myAnalysis.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print("Cx = " , myAnalysis.output.CXtot,
      "Cy = " , myAnalysis.output.CYtot,
      "Cz = " , myAnalysis.output.CZtot)