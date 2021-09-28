## [import]
# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run
## [import]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 and Tetgen Pytest Example',
                                 prog = 'su2_and_Tetgen_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "SU2TetgenAnalysisTest")

# -----------------------------------------------------------------
# Load CSM file and Change a design parameter - area in the geometry
# Any despmtr from the cfdMultiBody.csm  file are available inside the pyCAPS script.
# For example: thick, camber, area, aspect, taper, sweep, washout, dihedral
# -----------------------------------------------------------------
# Load CSM file
## [geometry]
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [geometry]

## [capsDespmtrs]
# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.area = 50
## [capsDespmtrs]

## [loadMeshAIM]
# Load EGADS Tess aim
mySurfMesh = myProblem.analysis.create(aim = "egadsTessAIM",
                                       name = "tess")

# Load Tetgen aim
myMesh = myProblem.analysis.create(aim = "tetgenAIM",
                                   name = "myMesh")
## [loadMeshAIM]

## [setMeshpmtrs]
# Set project name so a mesh file is generated
mySurfMesh.input.Proj_Name = "egadsTessMesh"

# Set new EGADS body tessellation parameters
mySurfMesh.input.Tess_Params = [0.5, 0.1, 20.0]

# Set output grid format since a project name is being supplied - Tecplot file
mySurfMesh.input.Mesh_Format = "Tecplot"

# Link surface mesh from EGADS to TetGen
myMesh.input["Surface_Mesh"].link(mySurfMesh.output["Surface_Mesh"])

# Preserve surface mesh while meshing
myMesh.input.Preserve_Surf_Mesh = True
## [setMeshpmtrs]

## [meshAnalysis]
##################################################
## egadsTess/TetGen are executed automatically  ##
##################################################
## [meshAnalysis]

## [loadSu2]
# Load SU2 aim - child of Tetgen AIM
myAnalysis = myProblem.analysis.create(aim = "su2AIM", name = "su2")
## [loadSu2]

## [setSu2Inputs]
# Link the Mesh to the TetGen Volume_Mesh
myAnalysis.input["Mesh"].link(myMesh.output["Volume_Mesh"])

# Set SU2 Version
myAnalysis.input.SU2_Version = "Blackbird"

# Set project name
myAnalysis.input.Proj_Name = "pyCAPS_SU2_Tetgen"

# Set AoA number
myAnalysis.input.Alpha = 1.0

# Set Mach number
myAnalysis.input.Mach = 0.5901

# Set equation type
myAnalysis.input.Equation_Type = "Compressible"

# Set number of iterations
myAnalysis.input.Num_Iter = 5

# Specifcy the boundares used to compute forces
myAnalysis.input.Surface_Monitor = ["Wing1", "Wing2"]
## [setSu2Inputs]

## [setSu2BCs]
# Set boundary conditions
inviscidBC1 = {"bcType" : "Inviscid"}
inviscidBC2 = {"bcType" : "Inviscid"}
myAnalysis.input.Boundary_Condition = {"Wing1"   : inviscidBC1,
                                       "Wing2"   : inviscidBC2,
                                       "Farfield": "farfield"}
## [setSu2BCs]

# Run AIM pre-analysis
## [su2PreAnalysis]
myAnalysis.preAnalysis()
## [su2PreAnalysis]

####### Run SU2 ####################
## [su2Run]
print ("\n\nRunning SU2......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myAnalysis.analysisDir) # Move into test directory

su2Run(myAnalysis.input.Proj_Name + ".cfg", args.numberProc) # Run SU2

os.chdir(currentDirectory) # Move back to top directory
## [su2Run]
#######################################

# Run AIM post-analysis
## [su2PostAnalysis]
myAnalysis.postAnalysis()
## [su2PostAnalysis]

## [su2Output]
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
## [su2Output]