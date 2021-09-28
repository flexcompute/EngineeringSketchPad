# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 and AFLR Pytest Example',
                                 prog = 'su2_and_AFLR4_AFLR3_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "su2_and_AFLR4_AFLR3")

# CSM file
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.area = 50

# Load AFLR4 aim
aflr4 = myProblem.analysis.create(aim = "aflr4AIM",
                                  name = "aflr4")

# Set project name so a mesh file is generated
aflr4.input.Proj_Name = "AFLR4_Mesh"

# Set mesh ascii flag
aflr4.input.Mesh_ASCII_Flag = True

# Set output mesh format
aflr4.input.Mesh_Format = "Tecplot"

# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Farfield growth factor
aflr4.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.5
aflr4.input.min_scale = 0.1

# Load AFLR3 aim - child of AFLR4 AIM
aflr3 = myProblem.analysis.create(aim = "aflr3AIM",
                                  name = "aflr3")

# Link the aflr4 output Surface_Mesh to the aflr3 input Surface_Mesh
aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

# Set project name so a mesh file is generated
aflr3.input.Proj_Name = "AFLR3_Mesh"

# Set mesh ascii flag
aflr3.input.Mesh_ASCII_Flag = True

# Set project name - so a mesh file is generated
aflr3.input.Mesh_Format = "Tecplot"

# Set AIM verbosity
aflr3.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False


# Load SU2 aim - child of AFLR3 AIM
su2 = myProblem.analysis.create(aim = "su2AIM",
                                name = "su2")

# Link the aflr3 output Volume_Mesh to the su2 input Mesh
su2.input["Mesh"].link(aflr3.output["Volume_Mesh"])

projectName = "su2AFLRTest"

# Set SU2 Version
su2.input.SU2_Version = "Blackbird"

# Set project name
su2.input.Proj_Name = projectName

# Set AoA number
su2.input.Alpha = 1.0

# Set Mach number
su2.input.Mach = 0.5901

# Set Reynolds number
#su2.input.Re = 10e5

# Set equation type
su2.input.Equation_Type = "compressible"

# Set equation type
su2.input.Physical_Problem = "Euler"

# Set number of iterations
su2.input.Num_Iter = 5

# Set number of inner iterations, testing Input_String
su2.input.Input_String = ["INNER_ITER= 10",
                          "CONV_FIELD= RMS_DENSITY",
                          "CONV_RESIDUAL_MINVAL= -8"]

# Set overwrite su2 cfg if not linking to Python library
su2.input.Overwrite_CFG = True

# Set boundary conditions
bc1 = {"bcType" : "Viscous", "wallHeatFlux" : 0}
bc2 = {"bcType" : "Inviscid"}
su2.input.Boundary_Condition = {"Wing1": bc1,
                                "Wing2": bc2,
                                "Farfield":"farfield"}

# Specifcy the boundares used to compute forces
su2.input.Surface_Monitor = ["Wing1", "Wing2"]


#################################################
##  AFLR4/AFLR3 are executed automatically      #
#################################################

# Run AIM pre-analysis
su2.preAnalysis()

####### Run su2 ####################
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

print ("Pressure Contribution")
# Get Lift and Drag coefficients
print ("Cl_p = " , su2.output.CLtot_p,
       "Cd_p = " , su2.output.CDtot_p)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_p = " , su2.output.CMXtot_p,
       "Cmy_p = " , su2.output.CMYtot_p,
       "Cmz_p = " , su2.output.CMZtot_p)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_p = " , su2.output.CXtot_p,
       "Cy_p = " , su2.output.CYtot_p,
       "Cz_p = " , su2.output.CZtot_p)

print ("Viscous Contribution")
# Get Lift and Drag coefficients
print ("Cl_v = " , su2.output.CLtot_v,
       "Cd_v = " , su2.output.CDtot_v)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_v = " , su2.output.CMXtot_v,
       "Cmy_v = " , su2.output.CMYtot_v,
       "Cmz_v = " , su2.output.CMZtot_v)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_v = " , su2.output.CXtot_v,
       "Cy_v = " , su2.output.CYtot_v,
       "Cz_v = " , su2.output.CZtot_v)
