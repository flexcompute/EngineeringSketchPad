# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D and AFLR Pytest Example',
                                 prog = 'fun3d_and_AFLR4_AFLR3_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "FUN3DandAFRLMeshTest")


# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AFLR4 aim
aflr4 = myProblem.analysis.create(aim = "aflr4AIM", name = "aflr4")

# Set project name - so a mesh file is generated
aflr4.input.Proj_Name = "pyCAPS_AFLR4_AFLR3"

# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set output grid format since a project name is being supplied - Tecplot  file
aflr4.input.Mesh_Format = "Tecplot"

# Farfield growth factor
aflr4.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.5
aflr4.input.min_scale = 0.05

#######################################
## Build volume mesh off of surface  ##
##  mesh(es) using AFLR3             ##
#######################################

# Load AFLR3 aim
aflr3 = myProblem.analysis.create(aim = "aflr3AIM", name = "aflr3")

aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

# Set AIM verbosity
aflr3.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set project name - so a mesh file is generated
aflr3.input.Proj_Name = "pyCAPS_AFLR4_AFLR3_VolMesh"

# Set output grid format since a project name is being supplied - Tecplot tetrahedral file
aflr3.input.Mesh_Format = "Tecplot"

#######################################
## Load FUN3D aim                    ##
#######################################

fun3d = myProblem.analysis.create(aim = "fun3dAIM", name = "fun3d")

fun3d.input["Mesh"].link(aflr3.output["Volume_Mesh"])

# Set project name
fun3d.input.Proj_Name = "fun3dAFLRTest"

# Reset the mesh ascii flag to false - This will generate a lb8.ugrid file
fun3d.input.Mesh_ASCII_Flag = False

# Set AoA number
fun3d.input.Alpha = 1.0

# Set Viscous term
fun3d.input.Viscous = "laminar"

# Set Mach number
fun3d.input.Mach = 0.5901

# Set Reynolds number
fun3d.input.Re = 10E5

# Set equation type
fun3d.input.Equation_Type = "compressible"

# Set number of iterations
fun3d.input.Num_Iter = 10

# Set CFL number schedule
fun3d.input.CFL_Schedule = [0.5, 3.0]

# Set read restart option
fun3d.input.Restart_Read = "off"

# Set CFL number iteration schedule
fun3d.input.CFL_Schedule_Iter = [1, 40]

# Set overwrite fun3d.nml if not linking to Python library
fun3d.input.Overwrite_NML = True

# Set boundary conditions
bc1 = {"bcType" : "Viscous", "wallTemperature" : 1}
bc2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
fun3d.input.Boundary_Condition = {"Wing1"   : bc1,
                                  "Wing2"   : bc2,
                                  "Farfield": "farfield"}

# Run AIM pre-analysis
fun3d.preAnalysis()

####### Run fun3d ####################
print ("\n\nRunning FUN3D......")
# Mesh is to large for a single core
if (args.noAnalysis == False):
    fun3d.system("mpirun -np " + str(args.numberProc) + " nodet_mpi  --animation_freq -1 --volume_animation_freq -1 > Info.out"); # Run fun3d via system call
#######################################

# Run AIM post-analysis
fun3d.postAnalysis()

# Get force results
print ("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print ("Cl = " , fun3d.output.CLtot,
       "Cd = " , fun3d.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx = " , fun3d.output.CMXtot,
       "Cmy = " , fun3d.output.CMYtot,
       "Cmz = " , fun3d.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print ("Cx = " , fun3d.output.CXtot,
       "Cy = " , fun3d.output.CYtot,
       "Cz = " , fun3d.output.CZtot)

print ("Pressure Contribution")
# Get Lift and Drag coefficients
print ("Cl_p = " , fun3d.output.CLtot_p,
       "Cd_p = " , fun3d.output.CDtot_p)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_p = " , fun3d.output.CMXtot_p,
       "Cmy_p = " , fun3d.output.CMYtot_p,
       "Cmz_p = " , fun3d.output.CMZtot_p)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_p = " , fun3d.output.CXtot_p,
       "Cy_p = " , fun3d.output.CYtot_p,
       "Cz_p = " , fun3d.output.CZtot_p)

print ("Viscous Contribution")
# Get Lift and Drag coefficients
print ("Cl_v = " , fun3d.output.CLtot_v,
       "Cd_v = " , fun3d.output.CDtot_v)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_v = " , fun3d.output.CMXtot_v,
       "Cmy_v = " , fun3d.output.CMYtot_v,
       "Cmz_v = " , fun3d.output.CMZtot_v)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_v = " , fun3d.output.CXtot_v,
       "Cy_v = " , fun3d.output.CYtot_v,
       "Cz_v = " , fun3d.output.CZtot_v)
