#------------------------------------------------------------------------------#
# Import pyCAPS module
import pyCAPS

# Import other modules
import os
import f90nml # f90nml is used to write fun3d inputs not available in the aim

# Import argparse module
import argparse
#------------------------------------------------------------------------------#


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D Mesh Morph w/ AFLR Pytest Example',
                                 prog = 'fun3d_Morph_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument('-cores', default = 1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

#------------------------------------------------------------------------------#

# Working directory
workDir = os.path.join(str(args.workDir[0]), "FUN3D_refine")

# Load CAPS Problem
capsProblem = pyCAPS.Problem(problemName = workDir,
                             capsFile=os.path.join("..","csmData","cfdSymmetry.csm"),
                             outLevel = args.outLevel);

#------------------------------------------------------------------------------#

# Create aflr4 AIM
aflr4 = capsProblem.analysis.create(aim  = "aflr4AIM",
                                    name = "aflr4")

# Make the mesh coarse
#aflr4.input.Mesh_Length_Factor = 10

# Farfield growth factor
aflr4.input.ff_cdfr = 1.8

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.2
aflr4.input.min_scale = 0.01

# Mesing boundary conditions
aflr4.input.Mesh_Sizing = {"Farfield": {"bcType":"Farfield"},
                           "Symmetry": {"bcType":"Symmetry"}}

#------------------------------------------------------------------------------#

# Create AFLR3 AIM to generate the volume mesh
aflr3 = capsProblem.analysis.create(aim  = "aflr3AIM",
                                    name = "alfr3")

# Link the aflr4 Surface_Mesh as input to aflr3
aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

#------------------------------------------------------------------------------#

# Create refine AIM
refine = capsProblem.analysis.create(aim = "refineAIM",
                                     name= "refine")

# Limit the number of passes to speed up this example (generally this should be the default)
refine.input.Passes = 2

# Define executeable string
refine.input.ref="mpiexec -n 4 refmpifull";

# Link the aflr3 mesh as the initial mesh for adaptation
refine.input["Mesh"].link(aflr3.output["Volume_Mesh"]);

#------------------------------------------------------------------------------#

# Create Fun3D AIM
fun3d = capsProblem.analysis.create(aim  = "fun3dAIM",
                                    name = "fun3d")

# Link the aflr3 Volume_Mesh as input to fun3d
fun3d.input["Mesh"].link(refine.output["Mesh"])

# Set project name. Files written to fun3d.analysisDir will have this name
fun3d.input.Proj_Name = "refineTest"

speedofSound = 340.0 # m/s
refVelocity  = 0.7*speedofSound # m/s
refDensity   = 1.2   # kg/m^3

fun3d.input.Alpha = 10.0                     # AoA
fun3d.input.Mach  = refVelocity/speedofSound # Mach number
fun3d.input.Equation_Type = "compressible"   # Equation type
fun3d.input.Num_Iter = 5                     # Number of iterations
fun3d.input.Restart_Read = 'off'             # Do not read restart
fun3d.input.Viscous = "inviscid"             # Inviscid calculation

# Set scaling factor for pressure data transfer
fun3d.input.Pressure_Scale_Factor = 0.5*refDensity*refVelocity**2
fun3d.input.Pressure_Scale_Offset = 0.0

# Set boundary conditions via capsGroup
inviscidBC = {"bcType" : "Inviscid"}

bc1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
fun3d.input.Boundary_Condition = {"Wing1"    : bc1,
                                  "Symmetry" : "SymmetryY",
                                  "Farfield" : "farfield"}

# Use python to add inputs to fun3d.nml file
fun3d.input.Use_Python_NML = True

# Write boundary output variables to the fun3d.nml file directly
fun3dnml = f90nml.Namelist()
fun3dnml['boundary_output_variables'] = f90nml.Namelist()
fun3dnml['boundary_output_variables']['mach'] = True
fun3dnml['boundary_output_variables']['cp'] = True
fun3dnml['boundary_output_variables']['average_velocity'] = True

# code runcontrol for FUN3D
fun3dnml['code_run_control'] = f90nml.Namelist()
fun3dnml['code_run_control']['stopping_tolerance'] = 5e-8

# HANIM Solver FUN3D
fun3dnml['nonlinear_solver_parameters'] = f90nml.Namelist()
fun3dnml['nonlinear_solver_parameters']['hanim'] = True

# additional FUN3D inputs to generate solb file with primitive variables for refine
fun3dnml['volume_output_variables'] = f90nml.Namelist()
fun3dnml['volume_output_variables']['mach'] = False
fun3dnml['volume_output_variables']['x']    = False
fun3dnml['volume_output_variables']['y']    = False
fun3dnml['volume_output_variables']['z']    = False
fun3dnml['volume_output_variables']['primitive_variables'] = True
fun3dnml['volume_output_variables']['export_to'] = 'solb'

fun3dnml['global'] = f90nml.Namelist()
fun3dnml['global']['volume_animation_freq'] = -1

fun3dnml.write(os.path.join(fun3d.analysisDir,"fun3d.nml"), force=True)

# Sepcify the file name for the primitive variable field use to compute multiscale metric
refine.input.ScalarFieldFile = os.path.join(fun3d.analysisDir, fun3d.input.Proj_Name + "_volume.solb");

# Indicate that solb file contains Fun3D primitive variables
# and generate a 'restart' file for Fun3d
refine.input.Fun3D = True

#------------------------------------------------------------------------------#

# Specficy number of adaptation iterations for each target complexity
nAdapt     = [   2,     2]
tar_Complex= [8000, 16000]
for ii,iComplex in enumerate(tar_Complex):
    for iadapt in range(0,nAdapt[ii]):
        # Set the current mesh size
        refine.input["Complexity"].value = iComplex

        if iadapt > 0 or iComplex > tar_Complex[0]:
            fun3dnml["flow_initialization"]=f90nml.Namelist()
            fun3dnml["flow_initialization"]['import_from'] = os.path.join(refine.analysisDir, 'refine_out-restart.solb')
            fun3dnml.write(os.path.join(fun3d.analysisDir,"fun3d.nml"), force=True)

        ####### Run abaqus ####################
        fun3d.preAnalysis()
        print("\n==> Running FUN3D, complexity = %6d, iadapt = %2d......" % (iComplex,iadapt))
        fun3d.system("mpiexec -np 4 nodet_mpi --write_aero_loads_to_file > Info.txt")
        fun3d.postAnalysis()
        #######################################

        # Get force results
        print ("\n==> Total Forces and Moments")
        # Get Lift and Drag coefficients
        print ("--> Cl = ", fun3d.output.CLtot,
                   "Cd = ", fun3d.output.CDtot)

        # Get Cmx, Cmy, and Cmz coefficients
        print ("--> Cmx = ", fun3d.output.CMXtot,
                   "Cmy = ", fun3d.output.CMYtot,
                   "Cmz = ", fun3d.output.CMZtot)

        # Get Cx, Cy, Cz coefficients
        print ("--> Cx = ", fun3d.output.CXtot,
                   "Cy = ", fun3d.output.CYtot,
                   "Cz = ", fun3d.output.CZtot)

        # Unlink the refine mesh input
        if iadapt == 0: refine.input["Mesh"].unlink()

capsProblem.closePhase()
