# Import pyCAPS module
## [importModules]
import pyCAPS

import os
import argparse
## [importModules]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D and Tetgen Pytest Example',
                                 prog = 'fun3d_and_Tetgen_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()


# Create working directory variable
## [localVariable]
workDir = os.path.join(str(args.workDir[0]), "FUN3DTetgenAnalysisTest")
## [localVariable]

# Load CSM file
## [loadGeom]
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [loadGeom]

# Change a design parameter - area in the geometry
## [setGeom]
myProblem.geometry.despmtr.area = 50
## [setGeom]

# Write out a new egads file for the geometry
myProblem.geometry.save("pyCAPS_EGADS",
                        directory = workDir,
                        extension="egads") # .egads extension added by default

## [loadMesh]
# Load egadsTess aim
myProblem.analysis.create(aim = "egadsTessAIM", name = "egadsTess")

# Load Tetgen aim
meshAIM = myProblem.analysis.create(aim = "tetgenAIM", name = "tetgen")
## [loadMesh]

## [setMesh]
# Set new EGADS body tessellation parameters
myProblem.analysis["egadsTess"].input.Tess_Params = [1.0, 0.01, 20.0]

# Preserve surface mesh while meshing
meshAIM.input.Preserve_Surf_Mesh = True

# Link Surface_Mesh
meshAIM.input["Surface_Mesh"].link(myProblem.analysis["egadsTess"].output["Surface_Mesh"])
## [setMesh]

# Load FUN3D aim
## [loadFUN3D]
fun3dAIM = myProblem.analysis.create(aim = "fun3dAIM",
                                     name = "fun3d")
## [loadFUN3D]

##[setFUN3D]
# Set project name
fun3dAIM.input.Proj_Name = "fun3dTetgenTest"

# Link the mesh
fun3dAIM.input["Mesh"].link(meshAIM.output["Volume_Mesh"])

fun3dAIM.input.Mesh_ASCII_Flag = False

# Set AoA number
myProblem.analysis["fun3d"].input.Alpha = 1.0

# Set Mach number
myProblem.analysis["fun3d"].input.Mach = 0.5901

# Set equation type
fun3dAIM.input.Equation_Type = "compressible"

# Set Viscous term
myProblem.analysis["fun3d"].input.Viscous = "inviscid"

# Set number of iterations
myProblem.analysis["fun3d"].input.Num_Iter = 10

# Set CFL number schedule
myProblem.analysis["fun3d"].input.CFL_Schedule = [0.5, 3.0]

# Set read restart option
fun3dAIM.input.Restart_Read = "off"

# Set CFL number iteration schedule
myProblem.analysis["fun3d"].input.CFL_Schedule_Iter = [1, 40]

# Set overwrite fun3d.nml if not linking to Python library
myProblem.analysis["fun3d"].input.Overwrite_NML = True
##[setFUN3D]

# Set boundary conditions
##[setFUN3DBC]
inviscidBC1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
inviscidBC2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
fun3dAIM.input.Boundary_Condition = {"Wing1"   : inviscidBC1,
                                     "Wing2"   : inviscidBC2,
                                     "Farfield": "farfield"}
##[setFUN3DBC]

# Run AIM pre-analysis
## [preAnalysisFUN3D]
fun3dAIM.preAnalysis()
## [preAnalysisFUN3D]

####### Run fun3d ####################
print ("\n\nRunning FUN3D......")
if (args.noAnalysis == False):
    fun3dAIM.system("nodet_mpi --animation_freq -1 --write_aero_loads_to_file> Info.out"); # Run fun3d via system call
#######################################

# Run AIM post-analysis
## [postAnalysiFUN3D]
fun3dAIM.postAnalysis()
## [postAnalysiFUN3D]

# Get force results
##[results]
print ("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print ("Cl = " , fun3dAIM.output.CLtot,
       "Cd = " , fun3dAIM.output.CDtot)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx = " , fun3dAIM.output.CMXtot,
       "Cmy = " , fun3dAIM.output.CMYtot,
       "Cmz = " , fun3dAIM.output.CMZtot)

# Get Cx, Cy, Cz coefficients
print ("Cx = " , fun3dAIM.output.CXtot,
       "Cy = " , fun3dAIM.output.CYtot,
       "Cz = " , fun3dAIM.output.CZtot)

print ("Pressure Contribution")
# Get Lift and Drag coefficients
print ("Cl_p = " , fun3dAIM.output.CLtot_p,
       "Cd_p = " , fun3dAIM.output.CDtot_p)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_p = " , fun3dAIM.output.CMXtot_p,
       "Cmy_p = " , fun3dAIM.output.CMYtot_p,
       "Cmz_p = " , fun3dAIM.output.CMZtot_p)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_p = " , fun3dAIM.output.CXtot_p,
       "Cy_p = " , fun3dAIM.output.CYtot_p,
       "Cz_p = " , fun3dAIM.output.CZtot_p)

print ("Viscous Contribution")
# Get Lift and Drag coefficients
print ("Cl_v = " , fun3dAIM.output.CLtot_v,
       "Cd_v = " , fun3dAIM.output.CDtot_v)

# Get Cmx, Cmy, and Cmz coefficients
print ("Cmx_v = " , fun3dAIM.output.CMXtot_v,
       "Cmy_v = " , fun3dAIM.output.CMYtot_v,
       "Cmz_v = " , fun3dAIM.output.CMZtot_v)

# Get Cx, Cy, and Cz, coefficients
print ("Cx_v = " , fun3dAIM.output.CXtot_v,
       "Cy_v = " , fun3dAIM.output.CYtot_v,
       "Cz_v = " , fun3dAIM.output.CZtot_v)
##[results]

# Print out dictionary of component forces
forces = fun3dAIM.output.Forces
for i in forces.keys():
    print ("Boundary name: ", i)
    for j in forces[i].keys():
        print (" ", j," = ", forces[i][j])
