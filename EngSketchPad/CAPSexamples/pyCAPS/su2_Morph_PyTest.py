# Import pyCAPS module
import pyCAPS

# Import os module
import os
import sys

# Import argparse module
import argparse

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 Mesh Morph w/ AFLR Pytest Example',
                                 prog = 'su2_Morph_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

def mesh_Setup(myProblem, writeMesh = False):

    # Load AFLR4 aim
    aflr4 = myProblem.analysis.create(aim = "aflr4AIM", name = "surfaceMesh")
    
    # Optional: Explicitly write mesh files
    if writeMesh:
        aflr4.input.Mesh_Format = "Tecplot"
    
    # Set AIM verbosity
    aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False
    
    # Farfield growth factor
    aflr4.input.ff_cdfr = 1.4
    
    # Set maximum and minimum edge lengths relative to capsMeshLength
    aflr4.input.max_scale = 0.5
    aflr4.input.min_scale = 0.05
    
    aflr4.input.Mesh_Length_Factor = 1

    #######################################
    ## Build volume mesh off of surface  ##
    ##  mesh(es) using AFLR3             ##
    #######################################
    
    # Load AFLR3 aim
    aflr3 = myProblem.analysis.create(aim = "aflr3AIM", name = "volumeMesh")
    
    aflr3.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])
    
    # Set AIM verbosity
    aflr3.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False
    
    # Optional: Explicitly write mesh files
    if writeMesh:
        aflr3.input.Mesh_Format = "Tecplot"
    
    return aflr4, aflr3

def flow_Setup(myProblem):
    
    su2 = myProblem.analysis.create(aim = "su2AIM", name = "su2")
    
    su2.input["Mesh"].link(myProblem.analysis["volumeMesh"].output["Volume_Mesh"])
    
    # Set project name
    su2.input.Proj_Name = "su2AFLRTest"

    # Set SU2 Version
    su2.input.SU2_Version = "Harrier"
    
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
    bc1 = {"bcType" : "Inviscid"}
    bc2 = {"bcType" : "Inviscid"}
    su2.input.Boundary_Condition = {"Wing1": bc1,
                                    "Wing2": bc2,
                                    "Farfield":"farfield",
                                    "Wake":{"bcType":"Internal"}}
    
    # Specifcy the boundares used to compute forces
    su2.input.Surface_Monitor = ["Wing1", "Wing2"]
    
    
    su2.input.Output_Format = "Paraview"
    
    # Set flag to indicate mesh morphing is desired
    su2.input.Mesh_Morph = True
    
    return su2


def flow_Run(su2, run=True, morph=False, redirect=False):
    
    su2.preAnalysis()

    currentDirectory = os.getcwd() # Get our current working directory
    
    os.chdir(su2.analysisDir) # Move into test directory

    if morph: 
        print ("\n\nMesh Morphing......") 
        # Run SU2 mesh deformation
        # Work around python 3 bug in su2Deform function
        if sys.version_info[0] < 3:
            su2Deform(su2.input.Proj_Name + ".cfg", args.numberProc)
        else:
            os.system("SU2_DEF " + su2.input.Proj_Name + ".cfg")

        
    ####### Run su2 ####################
    print ("\n\nRunning SU2......") 
    
    su2Run(su2.input.Proj_Name + ".cfg", args.numberProc) # Run SU2
    
    os.chdir(currentDirectory) # Move back to top directory
    
    su2.postAnalysis()
    
    # Get force results
    print ("Total Force - Pressure + Viscous")
    # Get Lift and Drag coefficients
    print ("Cl = " , su2.output.CLtot,
           "Cd = " , su2.output.CDtot)

    return (su2.output.CLtot,
            su2.output.CDtot)

#################################

# Working directory
workDir = os.path.join(str(args.workDir[0]), "SU2_Morphing")


# Load CSM file and build the geometry explicitly
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=os.path.join("..","csmData","cfdMultiBody.csm"),
                           outLevel=args.outLevel)

# Change geometry 
myProblem.geometry.despmtr.twist = 0

# Setup
mesh_Setup(myProblem, writeMesh = True)
su2 = flow_Setup(myProblem)
 
# Run once normally
cl1, cd1 = flow_Run(su2)

# Unlink the mesh to force a morph
su2.input["Mesh"].unlink() # Force a morph

# Change geometry 
myProblem.geometry.despmtr.twist = 16

# Re-run with deformed mesh
cl2, cd2 = flow_Run(su2, morph=True)

assert(cl2 > cl1)
assert(cd2 > cd1)

