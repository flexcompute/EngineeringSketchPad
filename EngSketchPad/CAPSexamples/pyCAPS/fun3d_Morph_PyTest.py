# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

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

def mesh_Setup(myProblem, writeMesh = False):

    # Load AFLR4 aim
    aflr4 = myProblem.analysis.create(aim = "aflr4AIM", name = "surfaceMesh")

    # Set AIM verbosity
    aflr4.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False
    
    # Optional: Explicitly write mesh files
    if writeMesh:
        aflr4.input.Mesh_Format = "Tecplot"
    
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
    
    fun3d = myProblem.analysis.create(aim = "fun3dAIM", name = "fun3d")
    
    fun3d.input["Mesh"].link(myProblem.analysis["volumeMesh"].output["Volume_Mesh"])
    
    # Set project name
    fun3d.input.Proj_Name = "fun3dAFLRTest"
    
    # Set AoA number
    fun3d.input.Alpha = 1.0
    
    # Set Viscous term
    fun3d.input.Viscous = "inviscid"
    
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
    bc1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
    bc2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
    fun3d.input.Boundary_Condition = {"Wing1"   : bc1,
                                      "Wing2"   : bc2,
                                      "Farfield": "farfield"}
    
    # Set flag to indicate mesh morphing is desired
    fun3d.input.Mesh_Morph = True
    fun3d.input.Design_Sensitivity = True
    
    fun3d.input.Design_Variable =  {"Alpha" : {"lowerBound": 0.0}, 
                                    "twist" : {"lowerBound": 0.0}}
    
    fun3d.input.Design_Functional = {"Cd2": {"function":"Cd", "power":2},
                                     "Const": [{"function":"Cl"}, {"function":"Cmz","power":2}]}
    
    return fun3d


def flow_Run(fun3d, cores, run=True, morph=False, redirect=False):
    
    fun3d.preAnalysis()
    
    if cores > 1:
        solver = "mpirun -np " + str(cores) + " nodet_mpi"
    else:
        solver = "nodet"
        
    if morph: 
        options =  " --animation_freq -1 --volume_animation_freq -1 --read_surface_from_file --write_mesh " + \
                   fun3d.input.Proj_Name 
    else: 
        options =  " --animation_freq -1 --volume_animation_freq -1"  
    
    if redirect: 
        options += " > Info.out"
        
    ####### Run fun3d ####################
    if run == True:
        print ("\n\nRunning FUN3D......") 
        
        fun3d.system (solver + options, "Flow")
        

def flow_Result(fun3d, ran=True):
    
    fun3d.postAnalysis()
    
    if ran == True:
        # Get force results
        print ("Total Force - Pressure + Viscous")
        # Get Lift and Drag coefficients
        print ("Cl = " , fun3d.output.CLtot,
               "Cd = " , fun3d.output.CDtot)

#################################

# Working directory
workDir = os.path.join(str(args.workDir[0]), "FUN3D_Morphing")


# Load CSM file and build the geometry explicitly
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=os.path.join("..","csmData","cfdMultiBody.csm"),
                           outLevel=args.outLevel, 
                           phaseName = "Scratch")

if args.noAnalysis:
    run = ran = False
else:
    run = ran = True

# Turn off the wake 
myProblem.geometry.cfgpmtr.wake = False

# Setup
mesh_Setup(myProblem, writeMesh = True)
fun3d = flow_Setup(myProblem)

# Put all surfaces into 1 body file 
fun3d.input.Mesh_Morph_Combine = True;
 
# Run once normally
flow_Run(fun3d, args.cores, run)
flow_Result(fun3d, ran)

# Unlink the mesh to force a morph
fun3d.input["Mesh"].unlink() # Force a morph

# Change geometry 
myProblem.geometry.despmtr.twist = 16

# Re-run with deformed mesh
flow_Run(fun3d, args.cores, run, morph=True)
flow_Result(fun3d, ran)


