# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Aeroelastic Displacement Fun3D and Mystran Example',
                                 prog = 'aeroelasticSimple_Displacement_Fun3D_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AeroelasticSimple_Displacement_fun3d_Mystran")

# Create projectName vairbale
projectName = "aeroelasticSimple_Displacement_Fun3D_Mystran"

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransferSimple.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
surfMesh = myProblem.analysis.create(aim = "egadsTessAIM", 
                                     name= "egads",
                                     capsIntent = "CFD")

mesh = myProblem.analysis.create(aim = "tetgenAIM", 
                                 name= "tetgen",
                                 capsIntent = "CFD")

mesh.input["Surface_Mesh"].link(surfMesh.output["Surface_Mesh"])

fun3d = myProblem.analysis.create(aim = "fun3dAIM", 
                                  name = "fun3d", 
                                  capsIntent = "CFD")

fun3d.input["Mesh"].link(mesh.output["Volume_Mesh"])

mystran = myProblem.analysis.create(aim = "mystranAIM",
                                    name = "mystran",
                                    capsIntent = "STRUCTURE",
                                    autoExec = False)

# Create the data transfer connections
boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)
    
    # Create the vertex sets on the bound for fun3d and mystran analysis
    fun3dVset   = bound.vertexSet.create(fun3d)
    mystranVset = bound.vertexSet.create(mystran)

    # Create displacement data sets
    fun3d_Displacement   = fun3dVset.dataSet.create("Displacement", pyCAPS.fType.FieldIn)
    mystran_Displacement = mystranVset.dataSet.create("Displacement", pyCAPS.fType.FieldOut)

    # Link the data set
    fun3d_Displacement.link(mystran_Displacement, "Interpolate")
    
    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()

# Set inputs for egads 
surfMesh.input.Tess_Params = [.05, 0.01, 20.0]

# Set inputs for tetgen 
mesh.input.Preserve_Surf_Mesh = True
mesh.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set inputs for fun3d
speedofSound = 340.0 # m/s
refVelocity = 100.0 # m/s
refDensity = 1.2 # kg/m^3

fun3d.input.Proj_Name = projectName
fun3d.input.Mesh_ASCII_Flag = False
fun3d.input.Mach = refVelocity/speedofSound
fun3d.input.Equation_Type = "compressible"
fun3d.input.Viscous = "inviscid"
fun3d.input.Num_Iter = 10
fun3d.input.CFL_Schedule = [1, 5.0]
fun3d.input.CFL_Schedule_Iter = [1, 40]
fun3d.input.Overwrite_NML = True
fun3d.input.Restart_Read = "off"

inviscid = {"bcType" : "Inviscid", "wallTemperature" : 1.1}
fun3d.input.Boundary_Condition = {"Skin"     : inviscid,
                                  "SymmPlane": "SymmetryY",
                                  "Farfield" : "farfield"}

# Set inputs for mystran
mystran.input.Proj_Name = projectName
mystran.input.Edge_Point_Max = 3
mystran.input.Edge_Point_Min = 3

mystran.input.Quad_Mesh = True
mystran.input.Tess_Params = [.5, .1, 15]
mystran.input.Analysis_Type = "Static"

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio" : 0.33,
                "density" : 2.8E3}
mystran.input.Material = {"Madeupium": madeupium}

skin  = {"propertyType" : "Shell",
         "membraneThickness" : 0.06,
         "material"        : "madeupium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

ribSpar  = {"propertyType" : "Shell",
            "membraneThickness" : 0.6,
            "material"        : "madeupium",
            "bendingInertiaRatio" : 1.0, # Default
            "shearMembraneRatio"  : 5.0/6.0} # Default

mystran.input.Property = {"Skin"    : skin,
                          "Rib_Root": ribSpar}

constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
mystran.input.Constraint = {"edgeConstraint": constraint}

# Static uniform pressure load
load = {"groupName": "Skin", "loadType" : "Pressure", "pressureForce": 1E8}
mystran.input.Load = {"pressureInternal": load}


####### Mystran #######################
# Run pre/post-analysis for mystran and execute
print ("\nRunning PreAnalysis ......", "mystran")
mystran.preAnalysis()

# Run mystran
print ("\n\nRunning Mystran......")
mystran.system("mystran.exe " + projectName +  ".dat > Info.out") # Run fun3d via system call

if os.path.getsize(os.path.join(mystran.analysisDir,"Info.out")) == 0:
    raise SystemError("Mystran excution failed\n")

print ("\nRunning PostAnalysis ......", "mystran")
mystran.postAnalysis()
#######################################

# Plot the dataTransfer
for boundName in boundNames:
    if args.noPlotData == False:
        try:
            print ("\tPlotting dataTransfer source......", boundName)
            myProblem.bound[boundName].vertexSet["mystran"].dataSet["Displacement"].view()
            print ("\tPlotting dataTransfer destination......")
            myProblem.bound[boundName].vertexSet["fun3d"].dataSet["Displacement"].view()
        except ImportError:
            print("Unable to plot data")

####### Run fun3d ####################
# Re-run the preAnalysis 
print ("\nRunning PreAnalysis ......", "fun3d")
fun3d.preAnalysis()

print ("\n\nRunning FUN3D......")  
cmdLineOpt = "--read_surface_from_file --animation_freq -1"

fun3d.system("mpirun -np 5 nodet_mpi " + cmdLineOpt + " > Info.out"); # Run fun3d via system call
    
if os.path.getsize(os.path.join(fun3d.analysisDir,"Info.out")) == 0: # 
    raise SystemError("FUN3D excution failed\n")

print ("\nRunning PostAnalysis ......", "fun3d")
# Run AIM post-analysis
fun3d.postAnalysis()
#######################################
