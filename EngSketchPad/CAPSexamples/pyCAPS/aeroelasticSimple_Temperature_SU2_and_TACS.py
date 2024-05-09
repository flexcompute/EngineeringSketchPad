# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import SU2 python environment
from parallel_computation import parallel_computation as su2Run

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Aeroelastic Temperature SU2 and TACS Example',
                                 prog = 'aeroelasticSimple_Temperature_SU2_and_TACS',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AeroelasticSimple_Temperature_SU2_and_TACS")

# Create projectName vairbale
projectName = "aeroelasticSimple_Temperature_SM"

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransferSimple.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
surfMesh = myProblem.analysis.create(aim = "aflr4AIM",
                                     name= "aflr4",
                                     capsIntent = "Aerodynamic")

mesh = myProblem.analysis.create(aim = "aflr3AIM",
                                 name= "aflr3",
                                 capsIntent = "Aerodynamic")

mesh.input["Surface_Mesh"].link(surfMesh.output["Surface_Mesh"])

su2 = myProblem.analysis.create(aim = "su2AIM",
                                name = "su2",
                                capsIntent = "Aerodynamic")

su2.input["Mesh"].link(mesh.output["Volume_Mesh"])

tacs = myProblem.analysis.create(aim = "tacsAIM",
                                    name = "tacs",
                                    capsIntent = "Structure",
                                    autoExec = False)

# Create EGADS tess aim
tess = myProblem.analysis.create(aim = "egadsTessAIM",
                                 name = "tess",
                                 capsIntent = "Structure")

# No Tess vertexes on edges (minimial mesh)
tess.input.Edge_Point_Min = 3
tess.input.Edge_Point_Max = 3

# Use regularized quads
tess.input.Mesh_Elements = "Mixed"

tacs.input["Mesh"].link(tess.output["Surface_Mesh"])

# Create the data transfer connections
boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)

    # Create the vertex sets on the bound for su2 and tacs analysis
    su2Vset     = bound.vertexSet.create(su2)
    tacsVset = bound.vertexSet.create(tacs)

    # Create displacement data sets
    su2_Temperature     = su2Vset.dataSet.create("Temperature")
    tacs_Temperature = tacsVset.dataSet.create("Temperature")

    # Link the data set
    tacs_Temperature.link(su2_Temperature, "Conserve")

    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()


# Farfield growth factor
surfMesh.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
surfMesh.input.max_scale = 0.75
surfMesh.input.min_scale = 0.1

# Set inputs for volume mesh
mesh.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set inputs for su2
speedofSound = 340.0 # m/s
refVelocity = 100.0 # m/s
refDensity = 1.2 # kg/m^3

su2.input.Proj_Name = projectName
su2.input.Mach = refVelocity/speedofSound
su2.input.Equation_Type = "compressible"
su2.input.Num_Iter = 3
su2.input.Output_Format = "Tecplot"
su2.input.SU2_Version = "Blackbird"
su2.input.Temperature_Scale_Factor = 1
su2.input.Surface_Monitor = ["Skin"]

inviscid = {"bcType" : "Inviscid"}
su2.input.Boundary_Condition = {"Skin"     : inviscid,
                                "SymmPlane": "SymmetryY",
                                "Farfield" : "farfield"}

# Set inputs for tacs
tacs.input.Proj_Name = projectName
tacs.input.Analysis_Type = "Static"

tacs.input.File_Format = "Free"

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio" : 0.33,
                "density" : 2.8E3}
tacs.input.Material = {"Madeupium": madeupium}

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

tacs.input.Property = {"Skin"    : skin,
                          "Rib_Root": ribSpar}

constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
tacs.input.Constraint = {"edgeConstraint": constraint}

# External Temperature load to tacs that we will inherited from su2
load = {"loadType" : "ThermalExternal"}
tacs.input.Load = {"Skin": load}


####### SU2 ###########################
# Re-run the preAnalysis
print ("\nRunning PreAnalysis ......", "su2")
su2.preAnalysis()

# Run su2
print ("\n\nRunning SU2......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(su2.analysisDir) # Move into test directory

try:
    su2Run(su2.input.Proj_Name + ".cfg", args.numberProc) # Run SU2 CFD solver
except:
    raise SystemError("SU2 FAILED!")

os.chdir(currentDirectory) # Move back to top directory

print ("\nRunning PostAnalysis ......", "su2")
# Run AIM post-analysis
su2.postAnalysis()
#######################################

# Plot the dataTransfer
for boundName in boundNames:
    if args.noPlotData == False:
        try:
            print ("\tPlotting dataTransfer source......", boundName)
            myProblem.bound[boundName].vertexSet["su2"].dataSet["Temperature"].view()
            print ("\tPlotting dataTransfer destination......")
            myProblem.bound[boundName].vertexSet["tacs"].dataSet["Temperature"].view()
        except ImportError:
            print("Unable to plot data")


####### TACS #######################
# Run pre/post-analysis for tacs and execute
print ("\nRunning PreAnalysis ......", "tacs")
tacs.preAnalysis()

# Run tacs...

print ("\nRunning PostAnalysis ......", "tacs")
tacs.postAnalysis()
#######################################
