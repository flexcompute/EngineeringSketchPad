# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import shutil module
import shutil

# Import SU2 python environment
from parallel_computation import parallel_computation as su2Run
from mesh_deformation import mesh_deformation as su2Deform

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Aeroelastic SU2 and Mystran Example',
                                 prog = 'aeroelastic_Iterative_SU2_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "Aeroelastic_Iterative_SM")

# Create projectName vairbale
projectName = "aeroelastic_Iterative"

# Set the number of transfer iterations
numTransferIteration = 4

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransfer.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
myProblem.analysis.create(aim = "egadsTessAIM",
                          name= "egads",
                          capsIntent = "Aerodynamic")

myProblem.analysis.create(aim = "tetgenAIM",
                          name= "tetgen",
                          capsIntent = "Aerodynamic")

myProblem.analysis["tetgen"].input["Surface_Mesh"].link(myProblem.analysis["egads"].output["Surface_Mesh"])

myProblem.analysis.create(aim = "su2AIM",
                          name = "su2",
                          capsIntent = "Aerodynamic")

myProblem.analysis["su2"].input["Mesh"].link(myProblem.analysis["tetgen"].output["Volume_Mesh"])

myProblem.analysis.create(aim = "mystranAIM",
                          name = "mystran",
                          capsIntent = "Structure",
                          autoExec = True)

# Create the data transfer connections
boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)
    
    # Create the vertex sets on the bound for su2 and mystran analysis
    su2Vset     = bound.vertexSet.create(myProblem.analysis["su2"])
    mystranVset = bound.vertexSet.create(myProblem.analysis["mystran"])
    
    # Create pressure data sets
    su2_Pressure     = su2Vset.dataSet.create("Pressure")
    mystran_Pressure = mystranVset.dataSet.create("Pressure")

    # Create displacement data sets
    su2_Displacement     = su2Vset.dataSet.create("Displacement", init=[0,0,0])
    mystran_Displacement = mystranVset.dataSet.create("Displacement")

    # Link the data sets
    mystran_Pressure.link(su2_Pressure, "Conserve")
    su2_Displacement.link(mystran_Displacement, "Interpolate")
    
    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()

# Set inputs for EGADS
myProblem.analysis["egads"].input.Tess_Params = [.6, 0.05, 20.0]
myProblem.analysis["egads"].input.Edge_Point_Max = 4

# Set inputs for tetgen
myProblem.analysis["tetgen"].input.Preserve_Surf_Mesh = True
myProblem.analysis["tetgen"].input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set inputs for su2
speedofSound = 340.0 # m/s
refVelocity = 100.0 # m/s
refDensity = 1.2 # kg/m^3

myProblem.analysis["su2"].input.Proj_Name             = projectName
myProblem.analysis["su2"].input.SU2_Version           = "Harrier"
myProblem.analysis["su2"].input.Mach                  = refVelocity/speedofSound
myProblem.analysis["su2"].input.Equation_Type         = "compressible"
myProblem.analysis["su2"].input.Num_Iter              = 3 # Way too few to converge the solver, but this is an example
myProblem.analysis["su2"].input.Output_Format         = "Tecplot"
myProblem.analysis["su2"].input.Overwrite_CFG         = True
myProblem.analysis["su2"].input.Pressure_Scale_Factor = 0.5*refDensity*refVelocity**2
myProblem.analysis["su2"].input.Surface_Monitor       = "Skin"

inviscid = {"bcType" : "Inviscid"}
myProblem.analysis["su2"].input.Boundary_Condition = {"Skin"     : inviscid,
                                                      "SymmPlane": "SymmetryY",
                                                      "Farfield" : "farfield"}

# Set inputs for mystran
myProblem.analysis["mystran"].input.Proj_Name = projectName
myProblem.analysis["mystran"].input.Edge_Point_Min = 3
myProblem.analysis["mystran"].input.Edge_Point_Max = 10

myProblem.analysis["mystran"].input.Quad_Mesh = True
myProblem.analysis["mystran"].input.Tess_Params = [.5, .1, 15]
myProblem.analysis["mystran"].input.Analysis_Type = "Static"

# External pressure load to mystran that we will inherited from su2
load = {"loadType" : "PressureExternal"}
myProblem.analysis["mystran"].input.Load = {"pressureAero": load}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.5E9, # lbf/ft^2
                "poissonRatio": 0.33,
                "density" : 5.5} #slug/ft^3

unobtainium  = {"materialType" : "isotropic",
                "youngModulus" : 1.0E10, # lbf/ft^2
                "poissonRatio": 0.33,
                "density" : 6.0} #slug/ft^3

myProblem.analysis["mystran"].input.Material = {"Madeupium": madeupium,
                                                "Unobtainium": unobtainium}

skin  = {"propertyType" : "Shell",
         "membraneThickness" : 0.06, # ft - 1/16 in
         "material"        : "madeupium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

ribSpar  = {"propertyType" : "Shell",
            "membraneThickness" : 0.6, # ft  - 0.25 in
            "material"        : "unobtainium",
            "bendingInertiaRatio" : 1.0, # Default
            "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["mystran"].input.Property = {"Skin"          : skin,
                                                "Ribs_and_Spars": ribSpar,
                                                "Rib_Root"      : ribSpar}

constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
myProblem.analysis["mystran"].input.Constraint = {"edgeConstraint": constraint}


# Aeroelastic iteration loop
for iter in range(numTransferIteration):

    ####### SU2 ###########################
    print ("\nRunning PreAnalysis ......", "su2")
    myProblem.analysis["su2"].preAnalysis()

    #------------------------------
    print ("\n\nRunning SU2......")
    currentDirectory = os.getcwd() # Get our current working directory

    os.chdir(myProblem.analysis["su2"].analysisDir) # Move into test directory

    # Run SU2 mesh deformation
    if iter > 0:
        su2Deform(projectName + ".cfg", args.numberProc)

    # Run SU2 solver
    su2Run(projectName + ".cfg", args.numberProc)

    shutil.copy("surface_flow_" + projectName + ".szplt", "surface_flow_" + projectName + "_Iteration_" + str(iter) + ".szplt")
    os.chdir(currentDirectory) # Move back to top directory
    #------------------------------

    print ("\nRunning PostAnalysis ......", "su2")
    myProblem.analysis["su2"].postAnalysis()
    #######################################


    ####### Mystran #######################
    #
    # Mystran executes automatically 
    #
    #######################################
