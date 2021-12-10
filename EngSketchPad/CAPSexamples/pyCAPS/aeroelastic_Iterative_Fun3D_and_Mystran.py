# Import pyCAPS module
import pyCAPS

# Import os module
import os
import shutil
import math

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Aeroelastic Fun3D and Mystran Example',
                                 prog = 'aeroelastic_Iterative_Fun2_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "Aeroelastic_Iterative_FM")

# Create projectName vairbale
projectName = "aeroelastic_Iterative"

# Set the number of transfer iterations
numTransferIteration = 2

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransfer.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
myProblem.analysis.create(aim = "egadsTessAIM",
                          name= "egads",
                          capsIntent = "CFD")

myProblem.analysis.create(aim = "tetgenAIM",
                          name= "tetgen",
                          capsIntent = "CFD")

myProblem.analysis["tetgen"].input["Surface_Mesh"].link(myProblem.analysis["egads"].output["Surface_Mesh"])

myProblem.analysis.create(aim = "fun3dAIM",
                          name = "fun3d",
                          capsIntent = "CFD")

myProblem.analysis["fun3d"].input["Mesh"].link(myProblem.analysis["tetgen"].output["Volume_Mesh"])

myProblem.analysis.create(aim = "mystranAIM",
                          name = "mystran",
                          capsIntent = "STRUCTURE",
                          autoExec = False)

boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)
    
    # Create the vertex sets on the bound for fun3d and mystran analysis
    fun3dVset   = bound.vertexSet.create(myProblem.analysis["fun3d"])
    mystranVset = bound.vertexSet.create(myProblem.analysis["mystran"])
    
    # Create pressure data sets
    fun3d_Pressure   = fun3dVset.dataSet.create("Pressure", pyCAPS.fType.FieldOut)
    mystran_Pressure = mystranVset.dataSet.create("Pressure", pyCAPS.fType.FieldIn)

    # Create displacement data sets
    fun3d_Displacement   = fun3dVset.dataSet.create("Displacement", pyCAPS.fType.FieldIn, init=[0,0,0])
    mystran_Displacement = mystranVset.dataSet.create("Displacement", pyCAPS.fType.FieldOut)

    # Link the data sets
    mystran_Pressure.link(fun3d_Pressure, "Conserve")
    fun3d_Displacement.link(mystran_Displacement, "Interpolate")
    
    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()

# Set inputs for tetgen
myProblem.analysis["tetgen"].input.Tess_Params = [.3, 0.01, 20.0]
myProblem.analysis["tetgen"].input.Preserve_Surf_Mesh = True

# Set inputs for fun3d
Mach = 0.3
Rair = 1716 #ft^2/ (s^2* R)
refTemperature = 491.4 # R
gamma = 1.4
speedofSound = math.sqrt(gamma*refTemperature*Rair) # ft/s

refVelocity = speedofSound*Mach # ft/s
refDensity = 0.0023366 # slug/ft^3
refPressure = (refDensity*speedofSound**2)/gamma

AoA = 10 # angle of attack

myProblem.analysis["fun3d"].input.Proj_Name = projectName
myProblem.analysis["fun3d"].input.Mesh_ASCII_Flag = False
myProblem.analysis["fun3d"].input.Mach = Mach
myProblem.analysis["fun3d"].input.Alpha = AoA
myProblem.analysis["fun3d"].input.Equation_Type = "compressible"
myProblem.analysis["fun3d"].input.Viscous = "inviscid"
myProblem.analysis["fun3d"].input.Num_Iter = 1000
myProblem.analysis["fun3d"].input.CFL_Schedule = [1, 5.0]
myProblem.analysis["fun3d"].input.CFL_Schedule_Iter = [1, 40]
myProblem.analysis["fun3d"].input.Overwrite_NML = True
myProblem.analysis["fun3d"].input.Restart_Read = "off"
myProblem.analysis["fun3d"].input.Pressure_Scale_Factor = 0.5*refDensity*refVelocity**2

inviscid = {"bcType" : "Inviscid", "wallTemperature" : 1.1}
myProblem.analysis["fun3d"].input.Boundary_Condition = {"Skin"      : inviscid,
                                                         "SymmPlane": "SymmetryY",
                                                         "Farfield" : "farfield"}

# Set inputs for mystran
myProblem.analysis["mystran"].input.Proj_Name = projectName
myProblem.analysis["mystran"].input.Edge_Point_Max = 3

myProblem.analysis["mystran"].input.Quad_Mesh = True
myProblem.analysis["mystran"].input.Tess_Params = [.5, .1, 15]
myProblem.analysis["mystran"].input.Analysis_Type = "Static"

# External pressure load to mystran that we will inherited from fun3d
load = {"loadType" : "PressureExternal", "loadScaleFactor" : -1.0}
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
         "membraneThickness" : 0.0052, # ft - 1/16 in
         "material"        : "madeupium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

ribSpar  = {"propertyType" : "Shell",
            "membraneThickness" : 0.0104, # ft  - 0.25 in
            "material"        : "unobtainium",
            "bendingInertiaRatio" : 1.0, # Default
            "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["mystran"].input.Property = {"Skin": skin,
                                                "Ribs_and_Spars": ribSpar,
                                                "Rib_Root": ribSpar}
constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
myProblem.analysis["mystran"].input.Constraint = {"edgeConstraint": constraint}


# Aeroelastic iteration loop
for iter in range(numTransferIteration):

    ####### FUN3D #########################
    print ("\nRunning PreAnalysis ......", "fun3d")
    myProblem.analysis["fun3d"].preAnalysis()

    #------------------------------
    print ("\n\nRunning FUN3D......")
    analysisDir = myProblem.analysis["fun3d"].analysisDir
    
    #--write_aero_loads_to_file --aeroelastic_external
    cmdLineOpt = "--write_aero_loads_to_file --animation_freq -1"
    if iter != 0:
        cmdLineOpt = cmdLineOpt + " --read_surface_from_file"

    myProblem.analysis["fun3d"].system("mpirun -np 10 nodet_mpi " + cmdLineOpt + " > Info.out"); # Run fun3d via system call

    if os.path.getsize(os.path.join(analysisDir, "Info.out")) == 0: #
        raise SystemError("FUN3D excution failed\n")

    shutil.copy(os.path.join(analysisDir, projectName + "_tec_boundary.plt"), os.path.join(analysisDir, projectName + "_tec_boundary" + "_Interation_" + str(iter) + ".plt"))
    #------------------------------

    print ("\nRunning PostAnalysis ......", "fun3d")
    # Run AIM post-analysis
    myProblem.analysis["fun3d"].postAnalysis()
    #######################################


    ####### Mystran #######################
    print ("\nRunning PreAnalysis ......", "mystran")
    myProblem.analysis["mystran"].preAnalysis()

    #------------------------------
    print ("\n\nRunning Mystran......")
    analysisDir = myProblem.analysis["mystran"].analysisDir

    myProblem.analysis["mystran"].system("mystran.exe " + projectName +  ".dat > Info.out") # Run mystran via system call

    if os.path.getsize(os.path.join(analysisDir, "Info.out")) == 0: #
        raise SystemError("Mystran excution failed\n")
    #------------------------------

    print ("\nRunning PostAnalysis ......", "mystran")
    myProblem.analysis["mystran"].postAnalysis()
    #######################################
