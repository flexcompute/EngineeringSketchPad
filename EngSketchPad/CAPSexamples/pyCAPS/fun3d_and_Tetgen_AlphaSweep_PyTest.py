# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D and Tetgen Alpha Sweep PyTest Example',
                                 prog = 'fun3d_and_Tetgen_AlphaSweep_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "FUN3DTetgenAlphaSweep")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdMultiBody.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load egadsTess aim
myProblem.analysis.create(aim = "egadsTessAIM", name = "egadsTess")

# Set new EGADS body tessellation parameters
myProblem.analysis["egadsTess"].input.Tess_Params = [1.0, 0.01, 20.0]


# Load Tetgen aim
myProblem.analysis.create(aim = "tetgenAIM", name = "tetgen")

# Link Surface_Mesh
myProblem.analysis["tetgen"].input["Surface_Mesh"].link(myProblem.analysis["egadsTess"].output["Surface_Mesh"])

# Preserve surface mesh while messhing
myProblem.analysis["tetgen"].input.Preserve_Surf_Mesh = True

# Get analysis info.
myProblem.analysis["tetgen"].info()

Alpha = [-5, 0, 5, 10]
Cl = []
Cd = []
for a in Alpha:
    analysisID = "Alpha_" + str(a)
    # Load FUN3D aim - child of Tetgen AIM
    myProblem.analysis.create(aim = "fun3dAIM",
                              name = analysisID)

    # Link the mesh
    myProblem.analysis[analysisID].input["Mesh"].link(myProblem.analysis["tetgen"].output["Volume_Mesh"])

    # Set project name
    myProblem.analysis[analysisID].input.Proj_Name = analysisID

    myProblem.analysis[analysisID].input.Mesh_ASCII_Flag = False

    # Set AoA number
    myProblem.analysis[analysisID].input.Alpha = a

    # Set Mach number
    myProblem.analysis[analysisID].input.Mach = 0.25

    # Set equation type
    myProblem.analysis[analysisID].input.Equation_Type = "compressible"

    # Set Viscous term
    myProblem.analysis[analysisID].input.Viscous = "inviscid"

    # Set number of iterations
    myProblem.analysis[analysisID].input.Num_Iter = 10

    # Set CFL number schedule
    myProblem.analysis[analysisID].input.CFL_Schedule = [0.5, 3.0]

    # Set read restart option
    myProblem.analysis[analysisID].input.Restart_Read = "off"

    # Set CFL number iteration schedule
    myProblem.analysis[analysisID].input.CFL_Schedule_Iter = [1, 40]

    # Set overwrite fun3d.nml if not linking to Python library
    myProblem.analysis[analysisID].input.Overwrite_NML = True

    # Set boundary conditions
    inviscidBC1 = {"bcType" : "Inviscid", "wallTemperature" : 1}
    inviscidBC2 = {"bcType" : "Inviscid", "wallTemperature" : 1.2}
    myProblem.analysis[analysisID].input.Boundary_Condition = {"Wing1"   : inviscidBC1,
                                                               "Wing2"   : inviscidBC2,
                                                               "Farfield": "farfield"}

    # Run AIM pre-analysis
    myProblem.analysis[analysisID].preAnalysis()

    ####### Run fun3d ####################
    print ("\n\nRunning FUN3D......")
    myProblem.analysis[analysisID].system("nodet_mpi --animation_freq -1000 > Info.out"); # Run fun3d via system call
    #######################################

    # Run AIM post-analysis
    myProblem.analysis[analysisID].postAnalysis()

    # Get force results
    if haveFUN3D:
        Cl.append(myProblem.analysis[analysisID].output.CLtot)
        Cd.append(myProblem.analysis[analysisID].output.CDtot)

print ("Alpha, CL, CD")
for i in range(len(Alpha)):
    print (Alpha[i], Cl[i], Cd[i])
