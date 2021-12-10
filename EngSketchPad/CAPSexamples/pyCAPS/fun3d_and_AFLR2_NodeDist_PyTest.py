# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D and AFLR2 Node Distance PyTest Example',
                                 prog = 'fun3d_and_AFLR2_NodeDist_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "FUN3DAFLR2NodeDistAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfd2D.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load aflr2 aim
myProblem.analysis.create(aim = "aflr2AIM", name = "aflr2")

airfoil = {"numEdgePoints" : 40, "edgeDistribution" : "Tanh", "initialNodeSpacing" : [0.01, 0.01]}

myProblem.analysis["aflr2"].input.Mesh_Sizing = {"Airfoil"   : airfoil,
                                                 "TunnelWall": {"numEdgePoints" : 50},
                                                 "InFlow"    : {"numEdgePoints" : 25},
                                                 "OutFlow"   : {"numEdgePoints" : 25},
                                                 "2DSlice"   : {"tessParams" : [0.50, .01, 45]}}

# Make quad/tri instead of just quads
#myProblem.analysis["aflr2AIM"].input.Mesh_Gen_Input_String = "mquad=1 mpp=3"

# Load FUN3D aim
myProblem.analysis.create(aim = "fun3dAIM",
                          name = "fun3d")

myProblem.analysis["fun3d"].input["Mesh"].link(myProblem.analysis["aflr2"].output["Area_Mesh"])

# Set project name
myProblem.analysis["fun3d"].input.Proj_Name = "pyCAPS_FUN3D_aflr2"

myProblem.analysis["fun3d"].input.Mesh_ASCII_Flag = True

# Set AoA number
myProblem.analysis["fun3d"].input.Alpha = 0.0

# Set equation type
myProblem.analysis["fun3d"].input.Equation_Type = "compressible"

# Set Viscous term
myProblem.analysis["fun3d"].input.Viscous = "laminar"

# Set the Reynolds number
myProblem.analysis["fun3d"].input.Re = 10

# Set number of iterations
myProblem.analysis["fun3d"].input.Num_Iter = 10

# Set CFL number schedule
myProblem.analysis["fun3d"].input.CFL_Schedule = [1, 40]

# Set read restart option
myProblem.analysis["fun3d"].input.Restart_Read = "off"

# Set CFL number iteration schedule
myProblem.analysis["fun3d"].input.CFL_Schedule_Iter = [1, 200]

# Set overwrite fun3d.nml if not linking to Python library
myProblem.analysis["fun3d"].input.Overwrite_NML = True

# Set the aim to 2D mode
myProblem.analysis["fun3d"].input.Two_Dimensional = True

# Set Mach number
Mach = 0.4
gamma = 1.4

myProblem.analysis["fun3d"].input.Mach = Mach

# Set boundary conditions
wallBC = {"bcType" : "Viscous", "wallTemperature" : 1}

backPressureBC = {"bcType" : "SubsonicOutflow",
                  "staticPressure" : 1}

inflow = {"bcType" : "SubsonicInflow",
          "totalPressure" : (1+ (gamma-1.0)/2.0*Mach**2)**(gamma/(gamma-1)), # Isentropic relation
          "totalTemperature" : (1+ (gamma-1.0)/2.0*Mach**2)}

myProblem.analysis["fun3d"].input.Boundary_Condition = {"Airfoil"   : wallBC,
                                                        "TunnelWall": wallBC,
                                                        "InFlow"    : inflow,
                                                        "OutFlow"   : backPressureBC,
                                                        "2DSlice"   : "SymmetryY"}

# Run AIM pre-analysis
myProblem.analysis["fun3d"].preAnalysis()

####### Run fun3d ####################
print ("\n\nRunning FUN3D......")
myProblem.analysis["fun3d"].system("nodet_mpi --animation_freq -1 > Info.out"); # Run fun3d via system call
#######################################

# Run AIM post-analysis
myProblem.analysis["fun3d"].postAnalysis()

# Get force results
if haveFUN3D:
    print ("Total Force - Pressure + Viscous")
    # Get Lift and Drag coefficients
    print ("Cl = " , myProblem.analysis["fun3d"].output.CLtot,
           "Cd = " , myProblem.analysis["fun3d"].output.CDtot)

    # Get Cmx, Cmy, and Cmz coefficients
    print ("Cmx = " , myProblem.analysis["fun3d"].output.CMXtot, 
           "Cmy = " , myProblem.analysis["fun3d"].output.CMYtot, 
           "Cmz = " , myProblem.analysis["fun3d"].output.CMZtot)

    # Get Cx, Cy, Cz coefficients
    print ("Cx = " , myProblem.analysis["fun3d"].output.CXtot,
           "Cy = " , myProblem.analysis["fun3d"].output.CYtot,
           "Cz = " , myProblem.analysis["fun3d"].output.CZtot)

    print ("Pressure Contribution")
    # Get Lift and Drag coefficients
    print ("Cl_p = " , myProblem.analysis["fun3d"].output.CLtot_p,
           "Cd_p = " , myProblem.analysis["fun3d"].output.CDtot_p)

    # Get Cmx, Cmy, and Cmz coefficients
    print ("Cmx_p = " , myProblem.analysis["fun3d"].output.CMXtot_p,
           "Cmy_p = " , myProblem.analysis["fun3d"].output.CMYtot_p,
           "Cmz_p = " , myProblem.analysis["fun3d"].output.CMZtot_p)

    # Get Cx, Cy, and Cz, coefficients
    print ("Cx_p = " , myProblem.analysis["fun3d"].output.CXtot_p,
           "Cy_p = " , myProblem.analysis["fun3d"].output.CYtot_p,
           "Cz_p = " , myProblem.analysis["fun3d"].output.CZtot_p)

    print ("Viscous Contribution")
    # Get Lift and Drag coefficients
    print ("Cl_v = " , myProblem.analysis["fun3d"].output.CLtot_v,
           "Cd_v = " , myProblem.analysis["fun3d"].output.CDtot_v)

    # Get Cmx, Cmy, and Cmz coefficients
    print ("Cmx_v = " , myProblem.analysis["fun3d"].output.CMXtot_v,
           "Cmy_v = " , myProblem.analysis["fun3d"].output.CMYtot_v,
           "Cmz_v = " , myProblem.analysis["fun3d"].output.CMZtot_v)

    # Get Cx, Cy, and Cz, coefficients
    print ("Cx_v = " , myProblem.analysis["fun3d"].output.CXtot_v,
           "Cy_v = " , myProblem.analysis["fun3d"].output.CYtot_v,
           "Cz_v = " , myProblem.analysis["fun3d"].output.CZtot_v)
