# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FUN3D and AFLR2 PyTest Example',
                                 prog = 'fun3d_and_AFLR2_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

#NREL's S809 Airfoil coordinates
#          x          y
S809 = [1.000000,  0.000000,
        0.996203,  0.000487,
        0.985190,  0.002373,
        0.967844,  0.005960,
        0.945073,  0.011024,
        0.917488,  0.017033,
        0.885293,  0.023458,
        0.848455,  0.030280,
        0.807470,  0.037766,
        0.763042,  0.045974,
        0.715952,  0.054872,
        0.667064,  0.064353,
        0.617331,  0.074214,
        0.567830,  0.084095,
        0.519832,  0.093268,
        0.474243,  0.099392,
        0.428461,  0.101760,
        0.382612,  0.101840,
        0.337260,  0.100070,
        0.292970,  0.096703,
        0.250247,  0.091908,
        0.209576,  0.085851,
        0.171409,  0.078687,
        0.136174,  0.070580,
        0.104263,  0.061697,
        0.076035,  0.052224,
        0.051823,  0.042352,
        0.031910,  0.032299,
        0.016590,  0.022290,
        0.006026,  0.012615,
        0.000658,  0.003723,
        0.000204,  0.001942,
        0.000000, -0.000020,
        0.000213, -0.001794,
        0.001045, -0.003477,
        0.001208, -0.003724,
        0.002398, -0.005266,
        0.009313, -0.011499,
        0.023230, -0.020399,
        0.042320, -0.030269,
        0.065877, -0.040821,
        0.093426, -0.051923,
        0.124111, -0.063082,
        0.157653, -0.073730,
        0.193738, -0.083567,
        0.231914, -0.092442,
        0.271438, -0.099905,
        0.311968, -0.105281,
        0.353370, -0.108181,
        0.395329, -0.108011,
        0.438273, -0.104552,
        0.481920, -0.097347,
        0.527928, -0.086571,
        0.576211, -0.073979,
        0.626092, -0.060644,
        0.676744, -0.047441,
        0.727211, -0.035100,
        0.776432, -0.024204,
        0.823285, -0.015163,
        0.866630, -0.008204,
        0.905365, -0.003363,
        0.938474, -0.000487,
        0.965086,  0.000743,
        0.984478,  0.000775,
        0.996141,  0.000290,
        1.000000,  0.000000]

# Shape design variable is 1x200 need to fill in the rest of the array with "unused values"
for i in range(0, 200-len(S809)):
    S809.append(1.0E6)

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "FUN3DAFLR2ArbShapeAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","cfd2DArbShape.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.xy = S809

print ("Saving geometry")
myProblem.geometry.save("GenericShape", directory = workDir, extension = "egads")

# Load aflr2 aim
myProblem.analysis.create(aim = "aflr2AIM", name = "aflr2")

myProblem.analysis["aflr2"].input.Mesh_Sizing = {"Airfoil"   : {"numEdgePoints" : 400},
                                                 "TunnelWall": {"numEdgePoints" : 50},
                                                 "InFlow"    : {"numEdgePoints" : 25},
                                                 "OutFlow"   : {"numEdgePoints" : 25},
                                                 "2DSlice"   : {"tessParams" : [0.50, .01, 45]}}

# Make quad/tri instead of just quads
myProblem.analysis["aflr2"].input.Mesh_Gen_Input_String = "mquad=1 mpp=3"

# Load FUN3D aim - child of aflr2 AIM
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
