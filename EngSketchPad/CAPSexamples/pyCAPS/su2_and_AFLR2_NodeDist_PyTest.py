# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Import SU2 Python interface module
from parallel_computation import parallel_computation as su2Run

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'SU2 and AFLR2 Node Distance PyTest Example',
                                 prog = 'su2_and_AFLR2_NodeDist_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Project name
projectName = "pyCAPS_SU2_aflr2"

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "SU2AFLR2NodeDistAnalysisTest")

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

# NO analysis is needed - AFLR2 executes automatically

# Load SU2 aim
su2 = myProblem.analysis.create(aim = "su2AIM",
                                name = "su2")

su2.input["Mesh"].link(myProblem.analysis["aflr2"].output["Area_Mesh"])

# Set project name
su2.input.Proj_Name = projectName

# Set SU2 Version
su2.input.SU2_Version = "Blackbird"

# Set AoA number
su2.input.Alpha = 0.0

# Set equation type
su2.input.Equation_Type = "Compressible"

# Set number of iterations
su2.input.Num_Iter = 10

# Set the aim to 2D mode
su2.input.Two_Dimensional = True

# Set output file format
su2.input.Output_Format = "Tecplot"

# Set boundary conditions
inviscidBC = {"bcType" : "Inviscid"}

# Set Mach number
Mach = 0.4
gamma = 1.4
Pinf = 101300.0
Tinf = 280

# Set Mach number
su2.input.Mach = Mach

# Set boundary conditions
wallBC = {"bcType" : "Inviscid"}

backPressureBC = {"bcType" : "SubsonicOutflow",
                  "staticPressure" : Pinf}

inflow = {"bcType" : "SubsonicInflow",
          "totalPressure" : Pinf*(1+ (gamma-1.0)/2.0*Mach**2)**(gamma/(gamma-1)), # Isentropic relation
          "totalTemperature" : Tinf*(1+ (gamma-1.0)/2.0*Mach**2),
          "uVelocity" : 1.0,
          "vVelocity" : 0.0,
          "wVelocity" : 0.0}

su2.input.Boundary_Condition = {"Airfoil"   : wallBC,
                                "TunnelWall": wallBC,
                                "InFlow"    : inflow,
                                "OutFlow"   : backPressureBC}

# Specifcy the boundares used to compute forces
su2.input.Surface_Monitor = ["Airfoil"]

# Run AIM pre-analysis
su2.preAnalysis()

####### Run su2 ####################
print ("\n\nRunning SU2......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(su2.analysisDir) # Move into test directory

su2Run(projectName + ".cfg", args.numberProc); # Run SU2

os.chdir(currentDirectory) # Move back to top directory
#######################################

# Run AIM post-analysis
su2.postAnalysis()

print ("Total Force - Pressure + Viscous")
# Get Lift and Drag coefficients
print ("Cl = " , myProblem.analysis["su2"].output.CLtot,
       "Cd = " , myProblem.analysis["su2"].output.CDtot)

# Get Cmz coefficient
print ("Cmz = " , myProblem.analysis["su2"].output.CMZtot)

# Get Cx and Cy coefficients
print ("Cx = " , myProblem.analysis["su2"].output.CXtot,
       "Cy = " , myProblem.analysis["su2"].output.CYtot)

print ("Pressure Contribution")
# Get Lift and Drag coefficients
print ("Cl_p = " , myProblem.analysis["su2"].output.CLtot_p,
       "Cd_p = " , myProblem.analysis["su2"].output.CDtot_p)

# Get Cmz coefficient
print ("Cmz_p = " , myProblem.analysis["su2"].output.CMZtot_p)

# Get Cx and Cy coefficients
print ("Cx_p = " , myProblem.analysis["su2"].output.CXtot_p,
       "Cy_p = " , myProblem.analysis["su2"].output.CYtot_p)

print ("Viscous Contribution")
# Get Lift and Drag coefficients
print ("Cl_v = " , myProblem.analysis["su2"].output.CLtot_v,
       "Cd_v = " , myProblem.analysis["su2"].output.CDtot_v)

# Get Cmz coefficient
print ("Cmz_v = " , myProblem.analysis["su2"].output.CMZtot_v)

# Get Cx and Cy coefficients
print ("Cx_v = " , myProblem.analysis["su2"].output.CXtot_v,
       "Cy_v = " , myProblem.analysis["su2"].output.CYtot_v)
