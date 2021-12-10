# Import pyCAPS module
import pyCAPS

# Import modules
import os
import sys
import shutil
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Masstran AGARD445.6 Pytest Example',
                                 prog = 'masstran_AGARD445_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "masstranAGARD445")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaAGARD445.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
myAnalysis = myProblem.analysis.create(aim = "masstranAIM",
                                       name = "masstran")

# Set meshing parameters
myAnalysis.input.Edge_Point_Max = 10
myAnalysis.input.Edge_Point_Min = 6

myAnalysis.input.Quad_Mesh = True

myAnalysis.input.Tess_Params = [.25,.01,15]

# Set materials
mahogany    = {"materialType"        : "orthotropic",
               "youngModulus"        : 0.457E6 ,
               "youngModulusLateral" : 0.0636E6,
               "poissonRatio"        : 0.31,
               "shearModulus"        : 0.0637E6,
               "shearModulusTrans1Z" : 0.00227E6,
               "shearModulusTrans2Z" : 0.00227E6,
               "density"             : 3.5742E-5}

myAnalysis.input.Material = {"Mahogany": mahogany}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.82,
          "material"        : "mahogany",
          "bendingInertiaRatio" : 1.0, # Default - not necesssary
          "shearMembraneRatio"  : 5.0/6.0} # Default - not necesssary

myAnalysis.input.Property = {"yatesPlate": shell}

#####################################
## masstran executes automatically ##
#####################################

# Get mass properties
print ("\nGetting results mass properties.....\n")
Area     = myAnalysis.output.Area
Mass     = myAnalysis.output.Mass
Centroid = myAnalysis.output.Centroid
CG       = myAnalysis.output.CG
Ixx      = myAnalysis.output.Ixx
Iyy      = myAnalysis.output.Iyy
Izz      = myAnalysis.output.Izz
Ixy      = myAnalysis.output.Ixy
Ixz      = myAnalysis.output.Ixz
Iyz      = myAnalysis.output.Iyz
I        = myAnalysis.output.I_Vector
II       = myAnalysis.output.I_Tensor

print("Area     ", Area)
print("Mass     ", Mass)
print("Centroid ", Centroid)
print("CG       ", CG)
print("Ixx      ", Ixx)
print("Iyy      ", Iyy)
print("Izz      ", Izz)
print("Ixy      ", Ixy)
print("Ixz      ", Ixz)
print("Iyz      ", Iyz)
print("I        ", I)
print("II       ", II)
