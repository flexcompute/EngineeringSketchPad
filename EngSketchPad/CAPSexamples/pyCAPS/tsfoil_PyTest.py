# Import pyCAPS module
import pyCAPS

import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'tsFoil Pytest Example',
                                 prog = 'tsfoil_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "tsFoilAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","airfoilSection.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change a design parameter - area in the geometry
myProblem.geometry.despmtr.camber = 0.05

# Load tsfoil aim
tsfoil = myProblem.analysis.create(aim = "tsfoilAIM",
                                   name = "tsfoil")

# Set Mach number, Reynolds number
tsfoil.input.Mach = .6 
#tsfoil.input.Re = 1.0E6 

# Set custom AoA
tsfoil.input.Alpha = 0.0

# tsfoil is automatically executed

print ("Getting results")
# Retrieve Cl and Cd
Cl = tsfoil.output.CL
print("Cl = ", Cl)

Cd = tsfoil.output.CD
print("Cd = ", Cd)

# Wave drag
Cd_Wave = tsfoil.output.CD_Wave
print("Cd Wave = ", Cd_Wave)

# Momement coeff
Cm = tsfoil.output.CM
print("Cm = ", Cm)

# Critical Cp
Cp = tsfoil.output.Cp_Critical
print("Cp _Critcal = ", Cp)
