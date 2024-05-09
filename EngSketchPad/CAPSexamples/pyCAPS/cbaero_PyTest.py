# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'CBAero Pytest Example',
                                 prog = 'cbaero_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
projectName = "CBAeroAnalysisTest"
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file 
geometryScript = os.path.join("..","csmData","cfdX43a.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Turn of the farfield domain
myProblem.geometry.cfgpmtr.makeBoundBox = 0

# Load egadsTess aim to create surface mesh  
meshAIM = myProblem.analysis.create(aim = "egadsTessAIM",
                                    name = "egadsTess" )

# Load CBAero aim 
cbaero = myProblem.analysis.create(aim = "cbaeroAIM", 
                                   name = "cbaero")

deg = pyCAPS.Unit("degree")
psi = pyCAPS.Unit("psi")

# Link mesh
cbaero.input["Surface_Mesh"].link(meshAIM.output["Surface_Mesh"])

# Set number of cases to run simultaneously
cbaero.input.NumParallelCase = 3;

# Set number of threads to use to solve each case
cbaero.input.NumThreadPerCase = 2;

# Set AoA number 
cbaero.input.Alpha = [1.0, 3.0, 5.0] * deg

# Set Mach number 
cbaero.input.Mach = [2.0, 2.5, 4.0]

# Set Dynamic pressure 
cbaero.input.Dynamic_Pressure = 10 * psi

# Explicitly run analysis (optional)
cbaero.runAnalysis()

# Get all ouputs
for i in cbaero.output.keys():
    print(str(i) + " = " + str(cbaero.output[i].value))
    