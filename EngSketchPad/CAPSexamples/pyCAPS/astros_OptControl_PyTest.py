# Import pyCAPS class file
from pyCAPS import capsProblem

# Import os module
import os
import sys
import shutil
import argparse

# ASTROS_ROOT should be the path containing ASTRO.D01 and ASTRO.IDX
try:
   ASTROS_ROOT = os.environ["ASTROS_ROOT"]
   os.putenv("PATH", ASTROS_ROOT + os.pathsep + os.getenv("PATH"))
except KeyError:
   print("Please set the environment variable ASTROS_ROOT")
   sys.exit(1)

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Astros Single Load Optimization Control Case PyTest Example',
                                 prog = 'astros_OptControl_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-verbosity", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Initialize capsProblem object
myProblem = capsProblem()

# Load CSM file
myProblem.loadCAPS("../csmData/feaSimplePlate.csm", verbosity=args.verbosity)

# Create project name
projectName = "AstrosSingleLoadPlateOptControl"

# Load Astros aim
myProblem.loadAIM( aim = "astrosAIM",
                   altName = "astros",
                   analysisDir = os.path.join(str(args.workDir[0]), projectName) )

# Set project name so a mesh file is generated
myProblem.analysis["astros"].setAnalysisVal("Proj_Name", projectName)

# Set meshing parameters
myProblem.analysis["astros"].setAnalysisVal("Edge_Point_Max", 4)

myProblem.analysis["astros"].setAnalysisVal("Quad_Mesh", True)

myProblem.analysis["astros"].setAnalysisVal("Tess_Params", [.25,.01,15])

# Set analysis type
myProblem.analysis["astros"].setAnalysisVal("Analysis_Type", "Static");

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}

myProblem.analysis["astros"].setAnalysisVal("Material", ("Madeupium", madeupium))

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["astros"].setAnalysisVal("Property", ("plate", shell))

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

myProblem.analysis["astros"].setAnalysisVal("Constraint", ("edgeConstraint", constraint))

# Set load
load = {"groupName" : "plate",
        "loadType" : "Pressure",
        "pressureForce" : 2.e6}

# Set loads
myProblem.analysis["astros"].setAnalysisVal("Load", ("appliedPressure", load ))

# Set Optimization Control
optControl = {"fullyStressedDesign" : 5}

myProblem.analysis["astros"].setAnalysisVal("OptimizationControl", optControl)

# Set analysis
# No analysis case information needs to be set for a single static load case

# Run AIM pre-analysis
myProblem.analysis["astros"].preAnalysis()

####### Run Astros ####################
print ("\n\nRunning Astros......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["astros"].analysisDir) # Move into test directory

# Copy files needed to run astros
astros_files = ["ASTRO.D01",  # *.DO1 file
                "ASTRO.IDX"]  # *.IDX file
for file in astros_files:
    if not os.path.isfile(file):
        shutil.copy(ASTROS_ROOT + os.sep + file, file)

if (args.noAnalysis == False):
    # Run Astros via system call
    os.system("astros.exe < " + projectName +  ".dat > " + projectName + ".out");


os.chdir(currentDirectory) # Move back to working directory
print ("Done running Astros!")
######################################

# Run AIM post-analysis
myProblem.analysis["astros"].postAnalysis()

# Close CAPS
myProblem.closeCAPS()
