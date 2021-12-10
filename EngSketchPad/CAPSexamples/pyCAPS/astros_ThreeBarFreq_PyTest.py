# Import pyCAPS module
import pyCAPS

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
parser = argparse.ArgumentParser(description = 'Astros Three Bar Frequency PyTest Example',
                                 prog = 'astros_ThreeBarFreq_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "AstrosThreeBarFreq")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaThreeBar.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
astrosAIM = myProblem.analysis.create(aim = "astrosAIM",
                                      name = "astros")

# Set project name so a mesh file is generated
projectName = "threebar_astros_Test"
astrosAIM.input.Proj_Name = projectName
astrosAIM.input.Mesh_File_Format = "Large"
astrosAIM.input.Edge_Point_Max = 2
astrosAIM.input.Edge_Point_Min = 2
astrosAIM.input.Analysis_Type = "Modal"
## [setInputs]

## [defineMaterials]
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1}

astrosAIM.input.Material = {"Madeupium": madeupium}

rod  =   {"propertyType"      : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 1.0}

rod2  =   {"propertyType"     : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 2.0}

astrosAIM.input.Property = {"bar1": rod,
                            "bar2": rod2,
                            "bar3": rod}

# Set constraints
constraint = {"groupName"         : ["boundary"],
              "dofConstraint"     : 123456}

astrosAIM.input.Constraint = {"BoundaryCondition": constraint}

eigen = { "extractionMethod"     : "MGIV",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormaliztion"    : "MASS"}

astrosAIM.input.Analysis = {"EigenAnalysis": eigen}

# astros is executed automatically just-in-time

# Get Eigen-frequencies
print ("\nGetting results for natural frequencies.....")
natrualFreq = myProblem.analysis["astros"].output.EigenFrequency

mode = 1
for i in natrualFreq:
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
    mode += 1
