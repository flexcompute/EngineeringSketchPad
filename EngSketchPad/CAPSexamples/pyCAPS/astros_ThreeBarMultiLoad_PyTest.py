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
parser = argparse.ArgumentParser(description = 'Astros Three Bar Multiple Load PyTest Example',
                                 prog = 'astros_ThreeBarMultiLoad_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "AstrosThreeBarMultiLoad")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaThreeBar.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
astrosAIM = myProblem.analysis.create( aim = "astrosAIM",
                                       name = "astros")

# Set project name so a mesh file is generated
projectName = "threebar_astros_Test"
astrosAIM.input.Proj_Name = projectName
astrosAIM.input.File_Format = "Free"
astrosAIM.input.Mesh_File_Format = "Large"
astrosAIM.input.Edge_Point_Max = 2
astrosAIM.input.Edge_Point_Min = 2
astrosAIM.input.Analysis_Type = "Static"

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
constraints = {}

constraint = {"groupName"         : ["boundary"],
              "dofConstraint"     : 123456}
constraints["conOne"] = constraint

constraint = {"groupName"         : ["boundary"],
              "dofConstraint"     : 123}
constraints["conTwo"] = constraint

astrosAIM.input.Constraint = constraints

loads = {}

load = {"groupName"         : "force",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 20000.0,
        "directionVector"   : [0.8, -0.6, 0.0]}
loads["loadOne"] = load

load = {"groupName"         : "force",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 20000.0,
        "directionVector"   : [-0.8, -0.6, 0.0]}
loads["loadTwo"] = load

astrosAIM.input.Load = loads

analysisCases = {}

value = {"analysisType"         : "Static",
         "analysisConstraint"   : "conOne",
         "analysisLoad"         : "loadOne"}
analysisCases["analysisOne"] = value

value = {"analysisType"         : "Static",
         "analysisConstraint"   : "conTwo",
         "analysisLoad"         : "loadTwo"}
analysisCases["analysisTwo"] = value

myProblem.analysis["astros"].input.Analysis = analysisCases

# Run AIM
astrosAIM.runAnalysis()
