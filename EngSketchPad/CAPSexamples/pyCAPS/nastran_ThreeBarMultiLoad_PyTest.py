## [importPrint]
from __future__ import print_function
## [importPrint]

## [import]
# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import argparse
## [import]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Three Bar Pytest Example',
                                 prog = 'nastran_ThreeBar_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()


workDir = os.path.join(str(args.workDir[0]), "NastranThreeBarMultiLoad")

## [initateProblem]
# Initialize CAPS Problem
geometryScript = os.path.join("..","csmData","feaThreeBar.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [initateProblem]

## [loadAIM]
# Load nastran aim
nastranAIM = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran",
                                       autoExec = False)
## [loadAIM]

## [setInputs]
# Set project name so a mesh file is generated
projectName = "threebar_nastran_Test"
nastranAIM.input.Proj_Name = projectName
nastranAIM.input.File_Format = "Free"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Edge_Point_Max = 2
nastranAIM.input.Edge_Point_Min = 2
nastranAIM.input.Analysis_Type = "Static"
## [setInputs]

## [defineMaterials]
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1}

nastranAIM.input.Material = {"Madeupium": madeupium}
## [defineMaterials]

## [defineProperties]
rod  =   {"propertyType"      : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 1.0}

rod2  =   {"propertyType"     : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 2.0}

nastranAIM.input.Property = {"bar1": rod,
                             "bar2": rod2,
                             "bar3": rod}
## [defineProperties]

# Set constraints
## [defineConstraints]
constraints = {}

conOne = {"groupName"         : ["boundary"],
          "dofConstraint"     : 123456}
conTwo = {"groupName"         : ["boundary"],
          "dofConstraint"     : 123}

nastranAIM.input.Constraint = {"conOne": conOne, "conTwo": conTwo}
## [defineConstraints]

## [defineLoad]
loadOne = {"groupName"         : "force",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 20000.0,
        "directionVector"   : [0.8, -0.6, 0.0]}

loadTwo = {"groupName"         : "force",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 20000.0,
        "directionVector"   : [-0.8, -0.6, 0.0]}

nastranAIM.input.Load = {"loadOne": loadOne, "loadTwo": loadTwo}

## [defineLoad]

## [defineAnalysis]
caseOne = {"analysisType"         : "Static",
         "analysisConstraint"     : "conOne",
         "analysisLoad"         : "loadOne"}

caseTwo = {"analysisType"         : "Static",
         "analysisConstraint"     : "conTwo",
         "analysisLoad"         : "loadTwo"}

nastranAIM.input.Analysis = {"caseOne": caseOne, "caseTwo": caseTwo}
## [defineAnalysis]

# Run AIM pre-analysis
## [preAnalysis]
nastranAIM.preAnalysis()
## [preAnalysis]

## [run]
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory
os.chdir(nastranAIM.analysisDir) # Move into test directory
if (args.noAnalysis == False): os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")
## [run]

## [postAnalysis]
nastranAIM.postAnalysis()
## [postAnalysis]
