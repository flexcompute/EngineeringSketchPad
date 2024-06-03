## [import]
# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import shutil
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

# Create project name
projectName = "NastranThreeBar"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

## [initateProblem]
# Initialize CAPS Problem
geometryScript = os.path.join("..","csmData","feaThreeBar.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [initateProblem]

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 2
egads.input.Edge_Point_Min = 2

egads.input.Tess_Params = [.05,.5,15]


## [loadAIM]
# Load nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")
## [loadAIM]

## [setInputs]
nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.Proj_Name = "threebar_nastran_Test"
nastran.input.File_Format = "Free"
nastran.input.Mesh_File_Format = "Large"
nastran.input.Analysis_Type = "Static"
## [setInputs]

## [defineMaterials]
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1}

nastran.input.Material = {"Madeupium": madeupium}
## [defineMaterials]

## [defineProperties]
rod  =   {"propertyType"      : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 1.0}

rod2  =   {"propertyType"     : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 2.0}

nastran.input.Property = {"bar1": rod,
                             "bar2": rod2,
                             "bar3": rod}
## [defineProperties]

# Set constraints
## [defineConstraints]
constraint = {"groupName"         : ["boundary"],
              "dofConstraint"     : 123456}

nastran.input.Constraint = {"BoundaryCondition": constraint}
## [defineConstraints]

## [defineLoad]
load = {"groupName"         : "force",
        "loadType"          : "GridForce",
        "forceScaleFactor"  : 20000.0,
        "directionVector"   : [0.8, -0.6, 0.0]}

nastran.input.Load = {"appliedForce": load }
## [defineLoad]

## [defineAnalysis]
value = {"analysisType"         : "Static",
         "analysisConstraint"   : "BoundaryCondition",
         "analysisLoad"         : "appliedForce"}

capsProblem.analysis["nastran"].input.Analysis = {"SingleLoadCase": value }
## [defineAnalysis]

# Run AIM pre-analysis
## [preAnalysis]
nastran.preAnalysis()
## [preAnalysis]

## [run]
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name + ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")
## [run]

## [postAnalysis]
nastran.postAnalysis()
## [postAnalysis]
