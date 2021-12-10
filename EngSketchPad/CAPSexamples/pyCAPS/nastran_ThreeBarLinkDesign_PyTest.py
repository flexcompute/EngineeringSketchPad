# Import pyCAPS module
import pyCAPS

# Import os module
import argparse
import os

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Three Bar Link Design Pytest Example',
                                 prog = 'nastran_ThreeBarLinkDesign_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "NastranThreeBarLinkDesign")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaThreeBar.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
nastranAIM = myProblem.analysis.create( aim = "nastranAIM",
                                        name = "nastran")

# Set project name so a mesh file is generated
projectName = "threebar_nastran_Test"
nastranAIM.input.Proj_Name = projectName
nastranAIM.input.File_Format = "Free"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Edge_Point_Max = 2
nastranAIM.input.Edge_Point_Min = 2
nastranAIM.input.Analysis_Type = "StaticOpt"


madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.0E7 ,
                "poissonRatio" : .33,
                "density"      : 0.1,
                "yieldAllow"   : 5.6E7}

nastranAIM.input.Material = {"Madeupium": madeupium}


rod  =   {"propertyType"      : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 1.0}

rod2  =   {"propertyType"      : "Rod",
          "material"          : "Madeupium",
          "crossSecArea"      : 2.0}

nastranAIM.input.Property = {"bar1": rod,
                             "bar2": rod2,
                             "bar3": rod}

DesVar1    = {"groupName" : "bar1",
              "initialValue" : rod["crossSecArea"],
              "lowerBound" : rod["crossSecArea"]*0.5,
              "upperBound" : rod["crossSecArea"]*1.5,
              "maxDelta"   : rod["crossSecArea"]*0.1,
              "fieldName" : "A"}

DesVar2    = {"groupName" : "bar2",
              "initialValue" : rod2["crossSecArea"],
              "lowerBound" : rod2["crossSecArea"]*0.5,
              "upperBound" : rod2["crossSecArea"]*1.5,
              "maxDelta"   : rod2["crossSecArea"]*0.1,
              "fieldName" : "A"}

DesVar3    = {"groupName" : "bar3",
              "initialValue" : rod["crossSecArea"],
              "lowerBound" : rod["crossSecArea"]*0.5,
              "upperBound" : rod["crossSecArea"]*1.5,
              "maxDelta"   : rod["crossSecArea"]*0.1,
              "variableWeight" : [0.0, 1.0],
              "independentVariableWeight" : 1.0,
              "independentVariable" : "Bar1A",
              "fieldName" : "A"}

myProblem.analysis["nastran"].input.Design_Variable = {"Bar1A": DesVar1,
                                                       "Bar2A": DesVar2,
                                                       "Bar3A": DesVar3}

# Set constraints
constraint = {"groupName" : ["boundary"],
              "dofConstraint" : 123456}

nastranAIM.input.Constraint = {"BoundaryCondition": constraint}


designConstraint1 = {"groupName" : "bar1",
                    "responseType" : "STRESS",
                    "lowerBound" : -madeupium["yieldAllow"],
                    "upperBound" :  madeupium["yieldAllow"]}

designConstraint2 = {"groupName" : "bar2",
                    "responseType" : "STRESS",
                    "lowerBound" : -madeupium["yieldAllow"],
                    "upperBound" :  madeupium["yieldAllow"]}

designConstraint3 = {"groupName" : "bar3",
                    "responseType" : "STRESS",
                    "lowerBound" : -madeupium["yieldAllow"],
                    "upperBound" :  madeupium["yieldAllow"]}

myProblem.analysis["nastran"].input.Design_Constraint = {"stress1": designConstraint1,
                                                         "stress2": designConstraint2,
                                                         "stress3": designConstraint3}

load = {"groupName" : "force",
        "loadType" : "GridForce",
        "forceScaleFactor" : 20000.0,
        "directionVector" : [0.8, -0.6, 0.0]}

nastranAIM.input.Load = {"appliedForce": load}


value = {"analysisType"         : "Static",
         "analysisConstraint"     : "BoundaryCondition",
         "analysisLoad"         : "appliedForce"}


myProblem.analysis["nastran"].input.Analysis = {"SingleLoadCase": value}

# Run AIM pre-analysis
nastranAIM.preAnalysis()

print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory
os.chdir(nastranAIM.analysisDir) # Move into test directory
if (args.noAnalysis == False): os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

nastranAIM.postAnalysis()
