# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Composite Wing Design Pytest Example',
                                 prog = 'nastran_CompositeWingDesign_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "NastranCompositeWing")

# Load CSM file
geometryScript = os.path.join("..","csmData","compositeWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load nastran aim
nastranAIM = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")

# Set project name so a mesh file is generated
projectName = "nastran_Composite_Wing"
nastranAIM.input.Proj_Name        = projectName
nastranAIM.input.Edge_Point_Max   = 40
nastranAIM.input.File_Format      = "Small"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Analysis_Type    = "StaticOpt"

Aluminum  = {"youngModulus" : 10.5E6 ,
             "poissonRatio" : 0.3,
             "density"      : 0.1/386,
             "shearModulus" : 4.04E6}

Graphite_epoxy = {"materialType"        : "Orthotropic",
                  "youngModulus"        : 20.8E6 ,
                  "youngModulusLateral" : 1.54E6,
                  "poissonRatio"        : 0.327,
                  "shearModulus"        : 0.80E6,
                  "density"             : 0.059/386,
                  "tensionAllow"        : 11.2e-3,
                  "tensionAllowLateral" : 4.7e-3,
                  "compressAllow"       : 11.2e-3,
                  "compressAllowLateral": 4.7e-3,
                  "shearAllow"          : 19.0e-3,
                  "allowType"           : 1}

nastranAIM.input.Material = {"Aluminum": Aluminum,
                             "Graphite_epoxy": Graphite_epoxy}

aluminum  = {"propertyType"        : "Shell",
             "material"            : "Aluminum",
             "bendingInertiaRatio" : 1.0, # Default - not necesssary
             "shearMembraneRatio"  : 0, # Turn of shear - no materialShear
             "membraneThickness"   : 0.125 }

composite  = {"propertyType"          : "Composite",
              "shearBondAllowable"      : 1.0e6,
              "bendingInertiaRatio"   : 1.0, # Default - not necesssary
              "shearMembraneRatio"    : 0, # Turn of shear - no materialShear
              "compositeMaterial"     : ["Graphite_epoxy"]*8,
              "compositeThickness"    : [0.00525]*8,
              "compositeOrientation"  : [0, 0, 0, 0, -45, 45, -45, 45],
              "symmetricLaminate"     : 1,
              "compositeFailureTheory": "STRN" }

#nastranAIM.input.Property = {"wing": aluminum}
nastranAIM.input.Property = {"wing": composite}

constraint = {"groupName" : "root",
              "dofConstraint" : 123456}

nastranAIM.input.Constraint = {"BoundaryCondition": constraint}

load = {"groupName" : "wing",
        "loadType" : "Pressure",
        "pressureForce" : 1.0}

nastranAIM.input.Load = {"appliedPressure": load}

value = {"analysisType"         : "Static",
         "analysisConstraint"     : "BoundaryCondition",
         "analysisLoad"         : "appliedPressure"}

nastranAIM.input.Analysis = {"StaticAnalysis": value}

DesVar1    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T1"}

DesVar2    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T2"}

DesVar3    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T3"}

DesVar4    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T4"}

DesVar5    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T5"}

DesVar6    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T6"}

DesVar7    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T7"}

DesVar8    = {"groupName" : "wing",
              "initialValue" : 0.00525,
              "lowerBound" : 0.00525*0.5,
              "upperBound" : 0.00525*1.5,
              "maxDelta"   : 0.00525*0.1,
              "fieldName" : "T8"}


myProblem.analysis["nastran"].input.Design_Variable = {"L1": DesVar1,
                                                       "L2": DesVar2,
                                                       "L3": DesVar3,
                                                       "L4": DesVar4,
                                                       "L5": DesVar5,
                                                       "L6": DesVar6,
                                                       "L7": DesVar7,
                                                       "L8": DesVar8}

designConstraint1 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA1"}

designConstraint2 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA2"}

designConstraint3 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA3"}

designConstraint4 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA4"}

designConstraint5 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA5"}

designConstraint6 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA6"}

designConstraint7 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA7"}

designConstraint8 = {"groupName" : "wing",
                    "responseType" : "CFAILURE",
                    "lowerBound" : 0.0,
                    "upperBound" :  0.9999,
                    "fieldName" : "LAMINA8"}

myProblem.analysis["nastran"].input.Design_Constraint = {"stress1": designConstraint1,
                                                         "stress2": designConstraint2,
                                                         "stress3": designConstraint3,
                                                         "stress4": designConstraint4,
                                                         "stress5": designConstraint5,
                                                         "stress6": designConstraint6,
                                                         "stress7": designConstraint7,
                                                         "stress8": designConstraint8}

# Run AIM pre-analysis
nastranAIM.preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......" )
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["nastran"].analysisDir) # Move into test directory

if (args.noAnalysis == False):
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

# Run AIM post-analysis
nastranAIM.postAnalysis()
