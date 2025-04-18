# Import pyCAPS module
import pyCAPS

# Import os module
import os
import shutil
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

# Create project name
projectName = "NastranCompositeWing_Design"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Initialize CAPS Problem
geometryScript = os.path.join("..","csmData","compositeWing.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 40
egads.input.Edge_Point_Min = 5

# All quads in the grid
egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.05,.5,15]

# Load nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")

nastran.input["Mesh"].link(egads.output["Surface_Mesh"])
nastran.input.File_Format      = "Small"
nastran.input.Mesh_File_Format = "Large"

nastran.input.Analysis_Type    = "StaticOpt"

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

nastran.input.Material = {"Aluminum": Aluminum,
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

#nastran.input.Property = {"wing": aluminum}
nastran.input.Property = {"wing": composite}

constraint = {"groupName" : "root",
              "dofConstraint" : 123456}

nastran.input.Constraint = {"BoundaryCondition": constraint}

load = {"groupName" : "wing",
        "loadType" : "Pressure",
        "pressureForce" : 1.0}

nastran.input.Load = {"appliedPressure": load}

value = {"analysisType"         : "Static",
         "analysisConstraint"   : "BoundaryCondition",
         "analysisLoad"         : "appliedPressure"}

nastran.input.Analysis = {"StaticAnalysis": value}

designVariables = {}
designVarRelations = {}
numT = 8
for i in range(1, numT+1):
    dvName = "L{}".format(i)
    dv =  {"initialValue" : 0.00525,
            "lowerBound" : 0.00525*0.5,
            "upperBound" : 0.00525*1.5,
            "maxDelta"   : 0.00525*0.1}
    designVariables[dvName] = dv

    dvRelName = "L{}R".format(i)
    dvRel = {"componentType": "Property",
             "componentName": "wing",
             "variableName": dvName,
             "fieldName": "T{}".format(i)}  # T1, T2, ..., T8
    designVarRelations[dvRelName] = dvRel
              
capsProblem.analysis["nastran"].input.Design_Variable = designVariables
capsProblem.analysis["nastran"].input.Design_Variable_Relation = designVarRelations

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

capsProblem.analysis["nastran"].input.Design_Constraint = {"stress1": designConstraint1,
                                                           "stress2": designConstraint2,
                                                           "stress3": designConstraint3,
                                                           "stress4": designConstraint4,
                                                           "stress5": designConstraint5,
                                                           "stress6": designConstraint6,
                                                           "stress7": designConstraint7,
                                                           "stress8": designConstraint8}

# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......" )

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name +  ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")

# Run AIM post-analysis
nastran.postAnalysis()
