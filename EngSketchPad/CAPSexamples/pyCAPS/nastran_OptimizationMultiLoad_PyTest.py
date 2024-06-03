# Import pyCAPS
import pyCAPS

# Import os module
import os
import shutil
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Multiload Opt Pytest Example',
                                 prog = 'nastran_OptimizationMultiLoad_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "NastranOptMultiLoadPlate"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Initialize capsProblem object
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)


# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 5
egads.input.Edge_Point_Min = 5

egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.25,.01,15]

# Load Nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM")

nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.File_Format = "Free"

nastran.input.Mesh_File_Format = "Large"

# Set analysis type
nastran.input.Analysis_Type = "StaticOpt";

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3,
                "yieldAllow" : 1E7}


nastran.input.Material = {"Madeupium": madeupium}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

nastran.input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

nastran.input.Constraint = {"edgeConstraint": constraint}

# Create multiple loads
numLoad = 2
loads = []
for i in range(numLoad):
    # Set load name
    name = "load_"+ str(i)

    # Set load values
    value =  {"groupName" : "plate",
              "loadType" : "Pressure",
              "pressureForce" : 1.e6*(i+1)}

    # Create temporary load tuple
    loadElement = (name, value)

    # Append loadElement to laod
    loads.append( loadElement)

# Set loads
nastran.input.Load = dict(loads)

# Set design variables and analysis - design variable relation
materialDV = {"initialValue"  : madeupium["youngModulus"],
              "discreteValue" : [madeupium["youngModulus"]*.5,
                                 madeupium["youngModulus"]*1.0,
                                 madeupium["youngModulus"]*1.5]}
materialDVR = {"componentType": "Material",
               "componentName": "madeupium",
               "variableName": "propertyDV",
               "fieldName" : "E"}

propertyDV = {"initialValue" : shell["membraneThickness"],
              "lowerBound" : shell["membraneThickness"]*0.5,
              "upperBound" : shell["membraneThickness"]*1.5,
              "maxDelta"   : shell["membraneThickness"]*0.1}
propertyDVR = {"componentType": "Property",
               "componentName": "plate",
               "variableName": "propertyDV",
               "fieldName" : "T"}


# Set design variable
nastran.input.Design_Variable = {"materialDV": materialDV,
                                 "propertyDV": propertyDV}
nastran.input.Design_Variable_Relation = {"materialDVR": materialDVR,
                                          "propertyDVR": propertyDVR}

# Set design constraints and responses

designConstraint1 = {"groupName" : "plate",
                    "responseType" : "STRESS",
                    "lowerBound" : -madeupium["yieldAllow"],
                    "upperBound" :  madeupium["yieldAllow"],}

designConstraint2 = {"groupName" : "plate",
                    "responseType" : "STRESS",
                    "lowerBound" : -0.5*madeupium["yieldAllow"],
                    "upperBound" :  0.5*madeupium["yieldAllow"],}


# Set design constraint
nastran.input.Design_Constraint = {"stress1": designConstraint1,
                                   "stress2": designConstraint2}

# Create multiple analysis cases
analysisCases = []
for i in range(numLoad):

    # Set analysis name
    name = loads[i][0]

    # set analysis value
    value = {"analysisType" : "Static",
             "analysisLoad" : name,
             "analysisDesignConstraint" : "stress" + str(i+1),
             }

    # Create temporary analysis tuple
    analysisElement = (name, value)

    # Append analysisElement to analysis cases
    analysisCases.append(analysisElement)

# Set analysis
nastran.input.Analysis = dict(analysisCases)

nastran.input.Objective_Min_Max = "Min"
nastran.input.Objective_Response_Type =  "WEIGHT"


# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name + ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running NASTRAN!")
######################################

# Run AIM post-analysis
nastran.postAnalysis()
