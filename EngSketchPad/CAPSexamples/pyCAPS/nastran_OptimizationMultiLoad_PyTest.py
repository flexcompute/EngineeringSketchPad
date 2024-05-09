from __future__ import print_function

# Import pyCAPS class file
from pyCAPS import capsProblem

# Import os module
import os
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

# Initialize capsProblem object
myProblem = capsProblem()

# Load CSM file
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem.loadCAPS(geometryScript, verbosity=args.outLevel)

# Create project name
projectName = "NastranOptMultiLoadPlate"

# Load Nastran aim
myMesh = myProblem.loadAIM(aim = "egadsTessAIM",
                           analysisDir = os.path.join(str(args.workDir[0]), projectName))

# Set project name so a mesh file is generated
myMesh.setAnalysisVal("Proj_Name", projectName)

# Set meshing parameters
myMesh.setAnalysisVal("Edge_Point_Max", 5)
myMesh.setAnalysisVal("Edge_Point_Min", 5)

myMesh.setAnalysisVal("Mesh_Elements", "Quad")

myMesh.setAnalysisVal("Tess_Params", [.25,.01,15])

# Generate mesh
myMesh.preAnalysis()
myMesh.postAnalysis()

# Load Nastran aim
myAnalysis = myProblem.loadAIM(aim = "nastranAIM",
                               parents = myMesh.aimName,
                               analysisDir = os.path.join(str(args.workDir[0]), projectName))

myAnalysis.setAnalysisVal("File_Format", "Free")

myAnalysis.setAnalysisVal("Mesh_File_Format", "Large")

# Set analysis type
myAnalysis.setAnalysisVal("Analysis_Type", "StaticOpt");

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3,
                "yieldAllow" : 1E7}


myAnalysis.setAnalysisVal("Material", ("Madeupium", madeupium))

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

myAnalysis.setAnalysisVal("Property", ("plate", shell))

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

myAnalysis.setAnalysisVal("Constraint", ("edgeConstraint", constraint))

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
myAnalysis.setAnalysisVal("Load", loads)

# Set design variables and analysis - design variable relation
materialDV = {"groupName" : "madeupium",
              "designVariableType" : "Material",
              "initialValue"  : madeupium["youngModulus"],
              "discreteValue" : [madeupium["youngModulus"]*.5,
                                 madeupium["youngModulus"]*1.0,
                                 madeupium["youngModulus"]*1.5],
              "fieldName" : "E"}

propertyDV = {"groupName" : "plate",
              "initialValue" : shell["membraneThickness"],
              "lowerBound" : shell["membraneThickness"]*0.5,
              "upperBound" : shell["membraneThickness"]*1.5,
              "maxDelta"   : shell["membraneThickness"]*0.1,
              "fieldName" : "T"}

# Set design variable
myAnalysis.setAnalysisVal("Design_Variable", [("youngModulus", materialDV),
                                                                 ("membranethickness", propertyDV)])

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
myAnalysis.setAnalysisVal("Design_Constraint",[("stress1", designConstraint1),
                                                                  ("stress2", designConstraint2)])

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
myAnalysis.setAnalysisVal("Analysis", analysisCases)

myAnalysis.setAnalysisVal("ObjectiveMinMax", "Min")
myAnalysis.setAnalysisVal("ObjectiveResponseType",  "WEIGHT")


# Run AIM pre-analysis
myAnalysis.preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myAnalysis.analysisDir) # Move into test directory

os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
#os.system("nastran.sh " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running NASTRAN!")
######################################

# Run AIM post-analysis
myAnalysis.postAnalysis()

# Close CAPS
myProblem.closeCAPS()
