# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Aeroelastic Pytest Example',
                                 prog = 'nastran_Aeroelastic_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "NastranOptMultiLoadPlate"

workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)


# Load nastran aim
myProblem.analysis.create(aim = "nastranAIM",
                          name = "nastran")

# Set project name so a mesh file is generated
myProblem.analysis["nastran"].input.Proj_Name = projectName

# Set meshing parameters
myProblem.analysis["nastran"].input.Edge_Point_Max = 5
myProblem.analysis["nastran"].input.Edge_Point_Min = 5

myProblem.analysis["nastran"].input.Quad_Mesh = True

myProblem.analysis["nastran"].input.Tess_Params = [.25,.01,15]

myProblem.analysis["nastran"].input.File_Format = "Free"

myProblem.analysis["nastran"].input.Mesh_File_Format = "Large"

# Set analysis type
myProblem.analysis["nastran"].input.Analysis_Type = "StaticOpt"

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3,
                "yieldAllow" : 1E7}


myProblem.analysis["nastran"].input.Material = {"Madeupium": madeupium}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["nastran"].input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

myProblem.analysis["nastran"].input.Constraint = {"edgeConstraint": constraint}

# Create multiple loads
numLoad = 2
loads = {}
for i in range(numLoad):
    # Set load name
    name = "load_"+ str(i)

    # Set load values
    value =  {"groupName" : "plate",
              "loadType" : "Pressure",
              "pressureForce" : 1.e6*(i+1)}

    # Add loadElement to laod
    loads[name] = value

# Set loads
myProblem.analysis["nastran"].input.Load = loads

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
myProblem.analysis["nastran"].input.Design_Variable = {"youngModulus": materialDV,
                                                       "membranethickness": propertyDV}

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
myProblem.analysis["nastran"].input.Design_Constraint = {"stress1": designConstraint1,
                                                         "stress2": designConstraint2}

# Create multiple analysis cases
analysisCases = {}
for name in loads.keys():

    # set analysis value
    value = {"analysisType" : "Static",
             "analysisLoad" : name,
             "analysisDesignConstraint" : "stress" + str(i+1),
             "objectiveMinMax" : "MIN",
             "objectiveResponseType" : "WEIGHT"}

    # Add analysisElement to analysis cases
    analysisCases[name] = value

# Set analysis
myProblem.analysis["nastran"].input.Analysis = analysisCases

# Run AIM pre-analysis
myProblem.analysis["nastran"].preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["nastran"].analysisDir) # Move into test directory

os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
#os.system("nastran.sh " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running NASTRAN!")
######################################

# Run AIM post-analysis
myProblem.analysis["nastran"].postAnalysis()
