# Import pyCAPS module
import pyCAPS

# Import modules
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
parser = argparse.ArgumentParser(description = 'Astros Multiple Load Pytest Example',
                                 prog = 'astros_MultipleLoadCase_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "AstrosMultipleLoadPlate"

workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)


# Load Astros aim
myAnalysis = myProblem.analysis.create( aim = "astrosAIM",
                                        name = "astros" )

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = projectName

# Set meshing parameters
myAnalysis.input.Edge_Point_Max = 4
myAnalysis.input.Edge_Point_Min = 4

myAnalysis.input.Quad_Mesh = True

myAnalysis.input.Tess_Params = [.25,.01,15]

# Set analysis type
myAnalysis.input.Analysis_Type = "Static"

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}

myAnalysis.input.Material = {"Madeupium": madeupium}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

myAnalysis.input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"edgeConstraint": constraint}

# Create multiple loads
numLoad = 5
loads = {}
for i in range(numLoad):
    # Set load name
    name = "load_"+ str(i)

    # Set load values
    value =  {"groupName" : "plate",
              "loadType" : "Pressure",
              "pressureForce" : 1.e6*(i+1)}

    # Add load to dict
    loads[name] = value

# Set loads
myAnalysis.input.Load = loads

# Create multiple analysis cases
analysisCases = {}
for name in loads.keys():

    # set analysis value s
    value = {"analysisType" : "Static",
             "analysisLoad" : name}

    # Add to analysis cases
    analysisCases[name] = value

# Set analysis
myAnalysis.input.Analysis = analysisCases

# Run AIM 
myAnalysis.runAnalysis()
