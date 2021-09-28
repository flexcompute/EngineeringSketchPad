# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'HSM Pytest Example',
                                 prog = 'hsm_SingleLoadCase_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "HSMSingleLoadPlate"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load mystran aim
myAnalysis = myProblem.analysis.create(aim = "hsmAIM")

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = projectName

# Set meshing parameters
myAnalysis.input.Edge_Point_Max = 10
myAnalysis.input.Edge_Point_Min = 9

myAnalysis.input.Quad_Mesh = True

myAnalysis.input.Tess_Params = [.25,.01,15]

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
          "shearMembraneRatio"  : 5.0/6.0} # Default

myAnalysis.input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName"     : "plateEdge",
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"edgeConstraint": constraint}

# Set load
load = {"groupName" : "plate",
        "loadType" : "Pressure",
        "pressureForce" : 2.e6}

# Set loads
myAnalysis.input.Load = {"appliedPressure": load}

# Run AIM
myAnalysis.runAnalysis()
