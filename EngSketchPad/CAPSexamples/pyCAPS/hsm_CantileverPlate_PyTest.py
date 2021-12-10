# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse
import math

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'HSM Pytest Example',
                                 prog = 'hsm_CantileverPlate_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "HSMCantileverPlate"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaCantileverPlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

plateLength = 1.0

# Modify the size of the plate
myProblem.geometry.despmtr.plateLength = plateLength
myProblem.geometry.despmtr.plateWidth  = 0.1

# Load mystran aim
myAnalysis = myProblem.analysis.create(aim = "hsmAIM")

# Set project name so a mesh file is generated
myAnalysis.input.Proj_Name = projectName

# Set meshing parameters
myAnalysis.input.Edge_Point_Max = 10
myAnalysis.input.Edge_Point_Min = 2

myAnalysis.input.Quad_Mesh = True

#myAnalysis.input.Tess_Params", [.25,.01,15])

youngModulus = 10000.0
poissonRatio = 0
shearModulus = youngModulus/(2*(1+poissonRatio))

# Set materials
madeupium    = {"materialType" : "isotropic",
                "massPerArea"  : 1.0,
                "youngModulus" : youngModulus ,
                "poissonRatio" : poissonRatio,
                "shearModulus" : shearModulus}

myAnalysis.input.Material = {"Madeupium": madeupium}

tshell = 0.1 * pow(1.2, 1.0/3.0);

# Set properties
shell  = {"propertyType"      : "Shell",
          "membraneThickness" : tshell,
          "material"          : "madeupium"}

myAnalysis.input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName"     : "plateEdge",
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"edgeConstraint": constraint}

# D = t^3 E / [ 12 (1-v^2) ]
dpar = tshell**3 * youngModulus / (12.0 * (1.0-poissonRatio*poissonRatio));

# end bending moment to bend beam into a circle segment
moment = -2.0*math.atan(1.)*4.*dpar/plateLength;

# Set load
load = {"groupName"         : "plateTip",
        "loadType"          : "LineMoment",
        "momentScaleFactor" : moment,
        "directionVector"   : [0.0, 1.0, 0.0]}

# Set loads
myAnalysis.input.Load = {"appliedLoad": load}

# Run AIM
myAnalysis.runAnalysis()
