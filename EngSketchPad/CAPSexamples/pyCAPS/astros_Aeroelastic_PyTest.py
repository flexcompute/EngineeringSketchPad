
# Import pyCAPS module
import pyCAPS

# Import os module
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
parser = argparse.ArgumentParser(description = 'Astros Aeroelastic PyTest Example',
                                 prog = 'astros_Aeroelastic_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir= os.path.join(str(args.workDir[0]), "AstrosAeroWingBEM")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaWingBEMAero.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load astros aim
myProblem.analysis.create(aim = "astrosAIM",
                          name = "astros")

# Set project name so a mesh file is generated
projectName = "astrosAero"
myProblem.analysis["astros"].input.Proj_Name = projectName

myProblem.analysis["astros"].input.Edge_Point_Max = 4

myProblem.analysis["astros"].input.Quad_Mesh = True

# Set analysis type
myProblem.analysis["astros"].input.Analysis_Type = "Aeroelastic"

# Set analysis
trim1 = { "analysisType" : "AeroelasticStatic",
          "trimSymmetry" : "SYM",
          "analysisConstraint" : ["ribConstraint"],
          "analysisSupport" : ["ribSupport"],
          "machNumber"     : 0.5,
          "dynamicPressure": 50000,
          "density" : 1.0,
          "rigidVariable"  : "ANGLEA",
          "rigidConstraint": ["URDD3"],
          "magRigidConstraint" : [-50],
          }


myProblem.analysis["astros"].input.Analysis = {"Trim1": trim1}

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

myProblem.analysis["astros"].input.Material = {"Unobtainium": unobtainium}

# Set property
shell  = {"propertyType"        : "Shell",
          "membraneThickness"   : 0.2,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

shell2  = {"propertyType"       : "Shell",
          "membraneThickness"   : 0.002,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

myProblem.analysis["astros"].input.Property = {"Ribs"    : shell,
                                               "Spar1"   : shell,
                                               "Spar2"   : shell,
                                               "Rib_Root": shell,
                                               "Skin"    : shell2}


# Defined Connections
connection = {"dofDependent"   : 123456,
              "connectionType" : "RigidBody"}

myProblem.analysis["astros"].input.Connect = {"Rib_Root": connection}


# Set constraints
constraint = {"groupName"     : ["Rib_Root_Point"],
              "dofConstraint" : 12456}

myProblem.analysis["astros"].input.Constraint = {"ribConstraint": constraint}

# Set supports
support = {"groupName" : ["Rib_Root_Point"],
           "dofSupport": 3}

myProblem.analysis["astros"].input.Support = {"ribSupport": support}


# Aero
wing = {"groupName"    : "Wing",
        "numChord"     : 8,
        "numSpanTotal" : 24}

# Note the surface name corresponds to the capsBound found in the *.csm file. This links
# the spline for the aerodynamic surface to the structural model
myProblem.analysis["astros"].input.VLM_Surface = {"Skin_Top": wing}

# astros is executed automatically just-in-time

# Run AIM 
myProblem.analysis["astros"].runAnalysis()
