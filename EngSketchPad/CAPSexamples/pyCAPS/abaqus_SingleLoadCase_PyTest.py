# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse
import time

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'ABAQUS Pytest Example',
                                 prog = 'abaqus_SingleLoadCase_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create project name
projectName = "AbaqusSingleLoadPlate"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Meshing
myMesh = myProblem.analysis.create(aim = "egadsTessAIM",
                                   name = "egadsTess" )

# Set meshing parameters
myMesh.input.Edge_Point_Max = 10
myMesh.input.Edge_Point_Min = 6

myMesh.input.Mesh_Elements = "Quad"

myMesh.input.Tess_Params = [.25,.01,15]

# Load abaqus aim
myProblem.analysis.create(aim = "abaqusAIM",
                          name = "abaqus")

# Set mesh
myProblem.analysis["abaqus"].input["Mesh"].link(myMesh.output["Surface_Mesh"])

# Set analysis type
myProblem.analysis["abaqus"].input.Analysis_Type = "Static"

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}

myProblem.analysis["abaqus"].input.Material = {"Madeupium": madeupium}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["abaqus"].input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

myProblem.analysis["abaqus"].input.Constraint = {"edgeConstraint": constraint}

# Set load
load = {"groupName" : "plate",
        "loadType" : "Pressure",
        "pressureForce" : 2.e6}

# Set loads
myProblem.analysis["abaqus"].input.Load = {"appliedPressure": load}

# Set analysis
# No analysis case information needs to be set for a single static load case

# Run AIM pre-analysis
myProblem.analysis["abaqus"].preAnalysis()

####### Run Abaqus ####################
print ("\n\nRunning Abaqus......")
# Use 'interactive' to force abaqus to finish writing files
# unset LD_PRELOAD for testing with LLVM address sanitizer
myProblem.analysis["abaqus"].system("unset LD_PRELOAD ; abaqus interactive job=" + myProblem.analysis["abaqus"].input.Proj_Name); # Run Abaqus via system call
print ("Done running Abaqus!")
########################################

# Run AIM post-analysis
myProblem.analysis["abaqus"].postAnalysis()

print("MISES = ", myProblem.analysis["abaqus"].dynout.vonMises_Grid)
print("Displacement = ", myProblem.analysis["abaqus"].dynout.Displacement)

