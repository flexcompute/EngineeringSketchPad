# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import argparse
import time

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Abaqus AGARD445.6 Pytest Example',
                                 prog = 'abaqus_AGARD445_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()


# Create project name
projectName = "AbaqusModalAGARD445"
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaAGARD445.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change the sweepAngle and span of the Geometry - Demo purposes
#myProblem.geometry.despmtr.sweepAngle = 5 # From 45 to 5 degrees
#myProblem.geometry.despmtr.semiSpan   = 5 # From 2.5 ft to 5 ft

# Meshing
myMesh = myProblem.analysis.create(aim = "egadsTessAIM",
                                   name = "egadsTess" )

# Set meshing parameters
myMesh.input.Edge_Point_Max = 10
myMesh.input.Edge_Point_Min = 6

myMesh.input.Mesh_Elements = "Quad"

myMesh.input.Tess_Params = [.25,.01,15]

# Load abaqus aim
myAnalysis = myProblem.analysis.create(aim = "abaqusAIM",
                                       name = "abaqus")

# Set mesh
myAnalysis.input["Mesh"].link(myMesh.output["Surface_Mesh"])

# Set analysis type
myAnalysis.input.Analysis_Type = "Modal"

# Set analysis inputs
eigen = { "extractionMethod"     : "Lanczos",
          "frequencyRange"       : [0.1, 200],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 2,
          "eigenNormalization"   : "MASS",
          "lanczosMode"          : 2,  # Default - not necesssary
          "lanczosType"          : "DPB"} # Default - not necesssary

myAnalysis.input.Analysis = {"EigenAnalysis": eigen,
                             "StaticCase": {"analysisType":"Static", "analysisLoad":["gravity", "force"]}}

# Set materials
mahogany    = {"materialType"        : "orthotropic",
               "youngModulus"        : 0.457E6 ,
               "youngModulusLateral" : 0.0636E6,
               "poissonRatio"        : 0.31,
               "shearModulus"        : 0.0637E6,
               "shearModulusTrans1Z" : 0.00227E6,
               "shearModulusTrans2Z" : 0.00227E6,
               "density"             : 3.5742E-5}

myAnalysis.input.Material = {"Mahogany": mahogany}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.82,
          "material"        : "mahogany",
          "bendingInertiaRatio" : 1.0, # Default - not necesssary
          "shearMembraneRatio"  : 5.0/6.0} # Default - not necesssary

myAnalysis.input.Property = {"yatesPlate": shell}

# Set constraints
constraint = {"groupName" : "constEdge",
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"edgeConstraint": constraint}

myAnalysis.input.Load = {"gravity": {"groupName":"force", "loadType": "Gravity", "directionVector" :[-10, 9, 0],  "gravityAcceleration": 5.1},
                         "force":   {"loadType": "Pressure", "pressureForce" :9}}

# Run AIM pre-analysis
myAnalysis.preAnalysis()

####### Run Abaqus ####################
print ("\n\nRunning Abaqus......")
myAnalysis.system("abaqus job=" + myAnalysis.input.Proj_Name); # Run Abaqus via system call

print ("Done running Abaqus!")
########################################

# Run AIM post-analysis
myAnalysis.postAnalysis()

# Get Eigen-frequencies
print ("\nGetting results natural frequencies.....")
natrualFreq = myAnalysis.output.EigenFrequency

for mode, i in enumerate(natrualFreq):
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.2f} '.format(i) + "(Hz)")
