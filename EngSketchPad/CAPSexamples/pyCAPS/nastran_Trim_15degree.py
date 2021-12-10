# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Trim 15degree Example',
                                 prog = 'nastran_trim_15degree',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "NastranTrim15Deg")

# Load CSM file
geometryScript = os.path.join("..","csmData","15degreeWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load nastran aim
myAnalysis = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")

# Set project name so a mesh file is generated
projectName = "nastran_trim_15degree"
myAnalysis.input.Proj_Name = projectName

myAnalysis.input.Mesh_File_Format = "Free"
myAnalysis.input.File_Format = "Small"

myAnalysis.input.Edge_Point_Max = 8

myAnalysis.input.Quad_Mesh = True

# Set analysis type
myAnalysis.input.Analysis_Type = "AeroelasticTrim"

# Set PARAM Inputs
myAnalysis.input.Parameter = {"COUPMASS": "1",
                              "WTMASS"  : "0.0025901",
                              "AUNITS"  : "0.0025901"}

# Set analysis
trim = { "analysisType" : "AeroelasticTrim",
          "aeroSymmetryXY" : "ASYM",
          "aeroSymmetryXZ" : "SYM",
          "analysisConstraint" : ["PointConstraint"],
          "analysisSupport" : ["PointSupport"],
          "machNumber"     : 0.45,
          "dynamicPressure": 2.0,
          "density" : 1.0,
          "rigidVariable"  : ["URDD3","URDD5"],
          "rigidConstraint": ["ANGLEA","PITCH"],
          "magRigidConstraint" : [0.17453, 0.0],
          }


myAnalysis.input.Analysis = {"Trim": trim}

# Set materials
aluminum  = {   "youngModulus" : 10.3E6 , #psi
                "shearModulus" : 3.9E6 , # psi
                "density"      : 0.1} # lb/in^3

myAnalysis.input.Material = {"aluminum": aluminum}

# Set property
shell  = {"propertyType"        : "Shell",
          "material"            : "aluminum",
          "membraneThickness"   : 0.041 }

shellEdge  = {  "propertyType"        : "Shell",
                "material"            : "aluminum",
                "membraneThickness"   : 0.041/2 }

mass    =   {   "propertyType"  :   "ConcentratedMass",
                "mass"          :   1.0E5,
                "massInertia"   :   [0.0, 0.0, 1.0E5, 0.0, 0.0, 0.0]}


myAnalysis.input.Property = {"Edge": shellEdge,
                             "Body": shell,
                             "Root": mass}


# Defined Connections
connection = {    "dofDependent" : 123456,
                "connectionType" : "RigidBody"}

myAnalysis.input.Connect = {"Root": connection}


# Set constraints
constraint = {"groupName" : ["Root_Point"],
              "dofConstraint" : 1246}

myAnalysis.input.Constraint = {"PointConstraint": constraint}

# Set supports
support = {"groupName" : ["Root_Point"],
           "dofSupport": 35}

myAnalysis.input.Support = {"PointSupport": support}

# Force & Gravity Loads
load = {"groupName" : "Root_Point",
        "loadType"  : "GridForce",
        "forceScaleFactor" : 1.0E5,
        "directionVector" : [0.0, 0.0, 1.0]}

grav = {"loadType"  : "Gravity",
        "gravityAcceleration" : 386.0,
        "directionVector" : [0.0, 0.0, -1.0]}

myAnalysis.input.Load = {"PointLoad": load, "GravityLoad": grav}


# Aero
wing = {"groupName"         : "Wing",
        "numChord"          : 4,
        "numSpanPerSection" : 6}

# Note the surface name corresponds to the capsBound found in the *.csm file. This links
# the spline for the aerodynamic surface to the structural model
myAnalysis.input.VLM_Surface = {"WingSurface": wing}

# Run AIM pre-analysis
myAnalysis.preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myAnalysis.analysisDir) # Move into test directory

if (args.noAnalysis == False):
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

# Run AIM post-analysis
myAnalysis.postAnalysis()
