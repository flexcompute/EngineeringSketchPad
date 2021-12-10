# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Flutter 15degree Example',
                                 prog = 'nastran_flutter_15degree',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "NastranFlutter15Deg")

# Load CSM file
geometryScript = os.path.join("..","csmData","15degreeWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load nastran aim
myAnalysis = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")

# Set project name so a mesh file is generated
projectName = "nastran_flutter_15degree"
myAnalysis.input.Proj_Name = projectName

myAnalysis.input.Mesh_File_Format = "Free"
myAnalysis.input.File_Format = "Small"

myAnalysis.input.Edge_Point_Max = 8

myAnalysis.input.Quad_Mesh = True

# Set analysis type
myAnalysis.input.Analysis_Type = "AeroelasticFlutter"
#myAnalysis.input.Analysis_Type = "Modal"

# Set PARAM Inputs
myAnalysis.input.Parameter = {"COUPMASS": "1",
                              "WTMASS"  :"0.0025901",
                              "AUNITS"  :"0.0025901",
                              "VREF"    :"12.0"}

# Set analysis
flutter = { "analysisType" : "AeroelasticFlutter",
            "extractionMethod"     : "Lanczos",
            "frequencyRange"       : [0.1, 300],
            "numEstEigenvalue"     : 4,
            "numDesiredEigenvalue" : 4,
            "eigenNormaliztion"    : "MASS",
            "aeroSymmetryXY" : "ASYM",
            "aeroSymmetryXZ" : "ANTISYM",
            "analysisConstraint" : ["PointConstraint"],
            "machNumber"     : 0.45,
            "dynamicPressure": 1.9,
            "density" : 0.967*1.1092e-7,
            "reducedFreq" : [0.001, 0.1, 0.12, 0.14, 0.16, 0.2],
          }


myAnalysis.input.Analysis = {"Flutter": flutter}

# Set materials
aluminum  = {   "youngModulus" : 9.2418E6 , #psi
                "shearModulus" : 3.4993E6 , # psi
                "density"      : 0.1} # lb/in^3

myAnalysis.input.Material = {"aluminum": aluminum}

# Set property
shell  = {"propertyType"        : "Shell",
          "material"            : "aluminum",
          "membraneThickness"   : 0.041 }

shellEdge  = {  "propertyType"        : "Shell",
                "material"            : "aluminum",
                "membraneThickness"   : 0.041/2 }

myAnalysis.input.Property = {"Edge": shellEdge,
                             "Body": shell}


# Defined Connections
connection = {    "dofDependent" : 123456,
                "connectionType" : "RigidBody"}

myAnalysis.input.Connect = {"Root": connection}


# Set constraints
constraint = {"groupName" : ["Root_Point"],
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"PointConstraint": constraint}

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
    #os.system("nastran.sh " + projectName +  ".dat"); # Run Nastran script to run nastran on mac va DREN
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

# Run AIM post-analysis
myAnalysis.postAnalysis()
