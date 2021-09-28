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

workDir = os.path.join(str(args.workDir[0]), "NastranAeroWingBEM")

# Load CSM file
geometryScript = os.path.join("..","csmData","feaWingBEMAero.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load nastran aim
myProblem.analysis.create(aim = "nastranAIM",
                          name = "nastran")

# Set project name so a mesh file is generated
projectName = "nastranAero"
myProblem.analysis["nastran"].input.Proj_Name = projectName

myProblem.analysis["nastran"].input.Edge_Point_Max = 4

myProblem.analysis["nastran"].input.Quad_Mesh = True

# Set analysis type
myProblem.analysis["nastran"].input.Analysis_Type = "Aeroelastic"

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


myProblem.analysis["nastran"].input.Analysis = {"Trim1": trim1}

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

myProblem.analysis["nastran"].input.Material = {"Unobtainium": unobtainium}

# Set property
shell  = {"propertyType"      : "Shell",
          "membraneThickness" : 0.2,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

shell2  = {"propertyType"      : "Shell",
          "membraneThickness" : 0.002,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

myProblem.analysis["nastran"].input.Property = {"Ribs"    : shell,
                                                "Spar1"   : shell,
                                                "Spar2"   : shell,
                                                "Rib_Root": shell,
                                                "Skin"    : shell2}


# Defined Connections
connection = {    "dofDependent" : 123456,
                "connectionType" : "RigidBody"}

myProblem.analysis["nastran"].input.Connect = {"Rib_Root": connection}


# Set constraints
constraint = {"groupName" : ["Rib_Root_Point"],
              "dofConstraint" : 12456}

myProblem.analysis["nastran"].input.Constraint = {"ribConstraint": constraint}

# Set supports
support = {"groupName" : ["Rib_Root_Point"],
           "dofSupport": 3}

myProblem.analysis["nastran"].input.Support = {"ribSupport": support}


# Aero
wing = {"groupName"         : "Wing",
        "numChord"          : 8,
        "numSpanPerSection" : 12}

# Note the surface name corresponds to the capsBound found in the *.csm file. This links
# the spline for the aerodynamic surface to the structural model
myProblem.analysis["nastran"].input.VLM_Surface = {"Skin_Top": wing}

# Run AIM pre-analysis
myProblem.analysis["nastran"].preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["nastran"].analysisDir) # Move into test directory

if (args.noAnalysis == False):
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")

# Run AIM post-analysis
myProblem.analysis["nastran"].postAnalysis()
