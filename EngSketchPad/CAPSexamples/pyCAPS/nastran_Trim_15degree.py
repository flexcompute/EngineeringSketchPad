# Import pyCAPS module
import pyCAPS

# Import os module
import os
import shutil
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

# Create project name
projectName = "NastranTrim15Deg"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","15degreeWing.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 8
egads.input.Edge_Point_Min = 2

# All quads in the grid
egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.05,.5,15]

# Loadnastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")

nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.Proj_Name = "nastran_trim_15degree"

nastran.input.Mesh_File_Format = "Free"
nastran.input.File_Format = "Small"

# Set analysis type
nastran.input.Analysis_Type = "AeroelasticTrim"

# Set PARAM Inputs
nastran.input.Parameter = {"COUPMASS": "1",
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


nastran.input.Analysis = {"Trim": trim}

# Set materials
aluminum  = {   "youngModulus" : 10.3E6 , #psi
                "shearModulus" : 3.9E6 , # psi
                "density"      : 0.1} # lb/in^3

nastran.input.Material = {"aluminum": aluminum}

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


nastran.input.Property = {"Edge": shellEdge,
                          "Body": shell,
                          "Root": mass}


# Defined Connections
connection = {    "dofDependent" : 123456,
                "connectionType" : "RigidBody"}

nastran.input.Connect = {"Root": connection}


# Set constraints
constraint = {"groupName" : ["Root_Point"],
              "dofConstraint" : 1246}

nastran.input.Constraint = {"PointConstraint": constraint}

# Set supports
support = {"groupName" : ["Root_Point"],
           "dofSupport": 35}

nastran.input.Support = {"PointSupport": support}

# Force & Gravity Loads
load = {"groupName" : "Root_Point",
        "loadType"  : "GridForce",
        "forceScaleFactor" : 1.0E5,
        "directionVector" : [0.0, 0.0, 1.0]}

grav = {"loadType"  : "Gravity",
        "gravityAcceleration" : 386.0,
        "directionVector" : [0.0, 0.0, -1.0]}

nastran.input.Load = {"PointLoad": load, "GravityLoad": grav}


# Aero
wing = {"groupName"         : "Wing",
        "numChord"          : 4,
        "numSpanPerSection" : 6}

# Note the surface name corresponds to the capsBound found in the *.csm file. This links
# the spline for the aerodynamic surface to the structural model
nastran.input.VLM_Surface = {"WingSurface": wing}

# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
   nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " +nastran.input.Proj_Name + ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")

# Run AIM post-analysis
nastran.postAnalysis()
