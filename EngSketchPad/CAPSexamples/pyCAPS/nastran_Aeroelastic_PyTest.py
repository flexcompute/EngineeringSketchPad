 # Import pyCAPS module
import pyCAPS

# Import os module
import os
import shutil
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

# Create project name
projectName = "NastranAeroWingBEM"

workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaWingBEMAero.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 5
egads.input.Edge_Point_Min = 5

# All quads in the grid
egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.05,.5,15]

# Load nastran aim
capsProblem.analysis.create(aim = "nastranAIM",
                          name = "nastran")

capsProblem.analysis["nastran"].input["Mesh"].link(egads.output["Surface_Mesh"])

capsProblem.analysis["nastran"].input.Proj_Name = "nastranAero"

# Set analysis type
capsProblem.analysis["nastran"].input.Analysis_Type = "Aeroelastic"

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


capsProblem.analysis["nastran"].input.Analysis = {"Trim1": trim1}

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

capsProblem.analysis["nastran"].input.Material = {"Unobtainium": unobtainium}

# Set property
shell  = {"propertyType"      : "Shell",
          "membraneThickness" : 0.2,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

shell2  = {"propertyType"      : "Shell",
          "membraneThickness" : 0.002,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

capsProblem.analysis["nastran"].input.Property = {"Ribs"    : shell,
                                                "Spar1"   : shell,
                                                "Spar2"   : shell,
                                                "Rib_Root": shell,
                                                "Skin"    : shell2}


# Defined Connections
connection = {    "dofDependent" : 123456,
                "connectionType" : "RigidBody"}

capsProblem.analysis["nastran"].input.Connect = {"Rib_Root": connection}


# Set constraints
constraint = {"groupName" : ["Rib_Root_Point"],
              "dofConstraint" : 12456}

capsProblem.analysis["nastran"].input.Constraint = {"ribConstraint": constraint}

# Set supports
support = {"groupName" : ["Rib_Root_Point"],
           "dofSupport": 3}

capsProblem.analysis["nastran"].input.Support = {"ribSupport": support}


# Aero
wing = {"groupName"         : "Wing",
        "numChord"          : 8,
        "numSpanPerSection" : 12}

# Note the surface name corresponds to the capsBound found in the *.csm file. This links
# the spline for the aerodynamic surface to the structural model
capsProblem.analysis["nastran"].input.VLM_Surface = {"Skin_Top": wing}

# Run AIM pre-analysis
capsProblem.analysis["nastran"].preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name +  ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(capsProblem.analysis["nastran"].analysisDir,capsProblem.analysis["nastran"].input.Proj_Name+".f06"))

print ("Done running Nastran!")

# Run AIM post-analysis
capsProblem.analysis["nastran"].postAnalysis()
