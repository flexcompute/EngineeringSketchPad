# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import shutil
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Miscellaneous PyTest Example',
                                 prog = 'nastran_PlateMisc_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()


# Create project name
projectName = "NastranMiscLoadPlate"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

# Initialize CAPS Problem
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
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

egads.input.Tess_Params = [.25,.01,15]


# Load nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")

nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.Proj_Name = projectName

# Set analysis type
nastran.input.Analysis_Type = "Static"

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}

nastran.input.Material = {"Madeupium": madeupium}

# Set properties
shear  = {"propertyType" : "Shear",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

nastran.input.Property = {"plate": shear}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

nastran.input.Constraint = {"edgeConstraint": constraint}

# Set grid moment load
momentLoad = {
    "groupName": "plate",
    "loadType": "GridMoment",
    "momentScaleFactor": 1.0,
    "directionVector": [0., 1., 0.]
}

pressLoad = {
    "groupName" : "plate",
    "loadType" : "Pressure",
    "pressureForce" : 1.e6
}

# Set thermal load
thermalLoad = {
    "groupName" : "plate",
    "loadType" : "Thermal",
    "temperate" : 99
}

# Set loads
loads = {
    "moment": momentLoad,
    "pressure": pressLoad,
    "thermal": thermalLoad
}
nastran.input.Load = loads

# Set design variable
designVar = {
    "upperBound": 10.0,
    "lowerBound": -10.0,
    "maxDelta": 0.5
}

# Set design variables
nastran.input.Design_Variable = {"pltvar": designVar}

designVarRel = {
    "componentType": "Property",
    "componentName": "plate",
    "fieldName": "TM",
    "variableName": "pltvar"
}

# Set design constraint
designCon = {
    "groupName": "plate",
    "upperBound": 10.0
}

# Set design constraints
nastran.input.Design_Constraint = {"pltcon": designCon}

# Set analysis
# No analysis case information needs to be set for a single static load case

# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name + ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")
######################################

# Run AIM post-analysis
nastran.postAnalysis()

