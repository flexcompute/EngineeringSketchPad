from __future__ import print_function

# Import pyCAPS class file
import pyCAPS

# Import os module
import os
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

workDir = os.path.join(str(args.workDir[0]), projectName)

# Initialize CAPS Problem
geometryScript = os.path.join("..","csmData","feaSimplePlate.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load nastran aim
nastranAIM = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran",
                                       autoExec = False)

# Set project name so a mesh file is generated
nastranAIM.input.Proj_Name = projectName

# Set meshing parameters
nastranAIM.input.Edge_Point_Min = 3
nastranAIM.input.Edge_Point_Max = 4

nastranAIM.input.Quad_Mesh = True

nastranAIM.input.Tess_Params = [.25,.01,15]

# Set analysis type
nastranAIM.input.Analysis_Type = "Static"

# Set materials
madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E9 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}

nastranAIM.input.Material = {"Madeupium": madeupium}

# Set properties
shear  = {"propertyType" : "Shear",
          "membraneThickness" : 0.006,
          "material"        : "madeupium",
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default

nastranAIM.input.Property = {"plate": shear}

# Set constraints
constraint = {"groupName" : "plateEdge",
              "dofConstraint" : 123456}

nastranAIM.input.Constraint = {"edgeConstraint": constraint}

# Set grid moment load
momentLoad = {
    "groupName": "plate",
    "loadType": "GridMoment",
    "momentScaleFactor": 1.0
}

pressLoad = {
    "groupName" : "plate",
    "loadType" : "PressureDistribute",
    "pressureDistributeForce" : [1e6, 2.e6, 3e6, 4e6]
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
nastranAIM.input.Load = loads

# Set design variable
designVar = {
    "upperBound": 10.0,
    "lowerBound": -10.0,
}

# Set design variables
nastranAIM.input.Design_Variable = {"pltvar": designVar}

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
nastranAIM.input.Design_Constraint = {"pltcon": designCon}

# Set analysis
# No analysis case information needs to be set for a single static load case

# Run AIM pre-analysis
myProblem.analysis["nastran"].preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["nastran"].analysisDir) # Move into test directory

if (args.noAnalysis == False):
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat") # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")
######################################

# Run AIM post-analysis
myProblem.analysis["nastran"].postAnalysis()

