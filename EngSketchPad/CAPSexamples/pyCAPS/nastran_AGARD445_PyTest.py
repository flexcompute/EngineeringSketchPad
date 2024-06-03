# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import shutil
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran AGARD445.6 Pytest Example',
                                 prog = 'nastran_AGARD445_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()


# Create project name
projectName = "NastranModalAGARD445"
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaAGARD445.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Change the sweepAngle and span of the Geometry - Demo purposes
#capsProblem.geometry.despmtr.sweepAngle = 5 # From 45 to 5 degrees
#capsProblem.geometry.despmtr.semiSpan   = 5 # From 2.5 ft to 5 ft

# Meshing 
myMesh = capsProblem.analysis.create(aim = "egadsTessAIM", 
                                     name = "egadsTess" )
 
# Set meshing parameters
myMesh.input.Edge_Point_Max = 10
myMesh.input.Edge_Point_Min = 6

myMesh.input.Mesh_Elements = "Quad"

myMesh.input.Tess_Params = [.25,.01,15]

# Load nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")

# Set mesh
nastran.input["Mesh"].link(myMesh.output["Surface_Mesh"])

# Set project name so a mesh file is generated
nastran.input.Proj_Name = projectName

# Set meshing parameters
nastran.input.File_Format = "Free"

# Set analysis type
nastran.input.Analysis_Type = "Modal"

# Set analysis inputs
eigen = { "extractionMethod"     : "MGIV", # "Lanczos",
          "frequencyRange"       : [0.1, 200],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 2,
          "eigenNormalization"   : "MASS",
          "lanczosMode"          : 2,  # Default - not necesssary
          "lanczosType"          : "DPB"} # Default - not necesssary

nastran.input.Analysis = {"EigenAnalysis": eigen}

# Set materials
mahogany    = {"materialType"        : "orthotropic",
               "youngModulus"        : 0.457E6 ,
               "youngModulusLateral" : 0.0636E6,
               "poissonRatio"        : 0.31,
               "shearModulus"        : 0.0637E6,
               "shearModulusTrans1Z" : 0.00227E6,
               "shearModulusTrans2Z" : 0.00227E6,
               "density"             : 3.5742E-5}

nastran.input.Material = {"Mahogany": mahogany}

# Set properties
shell  = {"propertyType" : "Shell",
          "membraneThickness" : 0.82,
          "material"        : "mahogany",
          "bendingInertiaRatio" : 1.0, # Default - not necesssary
          "shearMembraneRatio"  : 5.0/6.0} # Default - not necesssary

nastran.input.Property = {"yatesPlate": shell}

# Set constraints
constraint = {"groupName" : "constEdge",
              "dofConstraint" : 123456}

nastran.input.Constraint = {"edgeConstraint": constraint}

# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran ####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))


print ("Done running Nastran!")
########################################

# Run AIM post-analysis
nastran.postAnalysis()

# Get Eigen-frequencies
print ("\nGetting results natural frequencies.....")
natrualFreq = nastran.output.EigenFrequency

for mode, i in enumerate(natrualFreq):
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.2f} '.format(i) + "(Hz)")
