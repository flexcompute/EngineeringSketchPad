# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Pytest Example',
                                 prog = 'nastran_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
projectName = "NastranModalWingBEM"
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","feaWingBEM.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load egadsTess aim
myProblem.analysis.create( aim = "egadsTessAIM",
                           name = "tess" )

# No Tess vertexes on edges (minimial mesh)
myProblem.analysis["tess"].input.Edge_Point_Min = 2
myProblem.analysis["tess"].input.Edge_Point_Max = 2

# All quad grid
myProblem.analysis["tess"].input.Mesh_Elements = "Quad"

# Set global tessellation parameters
myProblem.analysis["tess"].input.Tess_Params = [1.,1.,30]

# Generate the surface mesh
myProblem.analysis["tess"].preAnalysis()
myProblem.analysis["tess"].postAnalysis()

# View the tessellation
#myProblem.analysis["tess"].geometry.view()

# Load nastran aim
myAnalysis = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")

# Link Surface_Mesh
myAnalysis.input["Mesh"].link(myProblem.analysis["tess"].output["Surface_Mesh"])

projectName = "pyCAPS_nastran_Test"
myAnalysis.input.Proj_Name = projectName

myAnalysis.input.Edge_Point_Max = 6

myAnalysis.input.Quad_Mesh = True

myAnalysis.input.File_Format = "Free"

myAnalysis.input.Mesh_File_Format = "Large"

# Set analysis
eigen = { "extractionMethod"     : "Lanczos",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormaliztion"    : "MASS"}
myAnalysis.input.Analysis = {"EigenAnalysis": eigen}

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.2E5 ,
                "poissonRatio" : .5,
                "density"      : 7850}
myAnalysis.input.Material = {"Unobtainium": unobtainium,
                             "Madeupium": madeupium}

# Set property
shell  = {"propertyType"      : "Shell",
          "material"          : "unobtainium",
          "bendingInertiaRatio" : 1.0, # Default - not necesssary
          "shearMembraneRatio"  : 0, # Turn of shear - no materialShear
          "membraneThickness" : 0.2}

myAnalysis.input.Property = {"Ribs_and_Spars": shell}

# Set constraints
constraint = {"groupName" : ["Rib_Constraint"],
              "dofConstraint" : 123456}

myAnalysis.input.Constraint = {"ribConstraint": constraint}


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

# Get Eigen-frequencies
print ("\nGetting results for natural frequencies.....")
natrualFreq = myAnalysis.output.EigenFrequency

for mode, i in enumerate(natrualFreq):
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
