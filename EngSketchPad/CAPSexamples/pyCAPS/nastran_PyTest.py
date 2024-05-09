from __future__ import print_function

# Import pyCAPS class file
from pyCAPS import capsProblem

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Pytest Example',
                                 prog = 'nastran_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Initialize capsProblem object
myProblem = capsProblem()

# Load CSM file
geometryScript = os.path.join("..","csmData","feaWingBEM.csm")
myProblem.loadCAPS(geometryScript, verbosity=args.outLevel)

# Load egadsTess aim
myProblem.loadAIM( aim = "egadsTessAIM",
                   altName = "tess" )


# All triangles in the grid
myProblem.analysis["tess"].setAnalysisVal("Mesh_Elements", "Quad")

# Set global tessellation parameters
myProblem.analysis["tess"].setAnalysisVal("Tess_Params", [.05,.5,15])

# Generate the surface mesh
myProblem.analysis["tess"].preAnalysis()
myProblem.analysis["tess"].postAnalysis()

# Load nastran aim
myAnalysis = myProblem.loadAIM(aim = "nastranAIM",
                               altName = "nastran",
                               analysisDir= os.path.join(args.workDir[0], "NastranModalWingBEM"),
                               parents = ["tess"] )

projectName = "pyCAPS_nastran_Test"
myAnalysis.setAnalysisVal("Proj_Name", projectName)

myAnalysis.setAnalysisVal("Edge_Point_Max", 6)

myAnalysis.setAnalysisVal("Quad_Mesh", True)

myAnalysis.setAnalysisVal("File_Format", "Free")

myAnalysis.setAnalysisVal("Mesh_File_Format", "Large")

# Set analysis
eigen = { "extractionMethod"     : "Lanczos",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormalization"   : "MASS"}
myAnalysis.setAnalysisVal("Analysis", ("EigenAnalysis", eigen))

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.2E5 ,
                "poissonRatio" : .5,
                "density"      : 7850}
myAnalysis.setAnalysisVal("Material", [("Unobtainium", unobtainium),
                                       ("Madeupium", madeupium)])

# Set property
shell  = {"propertyType"      : "Shell",
          "material"          : "unobtainium",
          "bendingInertiaRatio" : 1.0, # Default - not necesssary
          "shearMembraneRatio"  : 0, # Turn of shear - no materialShear
          "membraneThickness" : 0.2}

myAnalysis.setAnalysisVal("Property", ("Ribs_and_Spars", shell))

# Set constraints
constraint = {"groupName" : ["Rib_Constraint"],
              "dofConstraint" : 123456}

myAnalysis.setAnalysisVal("Constraint", ("ribConstraint", constraint))

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
natrualFreq = myAnalysis.getAnalysisOutVal("EigenFrequency")

mode = 1
for i in natrualFreq:
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
    mode += 1

# Close CAPS
myProblem.closeCAPS()
