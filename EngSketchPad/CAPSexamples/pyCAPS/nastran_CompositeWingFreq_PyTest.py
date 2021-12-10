## [import]
# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse
## [import]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Composite Wing Frequency Pytest Example',
                                 prog = 'nastran_CompositeWingFreq_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "NastranCompositeWing_Freq")

## [geometry]
# Load CSM file
geometryScript = os.path.join("..","csmData","compositeWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [geometry]

## [loadAIM]
# Load nastran aim
nastranAIM = myProblem.analysis.create(aim = "nastranAIM",
                                       name = "nastran")
## [loadAIM]

## [setInputs]
# Set project name so a mesh file is generated
projectName = "nastran_CompositeWing"
nastranAIM.input.Proj_Name        = projectName
nastranAIM.input.Edge_Point_Max   = 40
nastranAIM.input.File_Format      = "Small"
nastranAIM.input.Mesh_File_Format = "Large"
nastranAIM.input.Analysis_Type    = "Modal"
## [setInputs]

## [defineMaterials]
Aluminum  = {"youngModulus" : 10.5E6 ,
             "poissonRatio" : 0.3,
             "density"      : 0.1/386,
             "shearModulus" : 4.04E6}

Graphite_epoxy = {"materialType"        : "Orthotropic",
                  "youngModulus"        : 20.8E6 ,
                  "youngModulusLateral" : 1.54E6,
                  "poissonRatio"        : 0.327,
                  "shearModulus"        : 0.80E6,
                  "density"             : 0.059/386,
                  "tensionAllow"        : 11.2e-3,
                  "tensionAllowLateral" : 4.7e-3,
                  "compressAllow"       : 11.2e-3,
                  "compressAllowLateral": 4.7e-3,
                  "shearAllow"          : 19.0e-3,
                  "allowType"           : 1}

nastranAIM.input.Material = {"Aluminum": Aluminum,
                             "Graphite_epoxy": Graphite_epoxy}
## [defineMaterials]

## [defineProperties]
aluminum  = {"propertyType"         : "Shell",
             "material"             : "Aluminum",
             "bendingInertiaRatio"  : 1.0, # Default - not necesssary
             "shearMembraneRatio"   : 0, # Turn of shear - no materialShear
             "membraneThickness"    : 0.125 }

composite  = {"propertyType"           : "Composite",
              "shearBondAllowable"       : 1.0e6,
              "bendingInertiaRatio"    : 1.0, # Default - not necesssary
              "shearMembraneRatio"     : 0, # Turn off shear - no materialShear
              "compositeMaterial"      : ["Graphite_epoxy"]*8,
              "compositeThickness"     : [0.00525]*8,
              "compositeOrientation"   : [0, 0, 0, 0, -45, 45, -45, 45],
              "symmetricLaminate"      : True,
              "compositeFailureTheory" : "STRN" }

#nastranAIM.input.Property = {"wing": aluminum}
nastranAIM.input.Property = {"wing": composite}
## [defineProperties]

## [defineConstraints]
constraint = {"groupName" : "root",
              "dofConstraint" : 123456}

nastranAIM.input.Constraint = {"root": constraint}
## [defineConstraints]

## [defineAnalysis]
eigen = { "extractionMethod"     : "MGIV",#"Lanczos",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormaliztion"    : "MASS"}

nastranAIM.input.Analysis = {"EigenAnalysis": eigen}
## [defineAnalysis]

## [preAnalysis]
nastranAIM.preAnalysis()
## [preAnalysis]

## [run]
####### Run Nastran####################
print ("\n\nRunning Nastran......" )
currentDirectory = os.getcwd() # Get our current working directory

os.chdir(myProblem.analysis["nastran"].analysisDir) # Move into test directory

if (args.noAnalysis == False):
    os.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + projectName +  ".dat"); # Run Nastran via system call

os.chdir(currentDirectory) # Move back to working directory
print ("Done running Nastran!")
## [run]

## [postAnalysis]
# Run AIM post-analysis
nastranAIM.postAnalysis()
## [postAnalysis]

# Get Eigen-frequencies
print ("\nGetting results for natural frequencies.....")
natrualFreq = myProblem.analysis["nastran"].output.EigenFrequency

for mode, i in enumerate(natrualFreq):
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
