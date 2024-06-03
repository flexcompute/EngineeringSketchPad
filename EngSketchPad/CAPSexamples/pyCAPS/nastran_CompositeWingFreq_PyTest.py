## [import]
# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import shutil
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

# Create project name
projectName = "NastranCompositeWing_Freq"

# Working directory
workDir = os.path.join(str(args.workDir[0]), projectName)

## [geometry]
# Load CSM file
geometryScript = os.path.join("..","csmData","compositeWing.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)
## [geometry]

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 40
egads.input.Edge_Point_Min = 5

# All quads in the grid
egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.05,.5,15]


## [loadAIM]
# Load nastran aim
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")
## [loadAIM]



## [setInputs]
nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.Proj_Name        = "nastran_CompositeWing"
nastran.input.File_Format      = "Small"
nastran.input.Mesh_File_Format = "Large"
nastran.input.Analysis_Type    = "Modal"
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

nastran.input.Material = {"Aluminum": Aluminum,
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

#nastran.input.Property" = {"wing": aluminum}
nastran.input.Property = {"wing": composite}
## [defineProperties]

## [defineConstraints]
constraint = {"groupName" : "root",
              "dofConstraint" : 123456}

nastran.input.Constraint = {"root": constraint}
## [defineConstraints]

## [defineAnalysis]
eigen = { "extractionMethod"     : "MGIV",#"Lanczos",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormalization"   : "MASS"}

nastran.input.Analysis = {"EigenAnalysis": eigen}
## [defineAnalysis]

## [preAnalysis]
nastran.preAnalysis()
## [preAnalysis]

## [run]
####### Run Nastran####################
print ("\n\nRunning Nastran......" )

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name +  ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")
## [run]

## [postAnalysis]
# Run AIM post-analysis
nastran.postAnalysis()
## [postAnalysis]

# Get Eigen-frequencies
print ("\nGetting results for natural frequencies.....")
naturalFreq = capsProblem.analysis["nastran"].output.EigenFrequency

mode = 1
for i in naturalFreq:
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
    mode += 1

