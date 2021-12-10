# Import pyCAPS module
import pyCAPS

# Import modules
import os
import sys
import shutil
import argparse

# ASTROS_ROOT should be the path containing ASTRO.D01 and ASTRO.IDX
try:
   ASTROS_ROOT = os.environ["ASTROS_ROOT"]
   os.putenv("PATH", ASTROS_ROOT + os.pathsep + os.getenv("PATH"))
except KeyError:
   print("Please set the environment variable ASTROS_ROOT")
   sys.exit(1)

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Astros CompositeWingFreq PyTest Example',
                                 prog = 'astros_CompositeWingFreq_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "AstrosCompositeWing_Freq")

# Load CSM file
geometryScript = os.path.join("..","csmData","compositeWing.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Create astros aim
astros = myProblem.analysis.create( aim = "astrosAIM",
                                    name = "astros" )

# Setup bound information for visualization
numEigenVector = 5
eigenVectors = []
for i in range(numEigenVector):
    eigenVectors.append("EigenVector_" + str(i+1))

# Create the capsBounds for data transfer
boundNames = ["wing"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)
    
    # Create the vertex sets on the bound for astros analysis
    astrosVset = bound.vertexSet.create(astros)
    
    # Create eigenVector data sets
    for eigenVector in eigenVectors:
        astros_eigenVector = astrosVset.dataSet.create(eigenVector, pyCAPS.fType.FieldOut)

    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()

# Set project name so a mesh file is generated
projectName = "astros_CompositeWing"
astros.input.Proj_Name        = projectName
astros.input.Edge_Point_Max   = 10
astros.input.File_Format      = "Small"
astros.input.Mesh_File_Format = "Large"
astros.input.Analysis_Type    = "Modal"

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

astros.input.Material = {"Aluminum": Aluminum,
                         "Graphite_epoxy": Graphite_epoxy}

aluminum  = {"propertyType"         : "Shell",
             "material"             : "Aluminum",
             "bendingInertiaRatio"  : 1.0, # Default - not necesssary
             "shearMembraneRatio"   : 0, # Turn of shear - no materialShear
             "membraneThickness"    : 0.125 }

composite  = {"propertyType"           : "Composite",
              "shearBondAllowable"     : 1.0e6,
              "bendingInertiaRatio"    : 1.0, # Default - not necesssary
              "shearMembraneRatio"     : 0, # Turn off shear - no materialShear
              "compositeMaterial"      : ["Graphite_epoxy"]*8,
              "compositeThickness"     : [0.00525]*8,
              "compositeOrientation"   : [0, 0, 0, 0, -45, 45, -45, 45],
              "symmetricLaminate"      : True,
              "compositeFailureTheory" : "STRAIN" }

#astros.input.Property = {"wing": aluminum}
astros.input.Property = {"wing": composite}

constraint = {"groupName" : "root",
              "dofConstraint" : 123456}

astros.input.Constraint = {"root": constraint}

eigen = { "extractionMethod"     : "MGIV",
          "frequencyRange"       : [0, 1000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 10,
          "eigenNormaliztion"    : "MASS"}

astros.input.Analysis = {"EigenAnalysis": eigen}

# astros is executed automatically just-in-time

# Get Eigen-frequencies
print ("\nGetting results for natural frequencies.....")
natrualFreq = myProblem.analysis["astros"].output.EigenFrequency

mode = 1
for i in natrualFreq:
    print ("Natural freq (Mode {:d}) = ".format(mode) + '{:.5f} '.format(i) + "(Hz)")
    mode += 1

fp = open(astros.analysisDir + "/Data_All_Eigen.dat", "w")

for i in myProblem.bound.keys():
    for j in myProblem.bound[i].vertexSet.keys():
        for eigenVector in eigenVectors:
            myProblem.bound[i].vertexSet[j].dataSet[eigenVector].writeTecplot(file=fp)

fp.close()
