# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import shutil
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Nastran Composite Pytest Example',
                                 prog = 'nastran_Composite_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Initialize capsProblem object
projectName = "NastranCompositePlate"
workDir = os.path.join(str(args.workDir[0]), projectName)
geometryScript = os.path.join("..","csmData","feaCantileverPlate.csm")

capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
geometry = capsProblem.geometry

# Load egadsTess aim
egads = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set meshing parameters
egads.input.Edge_Point_Max = 5
egads.input.Edge_Point_Min = 5

# All quads in the grid
egads.input.Mesh_Elements = "Quad"

egads.input.Tess_Params = [.05,.5,15]

# Load CSM file
nastran = capsProblem.analysis.create(aim = "nastranAIM",
                                      name = "nastran")

nastran.input["Mesh"].link(egads.output["Surface_Mesh"])

nastran.input.Proj_Name = "pyCAPS_nastran_Test"
nastran.input.File_Format = "Small"
nastran.input.Mesh_File_Format = "Large"

# Set analysis
eigen = { "extractionMethod"     : "Lanczos",
          "frequencyRange"       : [0, 10000],
          "numEstEigenvalue"     : 1,
          "numDesiredEigenvalue" : 4,
          "eigenNormalization"   : "MASS"}
nastran.input.Analysis = {"EigenAnalysis": eigen}

# Set materials
unobtainium  = {"youngModulus" : 2.2E6 ,
                "poissonRatio" : .5,
                "density"      : 7850}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.2E5 ,
                "poissonRatio" : .5,
                "density"      : 7850}
nastran.input.Material = {"Unobtainium": unobtainium,
                       "Madeupium": madeupium}

# Set property
shell  = {"propertyType"           : "Composite",
          "shearBondAllowable"     : 1.0e6,
          "bendingInertiaRatio"    : 1.0, # Default - not necesssary
          "shearMembraneRatio"     : 0, # Turn off shear - no materialShear
          "membraneThickness"      : 0.2,
          "compositeMaterial"      : ["Unobtainium", "Madeupium", "Madeupium"],
          "compositeFailureTheory" : "STRN",
          "compositeThickness"     : [1.2, 0.5, 2.0],
          "compositeOrientation"   : [30.6, 45, 50.4]}

nastran.input.Property = {"plate": shell}

# Set constraints
constraint = {"groupName" : ["plateEdge"],
              "dofConstraint" : 123456}

nastran.input.Constraint = {"cantilever": constraint}


# Run AIM pre-analysis
nastran.preAnalysis()

####### Run Nastran####################
print ("\n\nRunning Nastran......")

if args.noAnalysis == False:
    nastran.system("nastran old=no notify=no batch=no scr=yes sdirectory=./ " + nastran.input.Proj_Name +  ".dat"); # Run Nastran via system call
else:
    # Copy old results if no analysis available
    shutil.copy2(os.path.join("..","analysisData","nastran",projectName+".f06"), 
                 os.path.join(nastran.analysisDir,nastran.input.Proj_Name+".f06"))

print ("Done running Nastran!")

# Run AIM post-analysis
nastran.postAnalysis()

