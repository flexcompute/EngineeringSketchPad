# Import pyCAPS module
## [importModules]
import pyCAPS

import os
import argparse
## [importModules]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Masstran Pytest Example',
                                 prog = 'masstran_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noAnalysis', action='store_true', default = False, help = "Don't run analysis code")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "masstranWingBEM")

# Load CSM file
## [loadGeom]
geometryScript = os.path.join("..","csmData","feaWingBEM.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [loadGeom]

## [structureMesh]
# Load egadsTess aim
myProblem.analysis.create( aim = "egadsTessAIM",
                           name = "tess" )

# All triangles in the grid
myProblem.analysis["tess"].input.Mesh_Elements = "Quad"

# Set global tessellation parameters
myProblem.analysis["tess"].input.Tess_Params = [.05,.5,15]

# The surfaces mesh is generated automatically just-in-time
## [structureMesh]


# Load masstran aim
## [loadMasstran]
masstranAIM = myProblem.analysis.create(aim = "masstranAIM",
                                        name = "masstran")
## [loadMasstran]

##[setInputs]
masstranAIM.input["Surface_Mesh"].link(myProblem.analysis["tess"].output["Surface_Mesh"])
##[setInputs]

##[setTuple]
# Set materials
unobtainium  = {"youngModulus" : 2.2E11 ,
                "poissonRatio" : .33,
                "density"      : 7850}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 1.2E9 ,
                "poissonRatio" : .5,
                "density"      : 7850}
masstranAIM.input.Material = {"Unobtainium": unobtainium,
                              "Madeupium"  : madeupium}

# Set property
shell  = {"propertyType"      : "Shell",
          "membraneThickness" : 0.2,
          "bendingInertiaRatio" : 1.0, # Default
          "shearMembraneRatio"  : 5.0/6.0} # Default }

masstranAIM.input.Property = {"Ribs_and_Spars": shell}
##[setTuple]

## [executeAnalysis]
#####################################
## masstran executes automatically ##
#####################################
## [executeAnalysis]

## [results]
# Get mass properties
print ("\nGetting results mass properties.....\n")
Area     = masstranAIM.output.Area
Mass     = masstranAIM.output.Mass
Centroid = masstranAIM.output.Centroid
CG       = masstranAIM.output.CG
Ixx      = masstranAIM.output.Ixx
Iyy      = masstranAIM.output.Iyy
Izz      = masstranAIM.output.Izz
Ixy      = masstranAIM.output.Ixy
Ixz      = masstranAIM.output.Ixz
Iyz      = masstranAIM.output.Iyz
I        = masstranAIM.output.I_Vector
II       = masstranAIM.output.I_Tensor

print("Area     ", Area)
print("Mass     ", Mass)
print("Centroid ", Centroid)
print("CG       ", CG)
print("Ixx      ", Ixx)
print("Iyy      ", Iyy)
print("Izz      ", Izz)
print("Ixy      ", Ixy)
print("Ixz      ", Ixz)
print("Iyz      ", Iyz)
print("I        ", I)
print("II       ", II)
## [results]
