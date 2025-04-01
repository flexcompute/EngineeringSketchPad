
# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'EGADS Curve Tess Spheres PyTest Example',
                                 prog = 'curveTess_Spheres_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "CurveTessSpheresTest")

# Load CSM file and build the geometry explicitly
geometryScript = os.path.join("..","csmData","spheres.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)

# Load egadsTess aim
egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM")

# Set new EGADS body tessellation parameters
egadsTess.input.Tess_Params = [.1, 0.1, 20.0]

#egadsTess.runAnalysis()
#egadsTess.geometry.view()

# Load curveTess aim
curveTess = capsProblem.analysis.create(aim = "curveTessAIM")

curveTess.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

# Set the element order and class
curveTess.input.Element_Order = 2
curveTess.input.Element_Class = 1

# Optional: Explicitly write mesh files
Mesh_Format = []
if os.path.isfile( os.path.join(os.environ["ESP_ROOT"],"lib","exodusWriter.so") ):
    Mesh_Format.append("exodus")

    curveTess.input.Mesh_Format = Mesh_Format

# Generate the curved mesh
curveTess.runAnalysis()
#curveTess.geometry.view()
