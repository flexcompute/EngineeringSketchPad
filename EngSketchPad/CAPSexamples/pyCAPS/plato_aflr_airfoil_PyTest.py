# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Plato Airfoil Pytest Example',
                                 prog = 'plato_aflr_airfoil_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Load CSM file
geometryScript = os.path.join("..","csmData","cfd2D.csm")
problem = pyCAPS.Problem(problemName = os.path.join(str(args.workDir[0]), "PlatoAirfoilTest"),
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use egads for surface tessellation
aflr2 = problem.analysis.create(aim='aflr2AIM', name='aflr2')

# Set new EGADS body tessellation parameters
aflr2.input.Tess_Params = [0.75, 0.001, 20.0]

# Turn on Quad meshing
aflr2.input.Mesh_Gen_Input_String = "mquad=1 mpp=3"

# Make the plato AIM
plato = problem.analysis.create(aim='platoAIM', name='plato')

plato.input["Mesh"].link(aflr2.output["Area_Mesh"])

plato.preAnalysis()
plato.postAnalysis()
