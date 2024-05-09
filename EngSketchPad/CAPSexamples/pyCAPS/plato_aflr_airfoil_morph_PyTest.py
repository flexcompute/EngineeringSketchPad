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
problem = pyCAPS.Problem(problemName = os.path.join(str(args.workDir[0]), "PlatoAirfoilMorphTest"),
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use egads for surface tessellation
aflr2 = problem.analysis.create(aim='aflr2AIM', name='aflr2')

# Set new EGADS body tessellation parameters
aflr2.input.Tess_Params = [0.5, 0.001, 20.0]

# Make the plato AIM
plato = problem.analysis.create(aim='platoAIM', name='plato')

plato.input["Mesh"].link(aflr2.output["Area_Mesh"])

plato.input.Mesh_Morph = True

plato.preAnalysis()
plato.postAnalysis()

plato.input["Mesh"].unlink()

problem.geometry.despmtr.tunnelHeight *= 1.1

plato.preAnalysis()
plato.postAnalysis()
