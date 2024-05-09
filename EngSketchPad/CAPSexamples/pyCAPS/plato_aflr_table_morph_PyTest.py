# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Plato Table Morph Pytest Example',
                                 prog = 'plato_aflr_table_morph_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Load CSM file
geometryScript = os.path.join("..","csmData","table.csm")
problem = pyCAPS.Problem(problemName = os.path.join(str(args.workDir[0]), "PlatoAFLRMorphMultiDomainTest"),
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use AFLR4 for surface tessellation
surface = problem.analysis.create(aim='aflr4AIM', name='aflr4')

surface.input.Mesh_Length_Factor = 10

# # Use AFLR3 for volume meshing
volume = problem.analysis.create(aim='aflr3AIM', name='aflr3')

# Link the surface mesh
volume.input["Surface_Mesh"].link(surface.output["Surface_Mesh"])

# Generate a mesh for each body
volume.input.Multiple_Mesh = 'MultiDomain'

# Generate the volume mesh
volume.runAnalysis()

# Make the plato AIM
plato = problem.analysis.create(aim='platoAIM', name='plato')

plato.input["Mesh"].link(volume.output["Volume_Mesh"])

plato.input.Mesh_Morph = True

plato.preAnalysis()
plato.postAnalysis()

plato.input["Mesh"].unlink()

problem.geometry.despmtr.L *= 1.5

plato.preAnalysis()
plato.postAnalysis()

