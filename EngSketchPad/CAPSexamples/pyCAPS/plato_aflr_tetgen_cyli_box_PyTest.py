# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Plato Cyli-Box Pytest Example',
                                 prog = 'plato_aflr_tetgen_cyli_box_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Load CSM file
geometryScript = os.path.join("..","csmData","cyli_box_shell.csm")
problem = pyCAPS.Problem(problemName = os.path.join(str(args.workDir[0]), "PlatoAFLRTetGenCyliBoxTest"),
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use egads for surface tessellation
surface = problem.analysis.create(aim='aflr4AIM', name='aflr4')

# Mesh each body independently
surface.input.Multiple_Mesh = "MultiDomain"

# Load TetGen aim with the surface mesh as the parent
volume = problem.analysis.create(aim='tetgenAIM', name='tetgen')

# Link the surface mesh
volume.input["Surface_Mesh"].link(surface.output["Surface_Mesh"])

# Optional: Explicitly write mesh files
volume.input.Mesh_Format = 'vtk'

# Sepecify a point in each region with the different ID's
regions = {'a': { 'id' : 1, 'seed' : [ problem.geometry.despmtr.r1/2 , 0.0, 0.0] },
           'b': { 'id' : 2, 'seed' : [(problem.geometry.despmtr.r1+
                                       problem.geometry.despmtr.r2)/2, 0.0, 0.0] },
           'c': { 'id' : 3, 'seed' : [(problem.geometry.despmtr.r2+
                                       problem.geometry.despmtr.Lx)/2, 0.0, 0.0] }}
volume.input.Regions = regions

# Generate the volume mesh
volume.runAnalysis()

# Make the plato AIM
plato = problem.analysis.create(aim='platoAIM', name='plato')

plato.input["Mesh"].link(volume.output["Volume_Mesh"])

plato.preAnalysis()
plato.postAnalysis()
