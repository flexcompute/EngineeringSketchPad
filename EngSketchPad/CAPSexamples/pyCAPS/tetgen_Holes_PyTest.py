# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Tetgen Holes Pytest Example',
                                 prog = 'tetgen_Holes_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "TetgenHolesTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","regions.csm")
problem = pyCAPS.Problem(problemName=workDir,
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use egads for surface tessellation
surface = problem.analysis.create(aim='egadsTessAIM', name='S')

# Set new EGADS body tessellation parameters
surface.input.Tess_Params = [.1, 0.001, 20.0]


# Load TetGen aim with the surface mesh as the parent
volume = problem.analysis.create(aim='tetgenAIM', name='V')

volume.input["Surface_Mesh"].link(surface.output["Surface_Mesh"])

# Set the volume analysis values
volume.input.Proj_Name = 'volume'
volume.input.Mesh_Format = 'Tecplot'

# Sepecify a hole point to remove a region in the mesh
holes = {'a': { 'seed' : [0.0, 0.0, -1.0] }}
volume.input.Holes = holes

# Generate the volume mesh
volume.runAnalysis()
