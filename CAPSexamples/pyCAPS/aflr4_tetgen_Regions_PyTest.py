# Import pyCAPS module
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 Tetgen Region Pytest Example',
                                 prog = 'aflr4_tetgen_Region_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLR4TetgenRegionTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","regions.csm")
problem = pyCAPS.Problem(problemName = workDir,
                         capsFile=geometryScript, 
                         outLevel=args.outLevel)

# Use egads for surface tessellation
surface = problem.analysis.create(aim='aflr4AIM', name='S')

# Set maximum and minimum edge lengths relative to capsMeshLength
surface.input.max_scale = 0.2
surface.input.min_scale = 0.01

# Load TetGen aim with the surface mesh as the parent
volume = problem.analysis.create(aim='tetgenAIM', name='V')

# Link the surface mesh
volume.input["Surface_Mesh"].link(surface.output["Surface_Mesh"])

# Set the volume analysis values
volume.input.Proj_Name = 'volume'
volume.input.Mesh_Format = 'VTK'

# Sepecify a point in each region with the different ID's
regions = {'a': { 'id' : 10, 'seed' : [0.0, 0.0, -1.0] },
           'b': { 'id' : 20, 'seed' : [0.0, 0.0,  1.0] }}
volume.input.Regions = regions

# Generate the volume mesh
volume.runAnalysis()
