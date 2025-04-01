# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 and AFLR3 Skip Volume Pytest Example',
                                 prog = 'aflr4_and_aflr3_SkipVolume_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["."+os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLRSkipVolumeTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","commonBox.csm")
problem = pyCAPS.Problem(problemName=workDir,
                         capsFile=geometryScript,
                         outLevel=args.outLevel)

# Use egads for surface tessellation
surface = problem.analysis.create(aim='aflr4AIM', name='S')

surface.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set maximum and minimum edge lengths relative to capsMeshLength
surface.input.max_scale = 0.1
surface.input.min_scale = 0.01

# Generate the surface mesh
surface.runAnalysis()

# Load TetGen aim with the surface mesh as the parent
volume = problem.analysis.create(aim='aflr3AIM', name='V')

# Link the surface mesh
volume.input["Surface_Mesh"].link(surface.output["Surface_Mesh"])

volume.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set the volume analysis values
volume.input.Proj_Name = 'volume'

Mesh_Format = ["VTK"]
if os.path.isfile( os.path.join(os.environ["ESP_ROOT"],"lib","exodusWriter.so") ):
    Mesh_Format.append("EXODUS")

volume.input.Mesh_Format = Mesh_Format

# Generate the volume mesh
volume.runAnalysis()
