# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AFLR4 Tip Treatment PyTest Example',
                                 prog = 'aflr4_TipTreat_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot surface meshes")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Working directory
workDir = os.path.join(str(args.workDir[0]), "AFLR4TipTreatAnalysisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","tiptreat.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript, 
                           outLevel=args.outLevel)

# Load AFLR4 aim
myAnalysis = myProblem.analysis.create(aim = "aflr4AIM")

# Set AIM outLevel
myAnalysis.input.Mesh_Quiet_Flag = True if args.outLevel == 0 else False

# Set output grid format since a project name is being supplied - Tecplot  file
myAnalysis.input.Mesh_Format = "Tecplot"

# Farfield growth factor
myAnalysis.input.ff_cdfr = 1.4

# Set maximum and minimum edge lengths relative to capsMeshLength
myAnalysis.input.max_scale = 0.2
myAnalysis.input.min_scale = 0.01

# Dissable curvature refinement when AFLR4 cannot generate a mesh
# myAnalysis.input.Mesh_Gen_Input_String = "auto_mode=0"

#myAnalysis.input.Mesh_Length_Factor = 0.25

# Run through a range of geometry configurations and generate surface meshes
BS = "BS"
begFacs  = [1] #[1, 10] # tip raduis ratio
sharptes = [0, 1] # sharp or blunt TE
conts    = [2, 0] # blend section continuity

for begFac in begFacs:
    for sharpte in sharptes:
        for cont in conts:

            myProblem.geometry.despmtr.begFac  = begFac
            myProblem.geometry.cfgpmtr.sharpte = sharpte
            myProblem.geometry.cfgpmtr.cont    = cont

            geom = str(begFac) + BS[sharpte] + str(cont)

            # Set project name so a mesh file is generated for each configuration
            myAnalysis.input.Proj_Name = "pyCAPS_AFLR4_" + geom + "_Test"

            # Save the gometry configuration
            filename = "wing" + geom + ".egads"
            myAnalysis.geometry.save(filename)

            # Run AIM
            myAnalysis.runAnalysis()

            if args.noPlotData == False:
                myAnalysis.geometry.view()
