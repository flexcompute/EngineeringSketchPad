# Import pyCAPS class file
import pyCAPS

# Import os module
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'FlightStream X43 Pytest Example',
                                 prog = 'cbaero_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Define units
deg  = pyCAPS.Unit("degree")
ft   = pyCAPS.Unit("feet")
s    = pyCAPS.Unit("s")
Pa   = pyCAPS.Unit("Pa")
mile = pyCAPS.Unit("mile")
hour = pyCAPS.Unit("hour")

# Create working directory variable
projectName = "FlightStreamX43Test"
workDir = os.path.join(str(args.workDir[0]), projectName)

# Load CSM file
geometryScript = os.path.join("..","csmData","cfdX43a.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)

# Turn of the farfield domain
capsProblem.geometry.cfgpmtr.makeBoundBox = 0

# Load egadsTess aim to create surface mesh
egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM",
                                        name = "egadsTess" )

# Set meshing parameters
egadsTess.input.Tess_Params = [0.025, 0.01, 10.0]

# Mixed quad/tri grid
egadsTess.input.Mesh_Elements = "Mixed"

# Load FlightStream aim
flightstream = capsProblem.analysis.create(aim = "flightstreamAIM",
                                           name = "FlightStream")

# Link mesh
flightstream.input["Surface_Mesh"].link(egadsTess.output["Surface_Mesh"])

# Specify the possibly full PATH to FlightStream executable
flightstream.input.FlightStream = "FlightStream"

# Set flow conditions
flightstream.input.Mach     = 4.0
flightstream.input.Alpha    = 4.0 * deg
flightstream.input.Altitude = 102100 * ft

# Reference velocity
flightstream.input.ReferenceVelocity = 4520 * mile/hour

# Explicitly run analysis (optional)
flightstream.runAnalysis()

# Get all ouputs
for i in flightstream.output.keys():
    print(str(i) + " = " + str(flightstream.output[i].value))
