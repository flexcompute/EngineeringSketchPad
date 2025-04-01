# Import pyCAPS module
import pyCAPS

# Import os module
import os
import sys

# Import shutil module
import shutil

# Import argparse module
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Aeroelastic FlightStream and Mystran Example',
                                 prog = 'aeroelasticSimple_Iterative_FlightStream_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AeroelasticSimple_Iterative_FlightStream_Mystran")

# Create projectName vairbale
projectName = "aeroelasticSimple_Iterative_FM"

# Set the number of transfer iterations
numTransferIteration = 4

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransferSimple.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                    name="aflr4",
                                    capsIntent = "AerodynamicOML")

flightstream = capsProblem.analysis.create(aim = "flightstreamAIM",
                                           name = "FlightStream",
                                           capsIntent = "AerodynamicOML",
                                           autoExec = False)

mystran = capsProblem.analysis.create(aim = "mystranAIM",
                                      name = "mystran",
                                      capsIntent = "Structure",
                                      autoExec = False)

#
# Create the data transfer connections
boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = capsProblem.bound.create(boundName)

    # Create the vertex sets on the bound for FlightStream and mystran analysis
    FlightStreamVset  = bound.vertexSet.create(flightstream)
    mystranVset = bound.vertexSet.create(mystran)

    # Create pressure data sets
    FlightStream_Pressure  = FlightStreamVset.dataSet.create("Pressure")
    mystran_Pressure = mystranVset.dataSet.create("Pressure")

    # Create displacement data sets
    FlightStream_Displacement  = FlightStreamVset.dataSet.create("Displacement", init=[0,0,0])
    mystran_Displacement = mystranVset.dataSet.create("Displacement")

    # Link the data sets
    mystran_Pressure.link(FlightStream_Pressure, "Conserve")
    FlightStream_Displacement.link(mystran_Displacement, "Interpolate")

    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()


# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.4
aflr4.input.min_scale = 0.1

# View the mesh
#aflr4.runAnalysis()
#aflr4.geometry.view()

deg = pyCAPS.Unit("degree")
ft  = pyCAPS.Unit("feet")
m   = pyCAPS.Unit("m")
s   = pyCAPS.Unit("s")
kg  = pyCAPS.Unit("kg")

# Set inputs for FlightStream
speedofSound = 340.0 * m/s
refVelocity = 100.0 * m/s
refDensity = 1.2 * kg/m**3

# Link the surface mesh
flightstream.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

flightstream.input.Mach                  = 7
flightstream.input.Alpha                 = 2.0 * deg
flightstream.input.Pressure_Scale_Factor = 1

# Set Altitude
flightstream.input.Altitude = 15000 * ft

flightstream.input.ReferenceVelocity = 2240.0 * m/s


# Set inputs for mystran
mystran.input.Proj_Name = projectName
mystran.input.Edge_Point_Min = 4
mystran.input.Edge_Point_Max = 10

mystran.input.Quad_Mesh = True
mystran.input.Tess_Params = [.05, .1, 15]
mystran.input.Analysis_Type = "Static"

# External pressure load to mystran that we will inherited from FlightStream
load = {"loadType" : "PressureExternal"}
mystran.input.Load = {"pressureAero": load}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E6 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}
mystran.input.Material = {"Madeupium": madeupium}

skin  = {"propertyType" : "Shell",
         "membraneThickness" : 0.6,
         "material"        : "madeupium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

ribSpar  = {"propertyType" : "Shell",
            "membraneThickness" : 6,
            "material"        : "madeupium",
            "bendingInertiaRatio" : 1.0, # Default
            "shearMembraneRatio"  : 5.0/6.0} # Default

mystran.input.Property = {"Skin"    : skin,
                                                "Rib_Root": ribSpar}

constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
mystran.input.Constraint = {"edgeConstraint": constraint}


# Aeroelastic iteration loop
for iter in range(numTransferIteration):

    if iter > 0:
        for boundName in boundNames:
            # Create the bound
            bound = capsProblem.bound[boundName]

            # Get the vertex sets on the bound for FlightStream analysis
            FlightStreamVset  = bound.vertexSet["FlightStream"]
            mystranVset = bound.vertexSet["mystran"]

            # Get data sets
            FlightStream_Displacement  = FlightStreamVset.dataSet["Displacement"]
            mystran_Displacement = mystranVset.dataSet["Displacement"]

            FlightStream_Displacement.writeVTK(os.path.join(mystran.analysisDir,str(iter) + "_FlightStream_Displacement"+boundName))
            mystran_Displacement.writeVTK(os.path.join(mystran.analysisDir,str(iter) + "_mystran_Displacement"+boundName))

    # Execute flightstream
    flightstream.runAnalysis()

    # Get Lift and Drag coefficients
    Cl = flightstream.output.CL
    Cd = flightstream.output.CDi

    # Print lift and drag
    print("Cl = ", Cl)
    print("Cd = ", Cd)

    for boundName in boundNames:
        # Create the bound
        bound = capsProblem.bound[boundName]

        # Get the vertex sets on the bound for FlightStream analysis
        FlightStreamVset  = bound.vertexSet["FlightStream"]
        mystranVset = bound.vertexSet["mystran"]

        # Get data sets
        FlightStream_Pressure  = FlightStreamVset.dataSet["Pressure"]
        mystran_Pressure = mystranVset.dataSet["Pressure"]

        FlightStream_Pressure.writeVTK(os.path.join(flightstream.analysisDir,str(iter) + "_FlightStream_Pressure_"+boundName))
        mystran_Pressure.writeVTK(os.path.join(flightstream.analysisDir,str(iter) + "_mystran_Pressure"+boundName))

    # Execute mystran
    mystran.runAnalysis()

