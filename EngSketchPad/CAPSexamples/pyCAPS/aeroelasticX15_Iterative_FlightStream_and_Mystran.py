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
                                 prog = 'aeroelasticX15_Iterative_FlightStream_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AeroelasticX15_Iterative_FlightStream_Mystran")

# Create projectName vairbale
projectName = "aeroelasticX15_Iterative_FM"

# Set the number of transfer iterations
numTransferIteration = 3

# Load CSM file
geometryScript = os.path.join("..","csmData","x15","x15.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Turn on both OML and IML for FSI
capsProblem.geometry.cfgpmtr.VIEW.Concept = 0
capsProblem.geometry.cfgpmtr.VIEW.Oml = 1
capsProblem.geometry.cfgpmtr.VIEW.Iml = 1

deg  = pyCAPS.Unit("degree")
ft   = pyCAPS.Unit("feet")
m    = pyCAPS.Unit("m")
cm   = pyCAPS.Unit("cm")
mm   = pyCAPS.Unit("mm")
s    = pyCAPS.Unit("s")
g    = pyCAPS.Unit("g")
kg   = pyCAPS.Unit("kg")
K    = pyCAPS.Unit("Kelvin")
Pa   = pyCAPS.Unit("Pa")
GPa  = pyCAPS.Unit("GPa")
mile = pyCAPS.Unit("mile")
hour = pyCAPS.Unit("hour")

# Load AIMs
aflr4 = capsProblem.analysis.create(aim = "aflr4AIM",
                                    name="aflr4",
                                    capsIntent = "AerodynamicOML")

flightstream = capsProblem.analysis.create(aim = "flightstreamAIM",
                                           name = "FlightStream",
                                           capsIntent = "AerodynamicOML",
                                           autoExec = False)

egadsTess = capsProblem.analysis.create(aim = "egadsTessAIM",
                                        name= "egadsTess",
                                        capsIntent = "Structure")

mystran = capsProblem.analysis.create(aim = "mystranAIM",
                                      name = "mystran",
                                      capsIntent = "Structure",
                                      autoExec = False,
                                      unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

#
# Create the data transfer connections
boundNames = ["upperLeft", "upperRight", "lowerLeft", "lowerRight"]
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
    mystran_Pressure.link(FlightStream_Pressure, "Interpolate")
    FlightStream_Displacement.link(mystran_Displacement, "Interpolate")

    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()


# Set AIM verbosity
aflr4.input.Mesh_Quiet_Flag = True

# Set maximum and minimum edge lengths relative to capsMeshLength
aflr4.input.max_scale = 0.3
aflr4.input.min_scale = 0.03

# View the mesh
#aflr4.runAnalysis()
#aflr4.geometry.view()

# Set inputs for FlightStream
speedofSound = 302 * m/s
refDensity = 0.015504 * kg/m**3
refPressure = 1013 * Pa
refSonicVelocity = 1013 * Pa
refTemperature = 227 * K
refViscosity = 0.000014812 * Pa*s
refVelocity = 4520 * mile/hour

# Link the surface mesh
flightstream.input["Surface_Mesh"].link(aflr4.output["Surface_Mesh"])

# Specify the possibly full PATH to FlightStream executable
flightstream.input.FlightStream = "FlightStream"

# Set additional FlightStream export files
flightstream.input.Export_Solver_Analysis = ["Tecplot"]

# Set flow conditions
flightstream.input.Mach                  = 6.7
flightstream.input.Alpha                 = 8.0 * deg
flightstream.input.Pressure_Scale_Factor = -1

# Set fluid properties
#flightstream.input.Altitude = 102100 * ft
flightstream.input.Fluid_Properties = {"density": refDensity,
                                       "pressure": refPressure,
                                       "sonic_velocity": speedofSound,
                                       "temperature": refTemperature,
                                       "viscosity": refViscosity}

# Reference velocity
flightstream.input.ReferenceVelocity = refVelocity


# No Tess vertexes on edges (minimial mesh)
egadsTess.input.Edge_Point_Min = 4
egadsTess.input.Edge_Point_Max = 10

# Mixed quad/tri grid
egadsTess.input.Mesh_Elements = "Mixed"

# Set global tessellation parameters
egadsTess.input.Tess_Params = [.05, .1, 15]


# Set inputs for mystran

# Link Surface_Mesh
mystran.input["Mesh"].link(egadsTess.output["Surface_Mesh"])

# Static analysis
mystran.input.Analysis_Type = "Static"

# External pressure load to mystran that we will inherited from FlightStream
load = {"loadType" : "PressureExternal"}
mystran.input.Load = {"pressureAero": load}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0 * GPa,
                "poissonRatio": 0.33,
                "density" : 1 * g/cm**3}
unobtainium  = {"materialType" : "isotropic",
                "youngModulus" : 20.0 * GPa,
                "poissonRatio": 0.33,
                "density" : 0.5 * g/cm**3}

mystran.input.Material = {"madeupium": madeupium,
                          "unobtainium":unobtainium}

skin  = {"propertyType" : "Shell",
         "membraneThickness" : 1 * mm,
         "material"        : "unobtainium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

rib  = {"propertyType" : "Shell",
        "membraneThickness" : 2 * mm,
        "material"        : "madeupium",
        "bendingInertiaRatio" : 1.0, # Default
        "shearMembraneRatio"  : 5.0/6.0} # Default

spar  = {"propertyType" : "Shell",
        "membraneThickness" : 3 * mm,
        "material"        : "madeupium",
        "bendingInertiaRatio" : 1.0, # Default
        "shearMembraneRatio"  : 5.0/6.0} # Default

mystran.input.Property = {"Skin" : skin,
                          "Rib"  : rib,
                          "Spar" : spar}

mystran.input.Constraint = {"root": {"dofConstraint" : 12}, # Constrain x-y
                            "SOB" : {"dofConstraint" : 3}}  # Constrain z


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

            FlightStream_Displacement.writeTecplot(os.path.join(mystran.analysisDir,str(iter) + "_FlightStream_Displacement"+boundName))
            mystran_Displacement.writeTecplot(os.path.join(mystran.analysisDir,str(iter) + "_mystran_Displacement"+boundName))

    flightstream.runAnalysis()

    # Save away the output file
    shutil.copy(os.path.join(flightstream.analysisDir,flightstream.input.Proj_Name+".tec"), os.path.join(flightstream.analysisDir,flightstream.input.Proj_Name+f"_{iter:03}.tec"))

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

        FlightStream_Pressure.writeTecplot(os.path.join(flightstream.analysisDir,str(iter) + "_FlightStream_Pressure_"+boundName))
        mystran_Pressure.writeTecplot(os.path.join(flightstream.analysisDir,str(iter) + "_mystran_Pressure"+boundName))

    mystran.runAnalysis()

