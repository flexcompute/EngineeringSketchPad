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
parser = argparse.ArgumentParser(description = 'Aeroelastic Cart3D and Mystran Example',
                                 prog = 'aeroelasticSimple_Iterative_Cart3D_and_Mystran',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = ["." + os.sep], nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-numberProc', default = 1, nargs=1, type=float, help = 'Number of processors')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AeroelasticSimple_Iterative_Cart3D_Mystran")

# Create projectName vairbale
projectName = "aeroelasticSimple_Iterative_CM"

# Set the number of transfer iterations
numTransferIteration = 4

# Load CSM file
geometryScript = os.path.join("..","csmData","aeroelasticDataTransferSimple.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# Load AIMs
myProblem.analysis.create(aim = "cart3dAIM",
                          name = "cart3d",
                          capsIntent = "Aerodynamic",
                          autoExec = False)

myProblem.analysis.create(aim = "mystranAIM",
                          name = "mystran",
                          capsIntent = "Structure",
                          autoExec = False)

# Create the data transfer connections
boundNames = ["Skin_Top", "Skin_Bottom", "Skin_Tip"]
for boundName in boundNames:
    # Create the bound
    bound = myProblem.bound.create(boundName)
    
    # Create the vertex sets on the bound for cart3d and mystran analysis
    cart3dVset  = bound.vertexSet.create(myProblem.analysis["cart3d"])
    mystranVset = bound.vertexSet.create(myProblem.analysis["mystran"])
    
    # Create pressure data sets
    cart3d_Pressure  = cart3dVset.dataSet.create("Pressure")
    mystran_Pressure = mystranVset.dataSet.create("Pressure")

    # Create displacement data sets
    cart3d_Displacement  = cart3dVset.dataSet.create("Displacement", init=[0,0,0])
    mystran_Displacement = mystranVset.dataSet.create("Displacement")

    # Link the data sets
    mystran_Pressure.link(cart3d_Pressure, "Conserve")
    cart3d_Displacement.link(mystran_Displacement, "Interpolate")
    
    # Close the bound as complete (cannot create more vertex or data sets)
    bound.close()

# Set inputs for cart3d
speedofSound = 340.0 # m/s
refVelocity = 100.0 # m/s
refDensity = 1.2 # kg/m^3

myProblem.analysis["cart3d"].input.Mach                  = refVelocity/speedofSound
myProblem.analysis["cart3d"].input.alpha                 = 2.0
myProblem.analysis["cart3d"].input.maxCycles             = 10
myProblem.analysis["cart3d"].input.nDiv                  = 6
myProblem.analysis["cart3d"].input.maxR                  = 9
myProblem.analysis["cart3d"].input.Pressure_Scale_Factor = 0.5*refDensity*refVelocity**2
myProblem.analysis["cart3d"].input.Tess_Params = [.05, .05, 15]

myProblem.analysis["cart3d"].input.y_is_spanwise = True
myProblem.analysis["cart3d"].input.Model_X_axis = "-Xb"
myProblem.analysis["cart3d"].input.Model_Y_axis = "Yb"
myProblem.analysis["cart3d"].input.Model_Z_axis = "-Zb"

# Set inputs for mystran
myProblem.analysis["mystran"].input.Proj_Name = projectName
myProblem.analysis["mystran"].input.Edge_Point_Min = 4
myProblem.analysis["mystran"].input.Edge_Point_Max = 10

myProblem.analysis["mystran"].input.Quad_Mesh = True
myProblem.analysis["mystran"].input.Tess_Params = [.05, .1, 15]
myProblem.analysis["mystran"].input.Analysis_Type = "Static"

# External pressure load to mystran that we will inherited from cart3d
load = {"loadType" : "PressureExternal"}
myProblem.analysis["mystran"].input.Load = {"pressureAero": load}

madeupium    = {"materialType" : "isotropic",
                "youngModulus" : 72.0E6 ,
                "poissonRatio": 0.33,
                "density" : 2.8E3}
myProblem.analysis["mystran"].input.Material = {"Madeupium": madeupium}

skin  = {"propertyType" : "Shell",
         "membraneThickness" : 0.06,
         "material"        : "madeupium",
         "bendingInertiaRatio" : 1.0, # Default
         "shearMembraneRatio"  : 5.0/6.0} # Default

ribSpar  = {"propertyType" : "Shell",
            "membraneThickness" : 0.6,
            "material"        : "madeupium",
            "bendingInertiaRatio" : 1.0, # Default
            "shearMembraneRatio"  : 5.0/6.0} # Default

myProblem.analysis["mystran"].input.Property = {"Skin"    : skin,
                                                "Rib_Root": ribSpar}

constraint = {"groupName" : "Rib_Root",
              "dofConstraint" : 123456}
myProblem.analysis["mystran"].input.Constraint = {"edgeConstraint": constraint}


# Aeroelastic iteration loop
for iter in range(numTransferIteration):

    if iter > 0:
        for boundName in boundNames:
            # Create the bound
            bound = myProblem.bound[boundName]
            
            # Get the vertex sets on the bound for cart3d analysis
            cart3dVset  = bound.vertexSet["cart3d"]
            mystranVset = bound.vertexSet["mystran"]
            
            # Get data sets
            cart3d_Displacement  = cart3dVset.dataSet["Displacement"]
            mystran_Displacement = mystranVset.dataSet["Displacement"]
            
            cart3d_Displacement.writeVTK(os.path.join(myProblem.analysis["mystran"].analysisDir,str(iter) + "_cart3d_Displacement"+boundName))
            mystran_Displacement.writeVTK(os.path.join(myProblem.analysis["mystran"].analysisDir,str(iter) + "_mystran_Displacement"+boundName))

    myProblem.analysis["cart3d"].runAnalysis()
    
    # Get Lift and Drag coefficients
    Cl = myProblem.analysis["cart3d"].output.C_L
    Cd = myProblem.analysis["cart3d"].output.C_D

    # Print lift and drag
    print("Cl = ", Cl)
    print("Cd = ", Cd)
    
    for boundName in boundNames:
        # Create the bound
        bound = myProblem.bound[boundName]
        
        # Get the vertex sets on the bound for cart3d analysis
        cart3dVset  = bound.vertexSet["cart3d"]
        mystranVset = bound.vertexSet["mystran"]
        
        # Get data sets
        cart3d_Pressure  = cart3dVset.dataSet["Pressure"]
        mystran_Pressure = mystranVset.dataSet["Pressure"]
        
        cart3d_Pressure.writeVTK(os.path.join(myProblem.analysis["cart3d"].analysisDir,str(iter) + "_cart3d_Pressure_"+boundName))
        mystran_Pressure.writeVTK(os.path.join(myProblem.analysis["cart3d"].analysisDir,str(iter) + "_mystran_Pressure"+boundName))

    myProblem.analysis["mystran"].runAnalysis()

