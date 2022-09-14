## [importModules]
# Import pyCAPS module
import pyCAPS

import os
import argparse
## [importModules]

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'MSES Pytest Example',
                                 prog = 'mses_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "./", nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

## [localVariable]
# Create working directory variable
workDir = "MSESAnalysisTest"
workDir = os.path.join(str(args.workDir[0]), workDir)
## [localVariable]

# Load CSM file
## [loadGeom]
geometryScript = os.path.join("..","csmData","airfoilSection.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)
## [loadGeom]

## [setGeom]
# Change a design parameter - camber in the geometry
myProblem.geometry.despmtr.camber = 0.0
## [setGeom]

## [loadMSES]
# Load mses aim
mses = myProblem.analysis.create(aim = "msesAIM")
## [loadMSES]

##[setMSES]
# Set Mach number, Reynolds number
mses.input.Mach = 0.5
mses.input.Re   = 5e6

# Set custom AoA
mses.input.Alpha = 3.0

# Trip the flow near the leading edge to get smooth gradient
mses.input.xTransition_Upper = 0.1
mses.input.xTransition_Lower = 0.1

# Set custom grid generation AoA
mses.input.GridAlpha = 3.0
##[setMSES]

## [results]
# Retrieve Cl, Cd, and Cm
Cl = mses.output.CL
print("Cl = ", Cl)

Cd = mses.output.CD
print("Cd = ", Cd)

Cm = mses.output.CM
print("Cm = ", Cm)
## [results]

if args.noPlotData == False:
    # Import pyplot module
    try:
        from matplotlib import pyplot as plt
        
        nmod = 4
        nplt = 3
        

        # Perform a grid convergence study
        Cl = []; Cl_Mach = []; Cl_mod = []; Cl_camber = []
        Cd = []; Cd_Mach = []; Cd_mod = []; Cd_camber = []
        Cm = []; Cm_Mach = []; Cm_mod = []; Cm_camber = []
        Airfoil_Points = [121, 131, 141, 151, 161, 171, 181, 191, 201, 209, 221, 231, 241, 251, 261, 271, 281, 291, 301, 311]
        for points in Airfoil_Points:
            print("Airfoil_Points =", points)
            mses.input.Airfoil_Points = points
            
            
            # Set camber as a design variable
            mses.input.Design_Variable = {"camber":{}}
            mses.input.Coarse_Iteration = 200
            mses.input.Cheby_Modes = None

            Cl.append(mses.output.CL); 
            Cd.append(mses.output.CD); 
            Cm.append(mses.output.CM); 

            Cl_camber.append(mses.output["CL"].deriv("camber"))
            Cd_camber.append(mses.output["CD"].deriv("camber"))
            Cm_camber.append(mses.output["CM"].deriv("camber"))

            Cl_Mach.append(mses.output["CL"].deriv("Mach"))
            Cd_Mach.append(mses.output["CD"].deriv("Mach"))
            Cm_Mach.append(mses.output["CM"].deriv("Mach"))
            
            mses.input.Coarse_Iteration = 0
            mses.input.Design_Variable = None
            mses.input.Cheby_Modes = [0]*nmod
            Cl_mod.append(mses.output["CL"].deriv("Cheby_Modes"))
            Cd_mod.append(mses.output["CD"].deriv("Cheby_Modes"))
            Cm_mod.append(mses.output["CM"].deriv("Cheby_Modes"))
            


        fig, axs = plt.subplots(3, nmod+nplt, sharex='col', figsize=(17,7))
        
        # Plot 
        axs[0,0].plot(Airfoil_Points, Cl, '-o'); 
        axs[0,0].set_ylabel("$C_L$")           ; 
        axs[0,1].plot(Airfoil_Points, Cl_Mach, '-o')
        axs[0,1].set_ylabel("$\partial C_L/ \partial \mathtt{Mach}$")
        axs[0,2].plot(Airfoil_Points, Cl_camber, '-o')
        axs[0,2].set_ylabel("$\partial C_L/ \partial \mathtt{camber}$")

        for i in range(nmod):
            mods = [None]*len(Airfoil_Points)
            for j in range(len(Airfoil_Points)):
                mods[j] = Cl_mod[j][i]
            axs[0,i+nplt].plot(Airfoil_Points, mods, '-o')
            axs[0,i+nplt].set_ylabel("$\partial C_L/ \partial \mathtt{mod"+str(i+1)+"}$")

        axs[1,0].plot(Airfoil_Points, Cd, '-o');
        axs[1,0].set_ylabel("$C_D$")           ;
        axs[1,1].plot(Airfoil_Points, Cd_Mach, '-o')
        axs[1,1].set_ylabel("$\partial C_D/ \partial \mathtt{Mach}$")
        axs[1,2].plot(Airfoil_Points, Cd_camber, '-o')
        axs[1,2].set_ylabel("$\partial C_D/ \partial \mathtt{camber}$")

        for i in range(nmod):
            mods = [None]*len(Airfoil_Points)
            for j in range(len(Airfoil_Points)):
                mods[j] = Cd_mod[j][i]
            axs[1,i+nplt].plot(Airfoil_Points, mods, '-o')
            axs[1,i+nplt].set_ylabel("$\partial C_D/ \partial \mathtt{mod"+str(i+1)+"}$")
         
        axs[2,0].plot(Airfoil_Points, Cm, '-o');
        axs[2,0].set_ylabel("$C_M$")           ; 
        axs[2,1].plot(Airfoil_Points, Cm_Mach, '-o')
        axs[2,1].set_ylabel("$\partial C_M/ \partial \mathtt{Mach}$")
        axs[2,2].plot(Airfoil_Points, Cm_camber, '-o')
        axs[2,2].set_ylabel("$\partial C_M/ \partial \mathtt{camber}$")

        for i in range(nmod):
            mods = [None]*len(Airfoil_Points)
            for j in range(len(Airfoil_Points)):
                mods[j] = Cd_mod[j][i]
            axs[2,i+nplt].plot(Airfoil_Points, mods, '-o')
            axs[2,i+nplt].set_ylabel("$\partial C_M/ \partial \mathtt{mod"+str(i+1)+"}$")
        
        for i in range(nmod+nplt):
            axs[2,i].set_xlabel("Airfoil Points")

        fig.suptitle("Grid Convergence of Coefficients and Derivatives\n (window must be closed to terminate Python script)")
        fig.tight_layout(pad=1.0)
        plt.show()

    except ImportError:
        print ("Unable to import matplotlib.pyplot module.")
