# Import pyCAPS and os module
import pyCAPS
import os
import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'Cart3D Nacelle Pytest Example',
                                 prog = 'cart3d_Nacelle_PyTest',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

workDir = os.path.join(str(args.workDir[0]), "cart3dNacelle")

geometryScript = os.path.join("..","csmData","cfdNacelle.csm")
myProblem = pyCAPS.Problem(problemName=workDir,
                           capsFile=geometryScript,
                           outLevel=args.outLevel)

# -----------------------------------------------------------------
# Load desired aim
# -----------------------------------------------------------------
print ("Loading AIM")
cart3d = myProblem.analysis.create(aim = "cart3dAIM")
# -----------------------------------------------------------------

# Set Cart3d Inputs
cart3d.input.Mach      = 0.8
cart3d.input.alpha     = 0.0
cart3d.input.maxR      = 7
cart3d.input.Model_X_axis = '-Xb'
cart3d.input.Model_Y_axis = 'Yb'
cart3d.input.Model_Z_axis = '-Zb'
cart3d.input.y_is_spanwise = True
cart3d.input.Xslices = [0.0, 0.8, 2.2, 3.5]
cart3d.input.Yslices = [0.0]

# Setup C_D as the adaptation functional
# cart3d.input.nAdaptCycles = 5
# cart3d.input.Adapt_Functional = {"CD" : {"function" : "C_D", "power" : 1}}


# POWERED BCS -----------------------------------------------------

# inlet boundary conditions
inletbc = [None]*3

# Cart3d BC opt. 1: "SurfBC 1  1.0 0.3 0. 0. 0.714285  # compID   Rho    Xvel     Yvel    Zvel   Press"
inletbc[0] = {"bcType" : "Riemann", "staticDensity" : 1.0, "uVelocity": 0.3, "vVelocity":0, "wVelocity":0, "staticPressure":0.714285}
# Cart3d BC opt. 2: "InletPressRatioBC   1   1.35       # compID  Pt/Ptinf"
inletbc[1] = {"bcType" : "InletPressRatio", "totalPressure" : 1.35}
# Cart3d BC opt. 3: "InletVelocityBC     1   0.6         # compID   Vnorm/c_inf"
inletbc[2] = {"bcType" : "InletVelocity", "machNumber" : 0.6}


# exhaust boundary conditions
exhaustbc = [None]*3

# Cart3d BC opt. 1: "SurfBC 1  2.0 1.5 0. 0. 1.5  # compID   Rho    Xvel     Yvel    Zvel   Press"
exhaustbc[0] = {"bcType" : "Riemann", "staticDensity" : 2.0, "uVelocity": 1.5, "vVelocity":0, "wVelocity":0, "staticPressure":1.5}
# Cart3d BC opt. 1: "PowerStagBC         1   1.4 1.05   # compID   Pt/Ptinf  Tt/Tttinf"
exhaustbc[1] = {"bcType" : "PowerStag", "totalPressure" : 1.4, "totalTemperature" : 1.05}
# Cart3d BC opt. 2: "PowerMFRBC          1   0.78 1.05   # compID   MFR  Tt/Tttinf"
exhaustbc[2] = {"bcType" : "PowerMFR", "massflow" : 0.78, "totalTemperature" : 1.05}


# Run though all possible power BCs (any inlet/exhust combiniation is valid) 
for ibc in range(3):

    # Apply the selected BCs
    cart3d.input.Boundary_Condition = {"inlet"   : inletbc[ibc],
                                       "exhaust" : exhaustbc[ibc]}

# -----------------------------------------------------------------
# Cart3D auto-executes
# -----------------------------------------------------------------
    
    print ("C_A    " + str(cart3d.output.C_A))
    print ("C_Y    " + str(cart3d.output.C_Y))
    print ("C_N    " + str(cart3d.output.C_N))
    print ("C_D    " + str(cart3d.output.C_D))
    print ("C_S    " + str(cart3d.output.C_S))
    print ("C_L    " + str(cart3d.output.C_L))
    print ("C_l    " + str(cart3d.output.C_l))
    print ("C_m    " + str(cart3d.output.C_m))
    print ("C_n    " + str(cart3d.output.C_n))
    print ("C_M_x  " + str(cart3d.output.C_M_x))
    print ("C_M_y  " + str(cart3d.output.C_M_y))
    print ("C_M_z  " + str(cart3d.output.C_M_z))


