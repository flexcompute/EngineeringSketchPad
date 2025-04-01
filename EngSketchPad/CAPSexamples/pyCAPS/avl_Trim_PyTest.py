# Import pyCAPS module
import pyCAPS

# Import os module
import os

# Import argparse
import argparse

from math import log, exp

#-----------------------------------------------------------
#
# Example of trimmed flight at a fixed Mach number using AVL
#
#-----------------------------------------------------------


# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(description = 'AVL Trim Pytest Example',
                                 prog = 'avl_Trim_PyTest.py',
                                 formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument('-workDir', default = "." + os.sep, nargs=1, type=str, help = 'Set working/run directory')
parser.add_argument("-outLevel", default = 1, type=int, choices=[0, 1, 2], help="Set output verbosity")
parser.add_argument('-noPlotData', action='store_true', default = False, help = "Don't plot data")
args = parser.parse_args()

# Create working directory variable
workDir = os.path.join(str(args.workDir[0]), "AVLTrimAnalsyisTest")

# Load CSM file
geometryScript = os.path.join("..","csmData","avlPlaneVanilla.csm")
capsProblem = pyCAPS.Problem(problemName=workDir,
                             capsFile=geometryScript,
                             outLevel=args.outLevel)

geometry = capsProblem.geometry


m    = pyCAPS.Unit("meter")
km   = pyCAPS.Unit("km")
kg   = pyCAPS.Unit("kg")
s    = pyCAPS.Unit("s")
K    = pyCAPS.Unit("Kelvin")
deg  = pyCAPS.Unit("degree")
ft   = pyCAPS.Unit("ft")
J    = pyCAPS.Unit("J")
Pa   = pyCAPS.Unit("Pa")
lbf  = pyCAPS.Unit("lbf")


def atmos(h):
#-----------------------------------------------------------------
#     Atmospheric functions  T(h), rho(h)...  valid to h = 47km
#
# Input:  h   = altitude above Sea Level (SL)  (km)
#
# Output: T   = temperature     (Kelvin)
#         p   = pressure        (Pa)
#         rho = density         (kg/m^3)
#         a   = speed of sound  (m/s)
#         mu  = viscosity       (kg/m-s)
#
#-    from Drela, "Flight Vehicle Aerodynamics", MIT Press 2014 
#-----------------------------------------------------------------

    pSL   = 1.0132e5 * Pa
    muSL  = 1.79e-5  * kg/(m*s)
    TSL   = 288.15   * K
    Tsuth = 110.0    * K
    
    cp   = 1004.0    * J/(kg*K)
    gam  = 1.4
    
    h = h/km
    
    T = (216.65 + 2.0*log( 1.0 + exp(35.75 - 3.25*h)
                              + exp(-3.0  + 0.0003*h**3) )) * K
    p = pSL*exp(-0.118*h - 0.0015*h**2/(1 - 0.018*h + 0.0011*h**2))
    
    rho = gam*p/((gam-1.0)*cp*T)
    a = (gam*p/rho)**(1/2)
    mu = muSL * (T/TSL)**(3/2) * (TSL+Tsuth)/(T+Tsuth)
    
    return (T,p,rho,a,mu)


def altp(p):
#------------------------------------------------------------------
#     Computes the inverse of p(h) function in subroutine atmos
#
# Input:  p = pressure        (Pa)
#
# Output: h = altitude above Sea Level (SL)  (km)
#
#------------------------------------------------------------------
#---- sea level pressure (1 bar = 1e5 Pa)
    pSL   = 1.0132e5 * Pa

#     p = pSL*exp(-0.118*h - 0.0015*h**2/(1 - 0.018*h + 0.0011*h**2))

    gprat = log(pSL/p)

#---- initial guess
    h = gprat/0.13

#---- Newton iteration to converge h
    for iter in range(10):
        res   = 0.118*h \
              + 0.0015*h**2 /(1 - 0.018*h + 0.0011*h**2) \
              - gprat
        res_h = 0.118 \
             + 0.0015*h*2.0/(1 - 0.018*h + 0.0011*h**2) \
             - 0.0015*h**2 /(1 - 0.018*h + 0.0011*h**2)**2 \
                          * (  - 0.018   + 0.0011*h*2.)
        
        dh = -res/res_h
        #print( res, dh)
        
        h = h + dh
        if (abs(dh) < 1.0e-8): 
            return h * km
    
    print('altp: pressure altitude convergence failed.  dh =', dh,' km')
    
    return h * km


# Load avl aim with units
print ("Loading AIM")
avl = capsProblem.analysis.create(aim = "avlAIM",
                                  unitSystem={"mass":kg, "length":m, "time":s, "temperature":K})

# Don't compute eigen values
avl.input.EigenValues = False

# Set Rotation rates
avl.input.RollRate  = 0.0
avl.input.PitchRate = 0.0
avl.input.YawRate   = 0.0

# Paracite drag
avl.input.CDp   = 0.022

# Gravity
g = 9.8 * m/s**2
avl.input.Gravity = g

# -----------------------------------------------------------------
# Set lifitng surface/slender body discretizations
# -----------------------------------------------------------------
fuse = {"groupName"    : "Fuselage",
        "numChord"     : 36}

wing = {"groupName"    : "Wing",
        "numChord"     : 8,
        "spaceChord"   : 1.0,
        "numSpanTotal" : 24}

htail = {"numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 16}

vtail = {"numChord"     : 8,
         "spaceChord"   : 1.0,
         "numSpanTotal" : 12}

avl.input.AVL_Surface = {"Wing": wing, "Htail": htail, "Vtail": vtail}


# -----------------------------------------------------------------
# Set control surface parameters
# -----------------------------------------------------------------
Aileron = {"deflectionAngle" : 0*deg,
           "leOrTe" : 1,
           "controlGain" : -1.0,
           "hingeLine" : [0, 0, 0],
           "deflectionDup"  : -1.0}

Rudder = {"deflectionAngle" : 0*deg,
           "leOrTe" : 1,
           "controlGain" : 1.0,
           "hingeLine" : [0, 0, 0],
           "deflectionDup"  : 1.0}

Stabilizer = {"deflectionAngle" : 0*deg,
              "leOrTe" : 0,
              "controlGain" : 1.0,
              "hingeLine" : [0, 1, 0],
              "deflectionDup"  : 1.0}

avl.input.AVL_Control = {"Aileron": Aileron, "Stabilizer":Stabilizer, "Rudder":Rudder}

# -----------------------------------------------------------------
# Trimmed analysis
# -----------------------------------------------------------------

W     = 80 * lbf # Plane Vanilla weight
CL0   = 1.1      # Reference CL
Mach0 = 0.2      # Reference Mach number
gam   = 1.4
S     = geometry.despmtr.wing.area * ft**2

p0 = W/(0.5*gam*Mach0**2*S*CL0)

h0 = altp(p0)

avl.input.Mach = Mach0
CGx = geometry.despmtr.wing.xroot + geometry.despmtr.wing.croot/4
Ixx = 1.350 * kg*m**2
Iyy = 0.7509 * kg*m**2
Izz = 2.095 * kg*m**2
Ixy = 0 * kg*m**2
Ixz = 0 * kg*m**2
Iyz = 0 * kg*m**2

CLmin = CL0-0.5
CLmax = CL0+0.5

if args.noPlotData == False:
    import matplotlib.pyplot as plt
    
    fig, axs = plt.subplots(5, sharex=True)
    lns = []

for dx in [-geometry.despmtr.wing.croot*0.2, 0, geometry.despmtr.wing.croot*0.2]:

    avl.input.MassProp = {"Aircraft":{"mass":W/g, "CG":[CGx+dx,0,0]*ft, "massInertia":[Ixx, Iyy, Izz, Ixy, Ixz, Iyz]}}

    CL = []
    LoD = []
    h = []
    CDi = []
    e = []
    stabIncidence = [] 
    
    for i in range(10):
        t = i/9.0
        CLi = CLmin*(1-t) + t*CLmax
        
        p = p0*CL0/CLi
        alt = altp(p)
        (T,p,rho,a,mu) = atmos(alt)
        
        # Set trim condition
        Operation = {"Alpha":{"CL":CLi},
                     "Stabilizer":{"Cm":0.0}
                     }
        avl.input.AVL_Operation = Operation
    
        avl.input.Density  = rho
        avl.input.Velocity = Mach0*a
        
        CLi = avl.output["CLtot"].value
        CD  = avl.output["CDtot"].value
        
        CL.append(CLi)
        h.append(alt/km)
        LoD.append(CLi/CD)
        CDi.append(CD)
        e.append(avl.output["e"].value)
        stabIncidence.append(avl.output["ControlDeflection"].value["Stabilizer"])

    if args.noPlotData == False:
        
        lns += axs[0].plot(CL, h, label = f"dGCx {dx:.2f}")
        axs[0].set_ylabel(r"Altitude [km]")
        
        axs[1].plot(CL, LoD)
        axs[1].set_ylabel(r"L/D")
        
        axs[2].plot(CL,CDi)
        axs[2].set_ylabel(r"$C_{Di}$")
        
        axs[3].plot(CL,e)
        axs[3].set_ylabel(r"e")
        
        axs[4].plot(CL,stabIncidence)
        axs[4].set_ylabel(r"$\delta_{stab}$ [deg]")
        
    
if args.noPlotData == False:
    axs[-1].set_xlabel(r"$C_L$")
    labs = [l.get_label() for l in lns]
    axs[0].legend(lns, labs, loc='lower center', bbox_to_anchor=(0.5, 1.1), fancybox=True, ncol=3)

    plt.show()

