import os
import numpy as np
pi = np.pi
import pyCAPS
from corsairlite.optimization.solve import solve
from corsairlite.optimization import Variable, RuntimeConstraint, Constant, Formulation
from corsairlite.core.data.standardAtmosphere import atm
from corsairlite import units
from fullModelNACA import fullModel

import argparse

# Setup and read command line options
parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)

# Set up available commandline options
parser.add_argument("-outLevel", default=0, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

# Get the pyCAPS Problem object
myProblem = pyCAPS.Problem("hoburg", 
                           capsFile=os.path.join("csm", "naca.csm"), 
                           phaseName="CMConstraint", 
                           phaseStart="MSES", 
                           phaseContinuation=True,
                           outLevel=args.outLevel)
myProblem.intentPhrase(["Add CM constraint to reduce aft camber"])

myProblem.geometry.cfgpmtr["view:MSES"].value = 1

# Define optimization formulation
N_segments = 3
# Constants
A_prop          = Constant(name = "A_prop",          value = 0.785,                 units = "m^2",         description = "Disk area of propeller")
CDA0            = Constant(name = "CDA0",            value = 0.05,                  units = "m^2",         description = "fuselage drag area")
C_Lmax          = Constant(name = "C_Lmax",          value = 1.5,                   units = "-",           description = "max CL with flaps down")
e               = Constant(name = "e",               value = 0.95,                  units = "-",           description = "Oswald efficiency factor")
eta_eng         = Constant(name = "eta_eng",         value = 0.35,                  units = "",            description = "Engine Efficiency")
eta_v           = Constant(name = "eta_v",           value = 0.85,                  units = "",            description = "Propeller Viscous Efficiency")
f_wadd          = Constant(name = "f_wadd",          value = 2.0,                   units = "",            description = "Added Weight Fraction")
g               = Constant(name = "g",               value = 9.81,                  units = "m/s^2",       description = "Gravitational constant")
h_fuel          = Constant(name = "h_fuel",          value = 46e6,                  units = "J/kg",        description = "Fuel Specific energy density")
k_ew            = Constant(name = "k_ew",            value = 0.0372,                units = "N/W^(0.803)", description = "Constant for engine weight")
mu              = Constant(name = "mu",              value = atm.mu(3000*units.m),  units = "kg/m/s",      description = "viscosity of air")
N_lift          = Constant(name = "N_lift",          value = 6.0,                   units = "-",           description = "Ultimate Load Factor" )
r_h             = Constant(name = "r_h",             value = 0.75,                  units = "-",           description = "Ratio of height at rear spar to maximum wing height")
rho             = Constant(name = "rho",             value = atm.rho(3000*units.m), units = "kg/m^3",      description = "density of air")
rho_cap         = Constant(name = "rho_cap",         value = 2700,                  units = "kg/m^3",      description = "Density of wing cap material (aluminum)") 
rho_SL          = Constant(name = "rho_SL",          value = atm.rho(0*units.m),    units = "kg/m^3",      description = "density of air, sea level")
rho_web         = Constant(name = "rho_web",         value = 2700,                  units = "kg/m^3",      description = "Density of wing web material (aluminum)") 
sigma_max       = Constant(name = "sigma_max",       value = 250,                   units = "MPa",         description = "Allowable tensile stress of aluminum" )
sigma_max_shear = Constant(name = "sigma_max_shear", value = 167,                   units = "MPa",         description = "Allowable shear stress of aluminum" )
w_bar           = Constant(name = "w_bar",           value = .5,                    units = "-",           description = "Ratio of spar box width to chord length" )
W_fixed         = Constant(name = "W_fixed",         value = 14700,                 units = "N",           description = "fixed weight")

# Vector Variables
V        = []
C_L      = []
C_D      = []
C_Dfuse  = []
C_Dp     = []
C_Di     = []
C_f      = []
T        = []
W        = []
Re       = []
eta_i    = []
eta_prop = []
eta_0    = []
z_bre    = []
for i in range(0,N_segments):
    V.append(        Variable(name="V_%d"%(i),        guess=1.0, units = "m/s", description="Velocity, segment %d"%(i)))
    C_L.append(      Variable(name="C_L_%d"%(i),      guess=1.0, units = "",    description="Lift Coefficient, segment %d"%(i)))
    C_D.append(      Variable(name="C_D_%d"%(i),      guess=1.0, units = "",    description="Drag Coefficient, segment %d"%(i)))
    C_Dfuse.append(  Variable(name="C_Dfuse_%d"%(i),  guess=1.0, units = "",    description="Fuselage Drag Coefficient, segment %d"%(i)))
    C_Dp.append(     Variable(name="C_Dp_%d"%(i),     guess=1.0, units = "",    description="Wing Profile Drag Coefficient, segment %d"%(i)))
    C_Di.append(     Variable(name="C_Di_%d"%(i),     guess=1.0, units = "",    description="Induced Drag Coefficient, segment %d"%(i)))
    C_f.append(      Variable(name="C_f_%d"%(i),      guess=1.0, units = "",    description="Friction Coefficient, segment %d"%(i)))
    T.append(        Variable(name="T_%d"%(i),        guess=1.0, units = "N",   description="Thrust, segment %d"%(i)))
    W.append(        Variable(name="W_%d"%(i),        guess=1.0, units = "N",   description="Weight, segment %d"%(i)))
    Re.append(       Variable(name="Re_%d"%(i),       guess=1e7, units = "",    description="Reynolds Number, segment %d"%(i)))
    eta_i.append(    Variable(name="eta_i_%d"%(i),    guess=1.0, units = "",    description="Inviscid Propeller Efficiency, segment %d"%(i)))
    eta_prop.append( Variable(name="eta_prop_%d"%(i), guess=1.0, units = "",    description="Propeller Efficiency, segment %d"%(i)))
    eta_0.append(    Variable(name="eta_0_%d"%(i),    guess=1.0, units = "",    description="Overall Efficiency, segment %d"%(i)))
    z_bre.append(    Variable(name="z_bre_%d"%(i),    guess=1.0, units = "",    description="Breguet Range Factor, segment %d"%(i)))

# Free Variables
AR         = Variable(name = "AR",         guess = 10.0,   units = "-",   description = "aspect ratio")
I_cap_bar  = Variable(name = "I_cap_bar",  guess = 1,      units = "-",   description = "Area moment of inertia of cap on 2D cross section, normalized by chord^4" )
M_r_bar    = Variable(name = "M_r_bar",    guess = 1000,   units = "N",   description = "Root Bending unit per unit chord" )
nu         = Variable(name = "nu",         guess = 0.5,    units = "-",   description = "Placeholder, (1+lam_w+lam_w**2)/(1+lam_w)**2" )
p          = Variable(name = "p",          guess = 3,      units = "-",   description = "Dummy Variable (1+2*lam_w)" )
P_max      = Variable(name = "P_max",      guess = 1000.0, units = "W",   description = "Maximum Engine Power")
q          = Variable(name = "q",          guess = 2,      units = "-",   description = "Dummy Variable (1+lam_w)" )
R          = Variable(name = "R",          guess = 5000,   units = "km",  description = "Single Segment range")
S          = Variable(name = "S",          guess = 10.0,   units = "m^2", description = "total wing area")
t_cap_bar  = Variable(name = "t_cap_bar",  guess = .05,    units = "-",   description = "Spar cap thickness per unit chord" )
t_web_bar  = Variable(name = "t_web_bar",  guess = .05,    units = "-",   description = "Spar web thickness per unit chord" )
tau        = Variable(name = "tau",        guess = 0.15,   units = "-",   description = "airfoil thickness to chord ratio")
V_stall    = Variable(name = "V_stall",    guess = 30.0,   units = "m/s", description = "stall speed")
W_cap      = Variable(name = "W_cap",      guess = 600,    units = "N",   description = "Weight of Wing Spar Cap" )
W_eng      = Variable(name = "W_eng",      guess = 5000,   units = "N",   description = "Engine Weight")
W_fuel_out = Variable(name = "W_fuel_out", guess = 2500,   units = "N",   description = "Weight of fuel, outbound")
W_fuel_ret = Variable(name = "W_fuel_ret", guess = 2500,   units = "N",   description = "Weight of fuel, return")
W_MTO      = Variable(name = "W_MTO",      guess = 2500,   units = "N",   description = "Maximum Takeoff Weight")
W_pay      = Variable(name = "W_pay",      guess = 5000,   units = "N",   description = "Payload Weight")
W_tilde    = Variable(name = "W_tilde",    guess = 5000,   units = "N",   description = "Dry Weight, no wing")
W_web      = Variable(name = "W_web",      guess = 400,    units = "N",   description = "Weight of Wing Spar Shear Web" )
W_wing     = Variable(name = "W_wing",     guess = 2500,   units = "N",   description = "wing weight")
W_zfw      = Variable(name = "W_zfw",      guess = 5000,   units = "N",   description = "Zero Fuel Weight")

# Add MSES and camber relevant variables
maxCamber      = Variable(name = "maxCamber",        guess = 0.02,  units = "-",   description = "maximum airfoil camber scaled by chord")
camberLocation = Variable(name = "camberLocation",   guess = 0.5,   units = "-",   description = "location of maximum airfoil camber scaled by chord")
Alpha_p20_MSES = []
C_L_MSES = []  
negCMp1 = []
for i in range(0,N_segments):
    C_L_MSES.append(       Variable(name="C_L_MSES_%d"%(i),        guess=1.0,   units = "-",     description="Lift Coefficient from MSES, segment %d"%(i)))
    Alpha_p20_MSES.append( Variable(name="Alpha_p20_MSES_%d"%(i),  guess=22.0,  units = "deg",   description="Angle of Attack + 20deg from MSES, segment %d"%(i)))
    negCMp1.append(        Variable(name="negCMp1_%d"%(i), guess=1.05, units = "-", description="Negative CM about 0.25c + 1, segment %d"%(i)))

objective = W_fuel_out + W_fuel_ret

constraints = []

# =====================================
# SLF
# =====================================
for i in range(0,N_segments):
    constraints += [ W[i] == 0.5 * rho * V[i]**2 * C_L[i] * S ]
    constraints += [ T[i] >= 0.5 * rho * V[i]**2 * C_D[i] * S ]
    constraints += [ Re[i] == rho * V[i] * S**0.5 / (AR**0.5 * mu)]    
# =====================================
# Landing
# =====================================
constraints += [
    W_MTO == 0.5 * rho_SL * V_stall**2 * C_Lmax * S,
    V_stall <= 38*units.m/units.s
] 
# =====================================
# Sprint
# =====================================
constraints += [
    P_max >= T[2] * V[2] / eta_0[2],
    V[2] >= 150*units.m/units.s
] 
# =====================================
# Drag Model
# =====================================
for i in range(0,N_segments):
    constraints += [ C_Dfuse[i] == CDA0/S ]
    constraints += [ C_Di[i] == C_L[i]**2/(np.pi * e * AR) ]
    constraints += [ C_D[i] >= C_Dfuse[i] + C_Dp[i] + C_Di[i] ]
# =====================================
# Propulsive Efficiency
# =====================================
for i in range(0,N_segments):
    constraints += [ eta_0[i] == eta_eng * eta_prop[i]]
    constraints += [ eta_prop[i] == eta_i[i] * eta_v]
    constraints += [ 4*eta_i[i] + T[i]*eta_i[i]**2 / (0.5 * rho * V[i]**2 * A_prop) <= 4]
# =====================================
# Range
# =====================================
constraints += [ R >= 5000 * units.km ]
for i in range(0,N_segments-1):
    constraints += [ z_bre[i] == g * R * T[i] / (h_fuel * eta_0[i] * W[i])]
constraints += [ W_fuel_out/W[0] >= z_bre[0] + z_bre[0]**2/2 + z_bre[0]**3/6 + z_bre[0]**4/24]
constraints += [ W_fuel_ret/W[1] >= z_bre[1] + z_bre[1]**2/2 + z_bre[1]**3/6 + z_bre[1]**4/24]
# =====================================
# Weight
# =====================================
constraints += [ 
    W_pay >= 500*units.kg * g,
    W_tilde >= W_fixed + W_pay + W_eng,
    W_zfw >= W_tilde + W_wing,
    W_eng >= k_ew * P_max**0.803,
    W_wing / f_wadd >= W_web + W_cap,
    W[0] >= W_zfw + W_fuel_ret,
    W_MTO >= W[0] + W_fuel_out,
    W[1] >= W_zfw,
    W[2] == W[0]
]
# =====================================
# Wing Structure
# =====================================
constraints += [ 
    2*q >= 1 + p,
    p >= 1.9,
    M_r_bar == W_tilde * AR * p / 24,
    0.92 * w_bar * tau * t_cap_bar**2 + I_cap_bar <= 0.92**2 / 2 * w_bar * tau**2 * t_cap_bar ,
    8 == N_lift * M_r_bar * AR * q**2 * tau / (S * I_cap_bar * sigma_max),
    12 == AR * W_tilde * N_lift * q**2 / (tau * S * t_web_bar * sigma_max_shear),
    nu**3.94 >= 0.86 * p**(-2.38) + 0.14*p**(0.56),
    W_cap >= 8 * rho_cap * g * w_bar * t_cap_bar * S**1.5 * nu / (3*AR**0.5),
    W_web >= 8 * rho_web * g * r_h * tau * t_web_bar * S**1.5 * nu / (3 * AR**0.5),
    tau <= 0.15,
    q <= 2
]
# =====================================
# MSES Runtime Drag Model
# =====================================
for i in range(0,N_segments):
    fm = fullModel()
    fm.inputMode = 5
    fm.outputMode = 3
    fm.Mach = 0.2
    fm.Coarse_Iteration  =  50 
    fm.Fine_Iteration    =  50
    fm.problemObj = myProblem
    constraints += [RuntimeConstraint([C_Dp[i], C_L_MSES[i], negCMp1[i]],['>=', '==', '>='],[Alpha_p20_MSES[i],Re[i],maxCamber,camberLocation,tau],fm)]
constraints += [ C_L_MSES[i] == C_L[i] for i in range(0,N_segments)]
# Impose moment constraint
constraints += [negCMp1[i] <= 1.06 for i in range(0,N_segments)]

formulation = Formulation(objective, constraints)

bdc = [vr >= 1e-12 * vr.units for vr in formulation.variables_only]
formulation.constraints.extend(bdc)

names = [vr.name for vr in formulation.variables_only]
AR = formulation.variables_only[names.index('AR')]
S = formulation.variables_only[names.index('S')]

phantomConstraints = [
    AR >= 17,
    S <= 40 * units.m**2
]
formulation.constraints.extend(phantomConstraints)

# Get initial guess from previous phase parameters
newX0 = {}
for vname in myProblem.parameter.keys():
    vl = myProblem.parameter[vname].value
    if isinstance(vl, pyCAPS.Quantity):
        vl = vl.value() * getattr(units, str(vl.unit()))
    else:
        vl = vl * units.dimensionless
    newX0[vname] = vl

# Set initial guesses for new variables
newX0['maxCamber']      = 0.02 * units.dimensionless
newX0['camberLocation'] = 0.4 * units.dimensionless
newX0['negCMp1_0']      = 1.05 * units.dimensionless
newX0['negCMp1_1']      = 1.05 * units.dimensionless
newX0['negCMp1_2']      = 1.05 * units.dimensionless

formulation.solverOptions.x0 = newX0
formulation.solverOptions.tau = 0.3
formulation.solverOptions.solver = 'cvxopt'
formulation.solverOptions.solveType = 'slcp'
formulation.solverOptions.relativeTolerance = 1e-5
formulation.solverOptions.baseStepSchedule = [1.0] * 15 + (1/np.linspace(1,100,100)**0.6).tolist()
formulation.solverOptions.progressFilename = None
rs = solve(formulation)

# Save results as CAPS parameters
for vname in rs.variables.keys():
    vl = rs.variables[vname].magnitude
    ut = '{:C}'.format(rs.variables[vname].units)
    capsVar = vl * pyCAPS.Unit(ut)
    if vname in myProblem.parameter:
        myProblem.parameter[vname].value = capsVar
    else:
        myProblem.parameter.create(vname, capsVar)

print(rs.result(10))

# Geometric quantities from optimization result
maxCamber = rs.variables['maxCamber'].to('').magnitude
maxLoc    = rs.variables['camberLocation'].to('').magnitude
thickness = rs.variables['tau'].to('').magnitude
area      = rs.variables['S'].to('m^2').magnitude
aspect    = rs.variables['AR'].to('').magnitude
taper     = rs.variables['q'].to('').magnitude - 1
# Fuel weight and drag from optimization result
fuelOut   = rs.variables['W_fuel_out'].to('N').magnitude
fuelRet   = rs.variables['W_fuel_ret'].to('N').magnitude
fuelTot   = fuelOut + fuelRet
CD0       = rs.variables['C_D_0'].to('').magnitude
CD1       = rs.variables['C_D_1'].to('').magnitude
CD_avg    = (CD0 + CD1) / 2

# Update geometry in ESP
myProblem.geometry.despmtr["camber"   ].value = maxCamber
myProblem.geometry.despmtr["maxloc"   ].value = maxLoc
myProblem.geometry.despmtr["thickness"].value = thickness
myProblem.geometry.despmtr["area"     ].value = area
myProblem.geometry.despmtr["aspect"   ].value = aspect
myProblem.geometry.despmtr["taper"    ].value = taper
myProblem.geometry.cfgpmtr["view:MSES"].value = 0

myProblem.closePhase()