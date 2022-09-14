import dill
import numpy as np
pi = np.pi
from corsairlite import units
from corsairlite.optimization import Variable
from corsairlite.optimization import Formulation, Constant
from corsairlite.core.data.standardAtmosphere import atm
from corsairlite.core.pickleFunctions import setPickle, unsetPickle
from corsairlite.optimization.solve import solve

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
r_h             = Constant(name = "r_h",             value = 0.75,                  units = "-",           description = "Ratio of height at rear spar to maximum wing height" )
rho             = Constant(name = "rho",             value = atm.rho(3000*units.m), units = "kg/m^3",      description = "density of air")
rho_cap         = Constant(name = "rho_cap",         value = 2700,                  units = "kg/m^3",      description = "Density of wing cap material (aluminum)") 
rho_SL          = Constant(name = "rho_SL",          value = atm.rho(0*units.m),    units = "kg/m^3",      description = "density of air, sea level")
rho_web         = Constant(name = "rho_web",         value = 2700,                  units = "kg/m^3",      description = "Density of wing web material (aluminum)") 
sigma_max       = Constant(name = "sigma_max",       value = 310,                   units = "MPa",         description = "Allowable tensile stress of aluminum" )
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

# Objective function
objective = W_fuel_out + W_fuel_ret

# Constraints
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
for i in range(N_segments):
    constraints += [ C_Dfuse[i] == CDA0/S ]
    constraints += [ C_Di[i] == C_L[i]**2/(np.pi * e * AR) ]
    constraints += [ C_D[i] >= C_Dfuse[i] + C_Dp[i] + C_Di[i] ]
    constraints += [ 1 >= (2.56    * C_L[i]**( 5.88) * tau**(-3.32) * Re[i]**(-1.54) * C_Dp[i]**(-2.26) + 
                            3.80e-9 * C_L[i]**(-0.92) * tau**( 6.23) * Re[i]**(-1.38) * C_Dp[i]**(-9.57) + 
                            2.20e-3 * C_L[i]**(-0.01) * tau**( 0.03) * Re[i]**( 0.14) * C_Dp[i]**(-0.73) + 
                            1.19e4  * C_L[i]**( 9.78) * tau**( 1.76) * Re[i]**(-1.00) * C_Dp[i]**(-0.91) + 
                            6.14e-6 * C_L[i]**( 6.53) * tau**(-0.52) * Re[i]**(-0.99) * C_Dp[i]**(-5.19) ) ]
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
    tau <= 0.15,
    M_r_bar == W_tilde * AR * p / 24,
    0.92 * w_bar * tau * t_cap_bar**2 + I_cap_bar <= 0.92**2 / 2 * w_bar * tau**2 * t_cap_bar ,
    8 == N_lift * M_r_bar * AR * q**2 * tau / (S * I_cap_bar * sigma_max),
    12 == AR * W_tilde * N_lift * q**2 / (tau * S * t_web_bar * sigma_max_shear),
    nu**3.94 >= 0.86 * p**(-2.38) + 0.14*p**(0.56),
    W_cap >= 8 * rho_cap * g * w_bar * t_cap_bar * S**1.5 * nu / (3*AR**0.5),
    W_web >= 8 * rho_web * g * r_h * tau * t_web_bar * S**1.5 * nu / (3 * AR**0.5)
]    

formulation = Formulation(objective, constraints)
bdc = [vr >= 1e-12 * vr.units for vr in formulation.variables_only]
formulation.constraints.extend(bdc)

formulation.solverOptions.solver = 'cvxopt'
formulation.solverOptions.solveType = 'gp'
rs = solve(formulation)
print(rs.result(5))

# Dump results in pickle file to use as initial guess for lsqp
setPickle()
pklfl = open('hoburg.pl','wb')
dill.dump(rs, pklfl)
pklfl.close()
unsetPickle()