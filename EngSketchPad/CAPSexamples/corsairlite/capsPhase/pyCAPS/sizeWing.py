import pyCAPS
import numpy as np
import os
from corsairlite.optimization.solve import solve
from corsairlite.optimization import Variable
from corsairlite.optimization import Formulation, Constant
from corsairlite import units

import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument("-nrock", default = 2, type=int, help="Number of rock payloads")
parser.add_argument("-outLevel", default = 0, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

myProblem = pyCAPS.Problem("ostrich",
                           capsFile=os.path.join("csm","wing1.csm"),
                           phaseName="sizeWing",
                           phaseContinuation=True,
                           outLevel=args.outLevel)

myProblem.intentPhrase(["initial sizing of rectangular wing"])

# get the number of rocks from the CAPS Problem
myProblem.parameter.create("nrock", args.nrock)
nrock = myProblem.parameter.nrock

print("nrock", nrock)

# Constants
rho      = Constant(name="rho",      value=1.225, units="kg/m^3", description="Density of air")
Wrock    = Constant(name="Wrock",    value=2.0,   units="N",      description="Weight of each rock")
Wfixed   = Constant(name="Wfixed",   value=0.1,   units="N",      description="Weight of fixed components")
rho_wing = Constant(name="rho_wing", value=1.0,   units="kg/m^3", description="Weight density of wing")
tau      = Constant(name="tau",      value=0.12,  units="",       description="Airfoil thickness as t/c")
CD0      = Constant(name="CD0",      value=0.04,  units="",       description="Profile drag")
e        = Constant(name="e",        value=0.90,  units="",       description="Oswawld efficiency factor")
nrock    = Constant(name="nrock",    value=nrock, units="",       description="Number of rocks to carry")
g        = Constant(name="g",        value=9.81,  units="m/s^2",  description="Gravitational constant")
CLmax    = Constant(name="CLmax",    value=1.0,   units="",       description="Max CL")
ARmax    = Constant(name="ARmax",    value=10.0,  units="",       description="Max AR")
E0       = Constant(name="E0",       value=10.0,  units="J",      description="Maximum elastic potential energy in rubber band")

# Variables
S      = Variable(name="S",     guess=1.0,  units="m^2", description="Wing area")
AR     = Variable(name="AR",    guess=10.0, units="",    description="Wing aspect ratio")
V      = Variable(name="V",     guess=1.0,  units="m/s", description="Airspeed")
V0     = Variable(name="V0",    guess=1.0,  units="m/s", description="Initial airspeed assuming rubber band launch")
W      = Variable(name="W",     guess=1.0,  units="N",   description="Aircraft weight")
Wwing  = Variable(name="Wwing", guess=0.1,  units="N",   description="Wing weight")
CL     = Variable(name="CL",    guess=1.0,  units="",    description="Lift Coefficient")
CD     = Variable(name="CD",    guess=1.0,  units="",    description="Drag Coefficient")

objective = CD / CL
constraints = []

# min cd / cl occurs at cl = 1.0635, cd/cl = 0.075225 in this formulation
constraints += [
    W     == 0.5 * rho * V**2 * CL * S,
    CL    <= CLmax,
    CD    >= CD0 + CL**2 / (np.pi * AR * e),
    Wwing == S * (S / AR)**0.5 * tau * rho_wing * g,
    W     >= Wfixed + nrock * Wrock + Wwing,
    AR    <= ARmax,
    V     <= (E0 / (0.5 * W/g)) ** 0.5
]

formulation = Formulation(objective, constraints)
formulation.solverOptions.solver = 'cvxopt'
formulation.solverOptions.solveType = 'gp'
rs = solve(formulation)

rs_S  = rs.variables['S'].to('m^2').magnitude
rs_AR = rs.variables['AR'].to('').magnitude
rs_V  = rs.variables['V'].to('m/s').magnitude
rs_W  = rs.variables['W'].to('N').magnitude
rs_CL = rs.variables['CL'].to('').magnitude
rs_CD = rs.variables['CD'].to('').magnitude

print(rs.result(3))
print("Optimized L/D:", rs_CL/rs_CD)

# put the best results into the CAPS Problem
myProblem.geometry.despmtr["wing:area"  ].value = rs_S
myProblem.geometry.despmtr["wing:aspect"].value = rs_AR

# Save optimal values as CAPS values.
myProblem.parameter.create("CLcruise", rs_CL)
myProblem.parameter.create("CLCD", rs_CL/rs_CD)
myProblem.parameter.create("W", rs_W)

# Close the phase
myProblem.closePhase()
