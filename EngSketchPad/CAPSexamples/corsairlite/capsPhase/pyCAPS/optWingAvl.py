import pyCAPS
import numpy as np
import os
from corsairlite.optimization.solve import solve
from corsairlite.optimization import Variable, RuntimeConstraint
from corsairlite.optimization import Formulation, Constant
from corsairlite.analysis.models.model import Model
from corsairlite import units

import argparse

# Setup and read command line options. Please note that this isn't required for pyCAPS
parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)

#Setup the available commandline options
parser.add_argument("-outLevel", default = 0, type=int, choices=[0, 1, 2], help="Set output verbosity")
args = parser.parse_args()

class AVLModel(Model):
    def __init__(self, problemObj):
        super(AVLModel, self).__init__()
        self.inputs.append( Variable("lam", 1.0, "", "Taper ratio"))
        self.inputs.append( Variable("AR", 1.0, "", "Aspect Ratio"))
        self.outputs.append( Variable("e", 1.0, "", "Oswald"))
        self.availableDerivative = 1
        self.problemObj = problemObj

    def BlackBox(*args, **kwargs):
        args = list(args)       # convert tuple to list
        self = args.pop(0)      # pop off the self argument

        # parse inputs to BlackBox from RuntimeConstraint
        runCases, returnMode, extra = self.parseInputs(*args, **kwargs)

        # get the AVL aim
        avl = self.problemObj.analysis["avl"]

        lam = runCases[0]['lam'].to('').magnitude
        AR = runCases[0]['AR'].to('').magnitude
        print("AVL run case - lam, AR:", lam, AR)
        self.problemObj.geometry.despmtr["wing:taper"].value=lam
        self.problemObj.geometry.despmtr["wing:aspect"].value=AR

        outputE = avl.output["e"].value
        y = outputE*units.dimensionless

        dydx = []

        self.problemObj.geometry.despmtr["wing:taper"].value = lam+0.01
        e_pos = avl.output["e"].value

        self.problemObj.geometry.despmtr["wing:taper"].value = lam-0.01
        e_neg= avl.output["e"].value
        de_dlam = (e_pos - e_neg) / 0.02
        dydx.append(de_dlam * units.dimensionless)

        self.problemObj.geometry.despmtr["wing:aspect"].value = AR+0.1
        e_pos = avl.output["e"].value

        self.problemObj.geometry.despmtr["wing:aspect"].value = AR-0.1
        e_neg= avl.output["e"].value
        de_dAR = (e_pos - e_neg) / 0.2
        dydx.append(de_dAR * units.dimensionless)

        print("y:", y)
        print("dydx:", dydx)
        opt = [y, dydx]

        return opt

myProblem = pyCAPS.Problem("ostrich",
                           capsFile=os.path.join("csm","wing2.csm"),
                           phaseName="optWingAvl",
                           phaseStart="sizeWing",
                           phaseContinuation=True,
                           outLevel=args.outLevel)

myProblem.intentPhrase(["optimizing wing"])

# Get/create the avl AIM and set inputs
avl = myProblem.analysis.create(aim  = "avlAIM", name = "avl")

# set the non-geometric AVL inputs

# Cruise CL from optimization of phase sizeWing
avl.input["CL"  ].value = myProblem.parameter["CLcruise"].value

avl.input["Mach"].value = 0
avl.input["Beta"].value = 0
avl.input.AVL_Surface   = {"wing" : {"numChord"     : 4,
                                     "numSpanTotal" : 24}}

# get the number of rocks from the CAPS Problem
nrock = myProblem.parameter.nrock

print("nrock", nrock)

# Constants
rho      = Constant(name="rho",      value=1.225, units="kg/m^3", description="Density of air")
Wrock    = Constant(name="Wrock",    value=2.0,   units="N",      description="Weight of each rock")
Wfixed   = Constant(name="Wfixed",   value=0.1,   units="N",      description="Weight of fixed components")
rho_wing = Constant(name="rho_wing", value=1.0,   units="kg/m^3", description="Weight density of wing")
tau      = Constant(name="tau",      value=0.12,  units="",       description="Airfoil thickness as t/c")
CD0      = Constant(name="CD0",      value=0.04,  units="",       description="Profile drag")
nrock    = Constant(name="nrock",    value=nrock, units="",       description="Number of rocks to carry")
g        = Constant(name="g",        value=9.81,  units="m/s^2",  description="Gravitational constant")
CLmax    = Constant(name="CLmax",    value=1.0,   units="",       description="Max CL")
ARmax    = Constant(name="ARmax",    value=10.0,  units="",       description="Max AR")
E0       = Constant(name="E0",       value=10.0,  units="J",      description="Maximum elastic potential energy in rubber band")

S      = Variable(name="S",     guess=1.0,  units="m^2", description="Wing area")
AR     = Variable(name="AR",    guess=9.0,  units="",    description="Wing aspect ratio")
V      = Variable(name="V",     guess=1.0,  units="m/s", description="Airspeed")
W      = Variable(name="W",     guess=10.0, units="N",   description="Aircraft weight")
Wwing  = Variable(name="Wwing", guess=0.1,  units="N",   description="Wing weight")
CL     = Variable(name="CL",    guess=1.0,  units="",    description="Lift Coefficient")
CD     = Variable(name="CD",    guess=0.1,  units="",    description="Drag Coefficient")
lam    = Variable(name="lam",   guess=0.3,  units="",    description="Wing taper ratio")
e      = Variable(name="e",     guess=1.0,  units="",    description="Oswald efficiency factor")

objective = CD / CL
constraints = []

# min cd / cl occurs at cl = 1.0635, cd/cl = 0.075225
constraints += [
    W     == 0.5 * rho * V**2 * CL * S,
    CL    <= CLmax,
    CD    >= CD0 + CL**2 / (np.pi * AR * e),
    Wwing == S * (S / AR)**0.5 * tau * rho_wing * g,
    W     >= Wfixed + nrock * Wrock + Wwing,
    AR    <= ARmax,
    V     <= (E0 / (0.5 * W/g)) ** 0.5,
]

avlModel = AVLModel(myProblem)
constraints += [RuntimeConstraint([e], ['=='], [lam, AR], avlModel)]

formulation = Formulation(objective, constraints)
formulation.solverOptions.solver = 'cvxopt'
formulation.solverOptions.solveType = 'slcp'
formulation.solverOptions.progressFilename = None
formulation.solverOptions.stepMagnitudeTolerance=1e-7
formulation.solverOptions.debugOutput = True
rs = solve(formulation)

rs_S   = rs.variables['S'].to('m^2').magnitude
rs_AR  = rs.variables['AR'].to('').magnitude
rs_CL  = rs.variables['CL'].to('').magnitude
rs_CD  = rs.variables['CD'].to('').magnitude
rs_W   = rs.variables['W'].to('N').magnitude
rs_lam = rs.variables['lam'].to('').magnitude
rs_e   = rs.variables['e'].to('').magnitude


print(rs.result(3))
print(rs_CL/rs_CD)

# put the best results into the CAPS Problem
myProblem.geometry.despmtr["wing:area"  ].value = rs_S
myProblem.geometry.despmtr["wing:aspect"].value = rs_AR
myProblem.geometry.despmtr["wing:taper" ].value = rs_lam

myProblem.parameter["CLcruise"].value = rs_CL
myProblem.parameter["CLCD"].value = rs_CL/rs_CD
myProblem.parameter["W"].value = rs_W

# Save optimized oswald efficiency as a CAPS value
myProblem.parameter.create("oswald", rs_e)

# Close the phase
myProblem.closePhase()
