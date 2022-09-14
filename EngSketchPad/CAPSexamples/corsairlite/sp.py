from corsairlite import units
from corsairlite.optimization.variable import Variable
from corsairlite.optimization import Formulation, Constant
from corsairlite.optimization.solve import solve

# Variables
x   = Variable(name = 'x', guess = 1, units = '-', description = 'The X variable')
y   = Variable(name = 'y', guess = 1, units = '-', description = 'The Y variable')
z   = Variable(name = 'z', guess = 1, units = '-', description = 'The Z variable')

# Constants
cst = Constant(name = 'c', value = 2.0, units = '-', description = 'A constant')

# Objective function
objective = x*z/y

# Constraints
constraints = [
    cst <= x,
    x <= 3,
    x**2 >= y,
    x/y == z/3,
    y == 4 + x,
]

formulation = Formulation(objective, constraints)
formulation.solverOptions.solveType = 'sp'
formulation.solverOptions.solver = 'cvxopt'
rs = solve(formulation)
print(rs.result(5))