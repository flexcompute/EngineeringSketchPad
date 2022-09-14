from corsairlite import units
from corsairlite.optimization import Variable, Formulation, Constant
from corsairlite.optimization.solve import solve

# Variables
x   = Variable(name = 'x', guess = 1.0, units = '-', description = 'The X variable')
y   = Variable(name = 'y', guess = 2.0, units = '-', description = 'The Y variable')
z   = Variable(name = 'z', guess = 3.0, units = '-', description = 'The Z variable')

# Constants
cst = Constant(name = 'c', value = 2.0, units = '-', description = 'A constant')

# Objective function
objective = x*z/y

# Constraints
constraints = [
    2.1 <= x,
    x <= 3,
    x**2 >= y**0.5,
    x/y == z**2
]

formulation = Formulation(objective, constraints)
formulation.solverOptions.solveType = 'gp'
rs = solve(formulation)
print(rs.result(5))
