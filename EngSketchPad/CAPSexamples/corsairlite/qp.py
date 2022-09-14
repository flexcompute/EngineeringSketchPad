from corsairlite import units
from corsairlite.optimization import Variable, Formulation, Constant
from corsairlite.optimization.solve import solve

# Constants
cst = Constant(name = 'c', value = 3.0, units = '-', description = 'A constant')
c2  = Constant(name = 'c2', value = 2, units = '-', description = 'A constant') 

# Variables
x   = Variable(name = 'x', guess = 1.0, units = '-', description = 'The X variable')
y   = Variable(name = 'y', guess = 2.0, units = '-', description = 'The Y variable')

# Objective function
obj = x**2 + 4*y**c2 - 8*x - 16*y +1 + x*y

# Constraints
cons = [
    x + y <= 5,
    x <= 3,
    x >= 0,
    y >= 0,
    x + y == cst,
]

formulation = Formulation(obj, cons)
formulation.solverOptions.solveType = 'qp'
formulation.solverOptions.solver = 'cvxopt'
rs = solve(formulation)
rs.computeSensitivities()
print(rs.result(5))