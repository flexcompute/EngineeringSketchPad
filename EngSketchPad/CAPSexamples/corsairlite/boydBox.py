from corsairlite import units
from corsairlite.optimization import Variable, Constant, Formulation
from corsairlite.optimization.solve import solve

# Variables
w   = Variable(name = 'w', guess = 1.0, units = 'm', description = 'The X variable')
h   = Variable(name = 'h', guess = 1.0, units = 'm', description = 'The Y variable')
d   = Variable(name = 'd', guess = 1.0, units = 'm', description = 'The Z variable')

# Constants
Aflr  = Constant(name = 'Aflr',  value = 1000  , units = 'm^2', description = 'Area of floor')
Awall = Constant(name = 'Awall', value =  100  , units = 'm^2', description = 'Area of wall')
alpha = Constant(name = 'alpha', value =    0.5, units = '-',   description = 'Aspect Ratio Constraint')
beta  = Constant(name = 'beta',  value =    2  , units = '-',   description = 'Aspect Ratio Constraint')
gamma = Constant(name = 'gamma', value =    0.5, units = '-',   description = 'Aspect Ratio Constraint')
delta = Constant(name = 'delta', value =    2  , units = '-',   description = 'Aspect Ratio Constraint')

# Objective function: maximize volume
objective = 1/(h*w*d)

# Constraints
constraints = [
    2*(h*w + h*d) <= Awall,
    w*d <= Aflr,
    alpha <= h/w,
    h/w <= beta,
    gamma <= d/w,
    d/w <= delta,
]

formulation = Formulation(objective, constraints)
formulation.solverOptions.solveType = 'gp'
rs = solve(formulation)
print(rs.result(5))