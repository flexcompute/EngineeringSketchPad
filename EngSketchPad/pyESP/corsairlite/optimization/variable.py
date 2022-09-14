import corsairlite
from corsairlite.core.dataTypes.variable import Variable as VD
from corsairlite.core.dataTypes.optimizationObject import OptimizationObject

class Variable(VD,OptimizationObject):
    def __init__(self, name, guess, units, description = '', type_ = float, bounds = None, size = 0, transformation = None):
        # Order matters
        variableUnits = units
        units = corsairlite.units
        super(Variable, self).__init__(name=name, guess=guess, units=variableUnits, description = description, 
                                       type_ = type_, bounds=bounds, size = size, transformation=transformation)
        self.guess = guess
        self.type_ = type_
        self.bounds = bounds
        self.transformation = transformation