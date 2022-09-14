# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
import numpy as np
import corsairlite
from corsairlite import units, Q_, _Unit
from corsairlite.optimization import Variable
from corsairlite.core.dataTypes.optimizationObject import OptimizationObject
# =====================================================================================================================
# Import the class which contains most of the functions for both Variable and Constant (look at this file if you can't
# find the function you are looking for)
# =====================================================================================================================
from corsairlite.core.dataTypes.sharedVariableConstant import SharedVariableConstant
# =====================================================================================================================
# Define the Constant Class
# =====================================================================================================================
class Constant(SharedVariableConstant, OptimizationObject):
    def __init__(self, name, value, units, description = '', size = 0):
        # Order matters
        constantUnits = units
        units = corsairlite.units
        super(Constant, self).__init__(name=name, units=constantUnits, description = description, size = size)
        self.value = value

# =====================================================================================================================
# The value of the constant
# =====================================================================================================================
    @property
    def value(self):
        if self._value is not None:
            return self._value.to(self.units)
        else:
            return None
    @value.setter
    def value(self,val):
        if val is None:
            self._value = None
        else:
            if isinstance(val, Q_):
                try:
                    val = val.to(self.units).magnitude
                except:
                    raise ValueError('Defined units do not match the units on the input value.  Please correct.')

            tmpArray = np.array(val)
            if tmpArray.shape == ():
                if self.size == 0:
                    self._value = val * self.units
                else:
                     raise ValueError('Defined size does not match the size of the input value.  Please correct.')
            else:
                inputSize = list(tmpArray.shape)
                val = tmpArray.tolist()
                if len(inputSize)==1:
                    inputSize = inputSize[0]
                if inputSize == self.size:
                    self._value = val * self.units
                else:
                    raise ValueError('Defined size does not match the size of the input value.  Please correct.')
# =====================================================================================================================
# Implements the C[0,0] syntax.  Returns scalar from an N dimensional constant
# =====================================================================================================================
    def __getitem__(self, key):
        from corsairlite.optimization.constant import Constant
        chk = []
        chk.append(isinstance(key, int))
        chk.append(isinstance(key, tuple))
        if not any(chk):
            raise ValueError('Index must be an integer or tuple.')
            
        if isinstance(self.size, int):
            if isinstance(key,tuple):
                raise ValueError('Too many indicies passed in')
            if self.size == 0:
                raise ValueError('Constant is scalar and be indexed.')

            res = Constant(name = self.name + '_%d'%(key),
                           value = self.value.__getitem__(key),
                           units = self.units, 
                           size = 0, 
                           description = self.description + ' [the value at index %d]'%(key),
                           )
            return res

        else:
            chk = []
            if isinstance(key, int):
                raise ValueError('Not enough indicies')
            if len(key) < len(self.size):
                raise ValueError('Not enough indicies')
            if len(key) > len(self.size):
                raise ValueError('Too many indicies')

            subscript = '_{'
            for ix in key:
                subscript += '%d'%(ix)
                subscript += ','
            subscript = subscript[0:-1]
            subscript += '}'
            res = Constant(name = self.name + subscript,
                           value = self.value.__getitem__(key),
                           units = self.units, 
                           size = 0, 
                           description = self.description + ' [the value at index %s]'%(subscript[1:])
                          )
            return res
# =====================================================================================================================
# Overload the right power, only applicable for ints and floats
# =====================================================================================================================
    def __rpow__(self,other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float, Q_)):
            if other <=0:
                raise ValueError('Negative number to constant power is not supported, derivative is not defined')
            tm = Term()
            tm.constant = 1.0
            tm.hiddenConstants.append(other)
            tm.hiddenConstantExponents.append(copy.deepcopy(self))
            return tm
# =====================================================================================================================
# Check if two constants are the same
# =====================================================================================================================
    def isEqual(self, other):
        from corsairlite.optimization.constant import Constant
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(other, (int, float)):
            return float(other) * units.dimensionless == self.value
        elif isinstance(other, Q_):
            return other == self.value
        elif isinstance(other, Variable):
            return False
        elif isinstance(other, Constant):
            checks = []
            checks.append(self.name == other.name)
            checks.append(self.value == other.value)
            checks.append(self.units == other.units)
            checks.append(self.size == other.size)
            checks.append(self.description == other.description)
            if all(checks):
                return True
            else:
                return False
        elif isinstance(other, Term):
            tm = self.toTerm()
            return tm.isEqual(other)
        elif isinstance(other, Expression):
            tm = self.toTerm()
            expr = tm.toExpression()
            return expr.isEqual(other)
        else:
            raise TypeError('Operation between a Constant and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# =====================================================================================================================
