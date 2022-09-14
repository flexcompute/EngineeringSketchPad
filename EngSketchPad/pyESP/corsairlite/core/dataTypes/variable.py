# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
import numpy as np
import corsairlite
from corsairlite import units, Q_, _Unit
# =====================================================================================================================
# Import the class which contains most of the functions for both Variable and Constant (look at this file if you can't
# find the function you are looking for)
# =====================================================================================================================
from corsairlite.core.dataTypes.sharedVariableConstant import SharedVariableConstant
# =====================================================================================================================
# Define the Variable Class
# =====================================================================================================================
class Variable(SharedVariableConstant):
    def __init__(self, name, guess, units, description = '', type_ = float, bounds = None, size = 0, transformation = None):
        # Order matters
        variableUnits = units
        units = corsairlite.units
        super(Variable, self).__init__(name=name, units=variableUnits, description = description, size = size)
        self.guess = guess
        self.type_ = type_
        self.bounds = bounds
        self.transformation = transformation
# =====================================================================================================================
# An initial guess for the variable, also used as normalization
# =====================================================================================================================
    @property
    def guess(self):
        if self._guess is not None:
            return self._guess.to(self.units)
        else:
            return None
    @guess.setter
    def guess(self,val):
        if val is None:
            self._guess = None
        else:
            if isinstance(val, Q_):
                try:
                    val = val.to(self.units).magnitude
                except:
                    raise ValueError('Defined units do not match the units on the input guess.  Please correct.')

            tmpArray = np.array(val)
            if tmpArray.shape == ():
                if self.size == 0:
                    self._guess = val * self.units
                else:
                     raise ValueError('Defined size does not match the size of the input guess.  Please correct.')
            else:
                inputSize = list(tmpArray.shape)
                val = tmpArray.tolist()
                if len(inputSize)==1:
                    inputSize = inputSize[0]
                if inputSize == self.size:
                    self._guess = val * self.units
                else:
                    raise ValueError('Defined size does not match the size of the input guess.  Please correct.')
# =====================================================================================================================
# Type of variable, can be float, int, or bool
# =====================================================================================================================
    @property
    def type_(self):
        return self._type_
    @type_.setter
    def type_(self, val):
        if val==bool or val==int or val==float:
            self._type_ = val
        else:
            raise ValueError('Invalid type_.  Must be a bool, int, or float.')
# =====================================================================================================================
# Bounds on the variable, should be two element iterable
# =====================================================================================================================
    @property
    def bounds(self):
        return self._bounds
    @bounds.setter
    def bounds(self, val):
        if val is None:
            self._bounds = val
        else:
            if not isinstance(val, (list,tuple)):
                raise ValueError('Invalid bound, must be a list or tuple')
            if len(val)!=2:
                raise ValueError('Invalid bound, must have two elements')
            if isinstance(val[0], (int, float)) and isinstance(val[1], (int, float)):
                v0 = val[0]*self.units
                v1 = val[1]*self.units
                self._bounds = [v0,v1]
            elif isinstance(val[0], Q_) and isinstance(val[1], Q_):
                v0 = val[0].to(self.units)
                v1 = val[1].to(self.units)
                self._bounds = [v0,v1]
            else:
                raise ValueError('bounds not assigned')


# =====================================================================================================================
# Optional transformation of the variable, useful for values known to be highly nonlinear (Reynolds, etc...)
# =====================================================================================================================
    @property
    def transformation(self):
        return self._transformation
    @transformation.setter
    def transformation(self, val):
        if val is None:
            self._transformation = val
        else:
            if type(val)==str:
                self._transformation = val
            else:
                raise ValueError('Invalid transformation.  Must be a string.')
# =====================================================================================================================
# Implements the V[0,0] syntax.  Returns scalar from an N dimensional variable
# =====================================================================================================================
    def __getitem__(self, key):
        chk = []
        chk.append(isinstance(key, int))
        chk.append(isinstance(key, tuple))
        if not any(chk):
            raise ValueError('Index must be an integer or tuple.')
            
        if isinstance(self.size, int):
            if isinstance(key,tuple):
                raise ValueError('Too many indicies passed in')
            if self.size == 0:
                raise ValueError('Variable is scalar and be indexed.')

            res = Variable(name = self.name + '_%d'%(key),
                           guess = self.guess.__getitem__(key),
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
            res = Variable(name = self.name + subscript,
                           guess = self.guess.__getitem__(key),
                           units = self.units, 
                           size = 0, 
                           description = self.description + ' [the value at index %s]'%(subscript[1:])
                          )
            return res
# =====================================================================================================================
# Check if two variables are the same
# =====================================================================================================================
    def isEqual(self, other):
        from corsairlite.optimization.constant import Constant      # Must be here to prevent circular import on init
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(other, Variable):
            checks = []
            checks.append(self.name == other.name)
            checks.append(self.units == other.units)
            checks.append(self.size == other.size)
            checks.append(self.description == other.description)
            if all(checks):
                return True
            else:
                return False
        elif isinstance(other, Constant):
            return False
        elif isinstance(other, Term):
            tm = self.toTerm()
            return tm.isEqual(other)
        elif isinstance(other, Expression):
            tm = self.toTerm()
            expr = tm.toExpression()
            return expr.isEqual(other)
        else:
            raise TypeError('Operation between a Variable and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Evaluates the variable to a quantity by taking in a design point
# =====================================================================================================================
    def evaluate(self, x_star):
        return x_star[self.name].to(self.units)
# =====================================================================================================================
# =====================================================================================================================
