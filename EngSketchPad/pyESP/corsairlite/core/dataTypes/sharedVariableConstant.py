# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
import numpy as np
import corsairlite
from corsairlite import units, Q_, _Unit
# =====================================================================================================================
# Define the Constant Class
# =====================================================================================================================
class SharedVariableConstant(object):
    def __init__(self, name, units, description = '', size = 0):
        # Order matters
        objectUnits = units
        units = corsairlite.units
        self.name = name
        self.units = objectUnits
        self.size = size
        self.description = description
# =====================================================================================================================
# The printing function
# =====================================================================================================================
    def __repr__(self):
        return self.name    
# =====================================================================================================================
# Define the constant name
# =====================================================================================================================
    @property
    def name(self):
        return self._name
    @name.setter
    def name(self,val):
        if isinstance(val, str):
            self._name = val
        else:
            raise ValueError('Invalid name.  Must be a string.')
# =====================================================================================================================
# Define the constant units
# =====================================================================================================================
    @property
    def units(self):
        return self._units
    @units.setter
    def units(self, val):
        # set dimensionless if a null string is passed in
        if isinstance(val, str):
            if val in ['-', '', ' ']:
                val = 'dimensionless'
        
        if isinstance(val, str):
            self._units = getattr(corsairlite.units, val)
        elif isinstance(val, _Unit):
            self._units = val
        else:
            raise ValueError('Invalid units.  Must be a string compatible with pint or a unit instance.')
# =====================================================================================================================
# Define the constant size
# =====================================================================================================================
    @property
    def size(self):
        return self._size
    @size.setter
    def size(self, val):
        invalid = False
        if isinstance(val,(list, tuple)):
            for x in val:
                if not isinstance(x,int):
                    raise ValueError('Invalid size.  Must be an integer or list/tuple of integers')
                if x == 1:
                    raise ValueError('A value of 1 is not valid for defining size.  Use fewer dimensions to represent your constant.')
            self._size = val
        else:
            if isinstance(val, int):
                if val == 1:
                    raise ValueError('A value of 1 is not valid for defining size.  Use 0 to indicate a scalar value.')
                else:
                    self._size = val
            else:
                raise ValueError('Invalid size.  Must be an integer or list/tuple of integers')
# =====================================================================================================================
# Define the constant description
# =====================================================================================================================
    @property
    def description(self):
        return self._description
    @description.setter
    def description(self, val):
        if isinstance(val,str):
            self._description = val
        else:
            raise ValueError('Invalid description.  Must be a string.')
# =====================================================================================================================
# Overload Addition
# =====================================================================================================================
    def __add__(self, other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint  # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float)):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            term1 = self.toTerm()
            term2 = Term()
            term2.constant = other
            res = term1 + term2
            res.simplify()
            return res
        elif isinstance(other, SharedVariableConstant):
            if other.size != 0:
                raise ValueError('Operation can only be performed with a scalar %s.'%(type(other).__name__))
            term1 = self.toTerm()
            term2 = Term()
            term2.variables.append(other)
            term2.constant = 1.0*units.dimensionless
            term2.exponents.append(1.0)
            res = term1 + term2
            res.simplify()
            return res
        elif isinstance(other, Term):
            term1 = self.toTerm()
            term2 = copy.deepcopy(other)
            res = term1 + term2
            res.simplify()
            return res
        elif isinstance(other, Expression):
            res = self.toTerm()
            for tm in other.terms:
                res += tm
            res.simplify()
            return res
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator += self * othr.denominator
            othr.simplify()
            return othr
        elif isinstance(other, Constraint):
            othr = copy.deepcopy(other)
            othr.RHS = self + other.RHS
            othr.RHS.simplify()
            othr.LHS = self + other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between a %s and an object of type %s is not defined.'%(type(self).__name__,type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __radd__(self, other):
        res = self.__add__(other)
        tms = copy.deepcopy(res.terms)
        newTerms = [tms[-1]] + tms[0:-1]
        res.terms = newTerms
        return res
# =====================================================================================================================
# Overload Subtraction
# =====================================================================================================================
    def __sub__(self, other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint  # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float)):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            term1 = self.toTerm()
            term2 = Term()
            term2.constant = other
            res = term1 - term2
            res.simplify()
            return res
        elif isinstance(other, SharedVariableConstant):
            if other.size != 0:
                raise ValueError('Operation can only be performed with a scalar %s.'%(type(other).__name__))
            term1 = self.toTerm()
            term2 = Term()
            term2.variables.append(other)
            term2.constant = 1.0*units.dimensionless
            term2.exponents.append(1.0)
            res = term1 - term2
            res.simplify()
            return res
        elif isinstance(other, Term):
            term1 = self.toTerm()
            term2 = copy.deepcopy(other)
            res = term1 - term2
            res.simplify()
            return res
        elif isinstance(other, Expression):
            res = self.toTerm()
            for tm in other.terms:
                res -= tm
            res.simplify()
            return res
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator -= self * othr.denominator
            othr.simplify()
            return othr
        elif isinstance(other, Constraint):
            othr = copy.deepcopy(other)
            othr.RHS = self - other.RHS
            othr.RHS.simplify()
            othr.LHS = self - other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between a %s and an object of type %s is not defined.'%(type(self).__name__,type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        expr = other + -1.0 * self
        newTerms = copy.deepcopy(expr.terms)
        expr.terms = newTerms # radd automatically switches 0 and 1
        return expr
# =====================================================================================================================
# Overload Multiplication
# =====================================================================================================================
    def __mul__(self, other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint  # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float)):
            other = other * units.dimensionless
        if isinstance(other, (Q_, _Unit)):
            res = self.toTerm()
            res.constant *= other
            return res
        elif isinstance(other, SharedVariableConstant):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar constants.')
            res = self.toTerm()
            res.variables.append(other)
            res.exponents.append(1.0)
            res.constant = 1.0*units.dimensionless
            res.simplify()
            return res
        elif isinstance(other, Term):
            res = self.toTerm()
            othr = copy.deepcopy(other)
            res.constant = othr.constant
            for i in range(0,len(othr.variables)):
                res.variables.append(othr.variables[i])
                res.exponents.append(othr.exponents[i])
            res.simplify()
            return res
        elif isinstance(other, Expression):
            res = Expression()
            slf = copy.deepcopy(self)
            for term in other.terms:
                x = slf * term
                res.terms.append(x)
            res.simplify()
            return res
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator *= self
            othr.simplify()
            return othr
        elif isinstance(other, Constraint):
            return NotImplemented
        else:
            raise TypeError('Operation between a %s and an object of type %s is not defined.'%(type(self).__name__,type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rmul__(self, other):
        return self.__mul__(other)
# =====================================================================================================================
# Overload Division
# =====================================================================================================================
    def __truediv__(self, other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float)):
            other = other * units.dimensionless
        if isinstance(other, (Q_, _Unit)):
            res = self.toTerm()
            res.constant /= other
            return res
        elif isinstance(other, SharedVariableConstant):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar constants.')
            res = self.toTerm()
            res.variables.append(other)
            res.exponents.append(-1.0)
            res.simplify()
            return res
        elif isinstance(other, Term):
            res = self.toTerm()
            for i in range(0,len(other.variables)):
                res.variables.append(other.variables[i])
                res.exponents.append(-1.0*other.exponents[i])
                res.constant /= other.constant
            return res
        elif isinstance(other, Expression):
            try:
                tm = other.toTerm()
                slf = copy.deepcopy(self)
                return slf/tm
            except:
                slf = copy.deepcopy(self)
                othr = copy.deepcopy(other)
                frac = Fraction()
                frac.numerator = slf
                frac.denominator = othr
                frac.simplify()
                return frac
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.flip()
            othr *= self
            othr.simplify()
            return othr
        else:
            raise TypeError('Operation between a %s and an object of type %s is not defined.'%(type(self).__name__,type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, other):
        term = self.toTerm(exponent = -1.0)
        return other * term
# =====================================================================================================================
# Overload Power (**)
# =====================================================================================================================
    def __pow__(self, other):
        from corsairlite.optimization.constant import Constant # Must be here to prevent circular import on init
        from corsairlite.optimization.term import Term # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        if isinstance(other, (int, float, Constant)):
            res = self.toTerm(exponent = other)
            return res
        elif isinstance(other, (Term, Expression)):
            if len(other.variables_only) == 0:
                res = self.toTerm(exponent = other)
                return res
            else:
                raise TypeError('Exponents cannot contain variables, only fixed values and constants.  The following variables are present: %s'%(str(other.variables)))    
        else:
            raise TypeError('Operation between a %s and an object of type %s is not defined.'%(type(self).__name__,type(other).__name__))
# =====================================================================================================================
# Overload negative sign
# =====================================================================================================================
    def __neg__(self):
        res = self.toTerm()
        res.constant = -1.0
        return res
# =====================================================================================================================
# Overload Constraint Creators
# =====================================================================================================================
    def generateConstraint(self, other, operator):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        return Constraint(LHS=self, RHS=other, operator=operator)
# ---------------------------------------------------------------------------------------------------------------------
    def __eq__(self, other):
        return self.generateConstraint(other, '==')
# ---------------------------------------------------------------------------------------------------------------------
    def __le__(self, other):
        return self.generateConstraint(other, '<=')
# ---------------------------------------------------------------------------------------------------------------------
    def __ge__(self, other):
        return self.generateConstraint(other, '>=')
# =====================================================================================================================
# Convert to a Term object
# =====================================================================================================================
    def toTerm(self, exponent = 1.0):
        from corsairlite.optimization.constant import Constant              # Must be here to prevent circular import on init
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression              # Must be here to prevent circular import on init
        if self.size != 0:
            raise ValueError('Operations can only be performed on scalar constants.')
        res = Term()
        res.variables.append(copy.deepcopy(self))
        res.constant = 1.0 * units.dimensionless
        if isinstance(exponent, (int, float, Constant)):
            res.exponents.append(exponent)
        elif isinstance(exponent, (Term, Expression)):
            if len(exponent.variables_only) == 0:
                res.exponents.append(exponent)
            else:
                raise TypeError('Exponents cannot contain variables, only fixed values and constants.  The following variables are present: %s'%(str(other.variables)))   
        else:
            raise ValueError('Invalid exponent')
        return res
# =====================================================================================================================
# =====================================================================================================================
