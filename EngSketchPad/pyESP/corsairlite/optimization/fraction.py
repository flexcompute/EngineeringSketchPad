# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
from itertools import compress
from corsairlite.core.dataTypes import TypeCheckedList
from corsairlite.optimization.constant import Constant
from corsairlite.optimization import Variable
from corsairlite.optimization.term import Term
from corsairlite.optimization.expression import Expression
from corsairlite import units, Q_, _Unit
from corsairlite.core.pint_custom.errors import DimensionalityError
from corsairlite.core.dataTypes.optimizationObject import OptimizationObject
# =====================================================================================================================
# Define the Fraction Class
# =====================================================================================================================
class Fraction(OptimizationObject):
    def __init__(self):
        self.numerator = Expression()
        self.denominator = Expression()
# =====================================================================================================================
# The printing function
# =====================================================================================================================    
    def __repr__(self):
        pstr = ''
        pstr += '( '
        pstr += self.numerator.__repr__()
        pstr += ' ) / ( '
        pstr += self.denominator.__repr__()
        pstr += ' )'
        return pstr
# =====================================================================================================================
# Flips numerator and denominator
# =====================================================================================================================  
    def flip(self):
        tempNum = copy.deepcopy(self.numerator)
        tempDem = copy.deepcopy(self.denominator)
        self.numerator = tempDem
        self.dem = tempNum
# =====================================================================================================================
# Getter and setter for the numerator
# =====================================================================================================================
    @property
    def numerator(self):
        return self._numerator
    @numerator.setter
    def numerator(self, val):
        if isinstance(val, Expression):
            self._numerator = val
        elif isinstance(val, Term):
            self._numerator = val.toExpression()
        elif isinstance(val, (Variable,Constant)):
            tm = val.toTerm()
            self._numerator = tm.toExpression()
        elif isinstance(val,(Q_,float,int)):
            tm = Term()
            cst = val*units.dimensionless
            tm.constant = cst
            self._numerator = tm.toExpression()
        else:
            raise ValueError('Numerator must be int, float, Quantity, Variable, Constant, Term, or Expression')
# =====================================================================================================================
# Getter and setter for the denominator
# =====================================================================================================================
    @property
    def denominator(self):
        return self._denominator
    @denominator.setter
    def denominator(self, val):
        if isinstance(val, Expression):
            self._denominator = val
        elif isinstance(val, Term):
            self._denominator = val.toExpression()
        elif isinstance(val, (Variable,Constant)):
            tm = val.toTerm()
            self._denominator = tm.toExpression()
        elif isinstance(val,(Q_,float,int)):
            tm = Term()
            cst = val*units.dimensionless
            tm.constant = cst
            self._denominator = tm.toExpression()
        else:
            raise ValueError('Numerator must be int, float, Quantity, Variable, Constant, Term, or Expression')
# =====================================================================================================================
# The variables in the expression
# =====================================================================================================================    
    @property
    def variables(self):
        nm_vars = self.numerator.variables
        dm_vars = self.denominator.variables
        compiledVariables = []
        for v in nm_vars + dm_vars:
            tfl = [v.isEqual(cv) for cv in compiledVariables]
            if not any(tfl):
                compiledVariables.append(v)
        return compiledVariables
# =====================================================================================================================
# The variables in the expression, variables only, no unsubbed constants
# =====================================================================================================================    
    @property
    def variables_only(self):
        nm_vars = self.numerator.variables_only
        dm_vars = self.denominator.variables_only
        compiledVariables = []
        for v in nm_vars + dm_vars:
            tfl = [v.isEqual(cv) for cv in compiledVariables]
            if not any(tfl):
                compiledVariables.append(v)
        return compiledVariables
# =====================================================================================================================
# The variables in the expression, all including exponents
# =====================================================================================================================    
    @property
    def variables_all(self):
        nm_vars = self.numerator.variables_all
        dm_vars = self.denominator.variables_all
        compiledVariables = []
        for v in nm_vars + dm_vars:
            tfl = [v.isEqual(cv) for cv in compiledVariables]
            if not any(tfl):
                compiledVariables.append(v)
        return compiledVariables
# =====================================================================================================================
# Gets and checks the units of the expression.  Returns only a valid unit type, not simplified or human intrepretable
# units.  To accomplish that, utilize the converter command: self.to(conversionUnits)
# =====================================================================================================================
    @property
    def units(self):
        nm_units = self.numerator.units
        dm_units = self.denominator.units
        return nm_units/dm_units
# =====================================================================================================================
# Pint style unit converter
# =====================================================================================================================
    def to(self, conversionUnits):
        try:
            tmp = 1*self.units
            tmp.to(conversionUnits)
        except:
            raise ValueError('Cannot convert from %s to %s'%(self.units, conversionUnits))
        
        self.numerator *= tmp.to(conversionUnits)/(1*self.units)
# =====================================================================================================================
# Pint style unit converter to base units
# =====================================================================================================================
    def to_base_units(self):
        uq = 1.0 * self.units
        uq = uq.to_base_units()
        self.to(uq.units)
# =====================================================================================================================
# Simplifies the numerator and denominator.  Multivariate posynomials are not generally factorable, so 
# no simplification is done between numerator and denominator.  Does check if denominator is a quantity
# or if the fraction is 1
# =====================================================================================================================
    def simplify(self):
        self.numerator.simplify()
        self.denominator.simplify()
        denm = copy.deepcopy(self.denominator)
        try:
            denm.toTerm()
        except:
            pass        
        
        if isinstance(denm, (Q_, Constant, Variable, Term)):
            tms = self.numerator.terms
            self.__class__ = Expression
            self.terms = tms
            self /= denm

        if self.denominator.isEqual(self.numerator):
            self.__class__ = Q_
            self._magnitude = 1.0
            self._units = units.dimensionless
# =====================================================================================================================
# Overload addition
# =====================================================================================================================
    def __add__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            try:
                other.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator + other*slf.denominator
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator + other*slf.denominator
            return slf
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator + other*slf.denominator
            return slf
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator + other*slf.denominator
            return slf
        elif isinstance(other,Fraction):
            othr = copy.deepcopy(other)
            slf = copy.deepcopy(self)
            if othr.denominator.isEqual(slf.denominator):
                newDem = slf.denominator
                newNum = slf.numerator + othr.numerator
            else:
                newDem = slf.denominator * othr.denominator
                newNum = othr.denominator*slf.numerator + slf.denominator*othr.numerator
            newFrac = Fraction()
            newFrac.numerator = newNum
            newFrac.denominator = newDem
            return newFrac
        elif isinstance(other, Constraint):
            othr = copy.deepcopy(other)
            othr.RHS = self + other.RHS
            othr.RHS.simplify()
            othr.LHS = self + other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __radd__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_)):
            res = self.__add__(other)
            return res
        elif isinstance(other, (Variable, Constant, Term, Expression, Constraint)):
            raise RuntimeError('The addition function of the %s object is failing'%(type(other).__name__))
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Overload subtraction
# =====================================================================================================================
    def __sub__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            try:
                other.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator - other*slf.denominator
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator - other*slf.denominator
            return slf
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator - other*slf.denominator
            return slf
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.numerator = slf.numerator - other*slf.denominator
            return slf
        elif isinstance(other,Fraction):
            othr = copy.deepcopy(other)
            slf = copy.deepcopy(self)
            if othr.denominator.isEqual(slf.denominator):
                newDem = slf.denominator
                newNum = slf.numerator - othr.numerator
            else:
                newDem = slf.denominator * othr.denominator
                newNum = othr.denominator*slf.numerator - slf.denominator*othr.numerator
            newFrac = Fraction()
            newFrac.numerator = newNum
            newFrac.denominator = newDem
            return newFrac
        elif isinstance(other, Constraint):
            othr = copy.deepcopy(other)
            othr.RHS = self - other.RHS
            othr.RHS.simplify()
            othr.LHS = self - other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# # ---------------------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        return other + -1.0 * self
# =====================================================================================================================
# Overload multiplication
# =====================================================================================================================
    def __mul__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_, Term, _Unit, Expression)):
            slf = copy.deepcopy(self)
            slf.numerator *= other
            return slf
        elif isinstance(other, (Constant, Variable)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            slf = copy.deepcopy(self)
            slf.numerator *= other
            return slf
        elif isinstance(other,Fraction):
            othr = copy.deepcopy(other)
            slf = copy.deepcopy(self)
            newNum = slf.numerator * othr.numerator
            newDem = slf.denominator * othr.denominator
            newFrac = Fraction()
            newFrac.numerator = newNum
            newFrac.denominator = newDem
            return newFrac
        elif isinstance(other, Constraint):
            return NotImplemented
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rmul__(self, other):
        return self.__mul__(other)
# =====================================================================================================================
# Overload division
# =====================================================================================================================
    def __truediv__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_, Term, _Unit, Expression)):
            slf = copy.deepcopy(self)
            slf.denominator *= other
            return slf
        elif isinstance(other, (Constant, Variable)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            slf = copy.deepcopy(self)
            slf.denominator *= other
            return slf
        elif isinstance(other,Fraction):
            othr = copy.deepcopy(other)
            slf = copy.deepcopy(self)
            newNum = slf.numerator * othr.denominator
            newDem = slf.denominator * othr.numerator
            newFrac = Fraction()
            newFrac.numerator = newNum
            newFrac.denominator = newDem
            return newFrac
        elif isinstance(other, Constraint):
            return NotImplemented
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, other):
        slf = copy.deepcopy(self)
        slf.flip()
        return other * slf
# =====================================================================================================================
# Overload power (**)
# =====================================================================================================================    
    def __pow__(self, other):
        slf = copy.deepcopy(self)
        slf.numerator = slf.numerator ** other
        slf.denominator = slf.denominator ** other
        return slf
# =====================================================================================================================
# Overload the right power, only applicable for ints and floats
# =====================================================================================================================
    def __rpow__(self,other):
        raise TypeError('Fraction is not a valid type in an exponent')
# =====================================================================================================================
# Overload negative sign
# =====================================================================================================================
    def __neg__(self):
        slf = copy.deepcopy(self)
        slf.numerator *= -1
        return slf
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
# Determines if two expressions are equal
# =====================================================================================================================   
    def isEqual(self, other):
        raise ValueError('Lack of general factorability makes this operation invalid.')
#           if isinstance(other, (float,int)):
#             tm = Term()
#             tm.constant = other * units.dimensionless
#             expr = tm.toExpression()
#             frac = Fraction()
#             frac.numerator = expr
#             frac.denominator = 1.0*units.dimensionless
#             return frac.isEqual(self)
#         elif isinstance(other, Q_):
#             tm = Term()
#             tm.constant = other
#             expr = tm.toExpression()
#             return expr.isEqual(self)
#         elif isinstance(other, (Variable, Constant)):
#             tm = other.toTerm()
#             expr = tm.toExpression()
#             return expr.isEqual(self)
#         elif isinstance(other,Term):
#             expr = other.toExpression()
#             return expr.isEqual(self)
#         elif isinstance(other, Expression):
#             slf = copy.deepcopy(self)
#             othr = copy.deepcopy(other)
#             slf.simplify()
#             othr.simplify()
#         elif isinstance(other, Fraction):
#             slf = copy.deepcopy(self)
#             othr = copy.deepcopy(other)
#             slf.simplify()
#             othr.simplify()

#             if len(slf.terms) != len(othr.terms):
#                 return False

#             for tm in slf.terms:
#                 tfl = []
#                 for checkTerm in othr.terms:
#                     tfl.append(tm.isEqual(checkTerm))
#                 if sum(tfl) !=1:
#                     return False

#             return True
#         else:
#             raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================    
# Substitutes in a dict of quantities
# =====================================================================================================================    
    def substitute(self, substitutions = {}, preserveType = False, skipConstants=False):
        numSub = self.numerator.substitute(substitutions = substitutions, preserveType = True, skipConstants=skipConstants)
        demSub = self.denominator.substitute(substitutions = substitutions, preserveType = True, skipConstants=skipConstants)
        
        newFrac = numSub/demSub
        newFrac.simplify()

        if isinstance(newFrac,Expression):
            newFrac = newFrac.substitute(preserveType = preserveType)

        return newFrac
# =====================================================================================================================
# Evaluates the variable to a quantity by taking in a design point
# =====================================================================================================================
    def evaluate(self, x_star):
        res = self.substitute(substitutions = x_star)
        if isinstance(res, Q_):
            return res
        else:
            raise ValueError('Insufficent variables defined in the solution dictionary')
# =====================================================================================================================
# Evaluates the sensitivity of a constant
# =====================================================================================================================
    def sensitivity(self, c, subsDict):
        self.simplify()
        if isinstance(c,Constant):
            val_top  = self.numerator.substitute(subsDict)
            val_bot  = self.denominator.substitute(subsDict)
            sens_top = self.numerator.sensitivity(c,subsDict)
            sens_bot = self.denominator.sensitivity(c,subsDict)
            
            sens_tot = (sens_top * val_bot - sens_bot * val_top) / (val_bot * val_bot)
            return sens_tot
        else:
            raise ValueError('Can only take sensitivity with respect to a constant')
# =====================================================================================================================
# Take derivative with respect to a variable
# =====================================================================================================================
    def derivative(self, x):
        if isinstance(x,Variable):
            self.simplify()
            
            top =  copy.deepcopy(self.numerator)
            bot =  copy.deepcopy(self.denominator)
            dtop = top.derivative(x)
            dbot = bot.derivative(x)
            
            newFrac = (dtop*bot - dbot*top) / (bot*bot)
            newFrac.simplify()
            return newFrac
        
        else:
            raise ValueError('Can only take derivative with respect to a Variable')
# =====================================================================================================================
# =====================================================================================================================