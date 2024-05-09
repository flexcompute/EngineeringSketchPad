# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
from itertools import compress
from corsairlite.core.dataTypes import TypeCheckedList
from corsairlite.optimization.constant import Constant
from corsairlite.optimization import Variable
from corsairlite.optimization.term import Term
from corsairlite import units, Q_, _Unit
from corsairlite.core.pint_custom.errors import DimensionalityError
from corsairlite.core.dataTypes.optimizationObject import OptimizationObject
# =====================================================================================================================
# Define the Expression Class
# =====================================================================================================================
class Expression(OptimizationObject):
    def __init__(self):
        self.terms = TypeCheckedList(Term, [])
# =====================================================================================================================
# The printing function
# =====================================================================================================================    
    def __repr__(self):
        pstr = ''
        if len(self.terms) == 1:
            self.simplify()
            pstr += self.terms[0].__repr__()
        elif len(self.terms) > 1:
            self.simplify()
            pstr += self.terms[0].__repr__()
            for i in range(1,len(self.terms)):
                tm = self.terms[i]
                if tm.constant < 0:
                    pstr += '  -  '
                    tmCopy = copy.deepcopy(tm)
                    tmCopy.constant = -1*tmCopy.constant
                    pstr += tmCopy.__repr__()
                else:
                    pstr += '  +  '
                    pstr += tm.__repr__()
        else:
            pstr += 'None'
        
        return pstr
# =====================================================================================================================
# The variables in the expression
# =====================================================================================================================    
    @property
    def variables(self):
        placeholdingTerm = Term()
        for i in range(0,len(self.terms)):
            tm = self.terms[i]
            placeholdingTerm.variables += tm.variables
            placeholdingTerm.exponents = [1 for ii in range(0,len(placeholdingTerm.variables))]
            placeholdingTerm.simplify()

        return placeholdingTerm.variables
# =====================================================================================================================
# The variables in the expression, variables only, no unsubbed constants
# =====================================================================================================================    
    @property
    def variables_only(self):
        placeholdingTerm = Term()
        for i in range(0,len(self.terms)):
            tm = self.terms[i]
            placeholdingTerm.variables += tm.variables_only
            placeholdingTerm.exponents = [1 for ii in range(0,len(placeholdingTerm.variables))]
            placeholdingTerm.simplify()

        return placeholdingTerm.variables
# =====================================================================================================================
# The variables in the expression, all including exponents
# =====================================================================================================================    
    @property
    def variables_all(self):
        placeholdingTerm = Term()
        for i in range(0,len(self.terms)):
            tm = self.terms[i]
            placeholdingTerm.variables += tm.variables_all
            placeholdingTerm.exponents = [1 for ii in range(0,len(placeholdingTerm.variables))]
            placeholdingTerm.simplify()

        return placeholdingTerm.variables
# =====================================================================================================================
# Gets and checks the units of the expression.  Returns only a valid unit type, not simplified or human intrepretable
# units.  To accomplish that, utilize the converter command: self.to(conversionUnits)
# =====================================================================================================================
    @property
    def units(self):
        termsCopy = copy.deepcopy(self.terms)
        uts = [term.units for term in self.terms]
        if len(uts) == 0:
            raise ValueError('No terms in the expression')
        if len(uts) == 1:
            return uts[0]
        if len(uts) > 1:
            rtnUnits = uts[0]
            for i in range(1,len(uts)):        
                try:
                    tm = 1.0*uts[i]
                    tm.to(rtnUnits)
                except DimensionalityError as e:
                    raise e
            # Returns units of first term, does no checking or simplification
            return rtnUnits
# =====================================================================================================================
# Pint style unit converter
# =====================================================================================================================
    def to(self, conversionUnits):
        for i in range(0,len(self.terms)):
            self.terms[i].to(conversionUnits)
# =====================================================================================================================
# Pint style unit converter to base units
# =====================================================================================================================
    def to_base_units(self):
        uq = 1.0 * self.units
        uq = uq.to_base_units()
        self.to(uq.units)
# =====================================================================================================================
# Simplifies the expression by combining like terms
# =====================================================================================================================
    def simplify(self):
        for tm in self.terms:
            try:
                tm.simplify()
            except:
                print(tm.__dict__)

        zeroConstants = [tm.constant.magnitude == 0.0 for tm in self.terms]
        if sum(zeroConstants) == len(zeroConstants):
            newTerms = TypeCheckedList(Term, [])
            tm = Term()
            tm.constant = 0.0 * self.units
            newTerms.append(tm)
            self.terms = newTerms
        elif sum(zeroConstants)>=1 and sum(zeroConstants)<len(zeroConstants):
            for i in sorted(range(0,len(self.terms)), reverse=True):
                if zeroConstants[i]:
                    del self.terms[i]
        else:
            # No zero terms
            pass

        simple = True
        for i in range(0,len(self.terms)):
            chks = []
            for j in range(0,len(self.terms)):
                if i != j:
                    if self.terms[i].similarTerms(self.terms[j]):
                        chks.append(True)
                        break
                    else:
                        chks.append(False)
            if any(chks):
                simple = False

        if not simple:
            oldTerms = self.terms
        
            newTerms = TypeCheckedList(Term, [])
            
            maxIters = len(oldTerms)
            for i in range(0,maxIters):
                tm = oldTerms[0]
                # Find terms which share all but constant
                tfl = [tm.similarTerms(oldTerms[ii]) for ii in range(0,len(oldTerms))]
                # Extract similar terms
                similarTerms = list(compress(oldTerms, tfl))
                
                newTm = copy.deepcopy(tm)
                newTm_units = newTm.units
                newTm_stripped = copy.deepcopy(newTm)
                newTm_stripped.constant = 1.0*units.dimensionless
                strippedUnits = newTm_stripped.units
                constantUnits = newTm.units/strippedUnits
                # Sum the constants of all the similar terms
                for smt in similarTerms:
                    smt.to(newTm_units)
                ncArr = [similarTerms[ii].constant for ii in range(0,len(similarTerms))]
                for ii in range(0,len(ncArr)):
                    if ii == 0:
                        newCst = ncArr[0].magnitude * constantUnits
                    else:
                        newCst += ncArr[ii].magnitude * constantUnits
                # Set the constant for the new term
                newTm.constant = newCst
                
                # Append to expression if the new term has a non-zero constant
                if newCst != 0:
                    newTerms.append(newTm)
                    
                # Strip out all the terms that were considered from the oldTerms list
                oldTerms = list(compress(oldTerms, [not ii for ii in tfl]))
                    
                # No terms remain
                if len(oldTerms) == 0:
                    break
                    
            # Terms all canceled out, none remain
            if len(newTerms) == 0:
                tm = Term()
                tm.constant = 0.0 * self.units
                newTerms.append(tm)
                    
            self.terms = newTerms
# =====================================================================================================================
# Getter and setter for the terms in the expression
# =====================================================================================================================
    @property
    def terms(self):
        return self._terms
    @terms.setter
    def terms(self, val):
        if isinstance(val, (list, tuple)):
            if all(isinstance(x, Term) for x in val):
                self._terms = list(val)
        else:
            raise ValueError('Term list must be a list or tuple.  TypeCheckedList and list types are preferred.')
# =====================================================================================================================
# Overload addition
# =====================================================================================================================
    def __add__(self, other):
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            try:
                other.to(self.units)
            except DimensionalityError as e:
                raise e
            cst = Term()
            cst.constant = other.to(self.units)
            slf = copy.deepcopy(self)
            slf.terms.append(cst)
            slf.simplify()
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            newTerm = other.toTerm()
            slf = copy.deepcopy(self)
            slf.terms.append(newTerm)
            slf.simplify()
            return slf
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.terms.append(other)
            slf.simplify()
            return slf
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.terms += other.terms
            slf.simplify()
            return slf
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
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __radd__(self, other):
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_)):
            res = self.__add__(other)
            tms = copy.deepcopy(res.terms)
            newTerms = [tms[-1]] + tms[0:-1]
            res.terms = newTerms
            return res
        elif isinstance(other, (Variable, Constant, Term, Expression, Constraint)):
            raise RuntimeError('The addition function of the %s object is failing'%(type(other).__name__))
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Overload subtraction
# =====================================================================================================================
    def __sub__(self, other):
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = other * units.dimensionless
        if isinstance(other, Q_):
            try:
                other.to(self.units)
            except DimensionalityError as e:
                raise e
            cst = Term()
            cst.constant = -1.0*other.to(self.units)
            slf = copy.deepcopy(self)
            slf.terms.append(cst)
            slf.simplify()
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            newTerm = other.toTerm()
            newTerm.constant = -1.0*units.dimensionless
            slf = copy.deepcopy(self)
            slf.terms.append(newTerm)
            slf.simplify()
            return slf
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            slf.terms.append(-1.0*other)
            slf.simplify()
            return slf
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            ot = copy.deepcopy(other.terms)
            negOT = [-1.0 * tm for tm in ot]
            slf.terms += negOT
            slf.simplify()
            return slf
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
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        return other + -1.0 * self
# =====================================================================================================================
# Overload multiplication
# =====================================================================================================================
    def __mul__(self, other):
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_, Term, _Unit)):
            slf = copy.deepcopy(self)
            newTerms = []
            for tm in slf.terms:
                if isinstance(other, int):
                    tm *= float(other)
                    newTerms.append(tm)
                else:
                    tm *= other
                    newTerms.append(tm)
            slf.terms = newTerms
            slf.simplify()
            return slf
        elif isinstance(other, (Constant, Variable)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            newTerms = []
            slf = copy.deepcopy(self)
            for tm in slf.terms:
                tm *= other
                newTerms.append(tm)
            slf.terms = newTerms
            slf.simplify()
            return slf
        elif isinstance(other, Expression):
            newExpression = Expression()
            for i in range(0,len(self.terms)):
                for j in range(0,len(other.terms)):
                    newExpression.terms.append(self.terms[i]*other.terms[j])
            newExpression.simplify()
            return newExpression
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator *= self
            othr.simplify()
            return othr
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
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_, _Unit)):
            slf = copy.deepcopy(self)
            return 1.0/other * slf
        elif isinstance(other, (Constant, Variable)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            newTerms = []
            slf = copy.deepcopy(self)
            for tm in slf.terms:
                tm /= other
                newTerms.append(tm)
            slf.terms = newTerms
            slf.simplify()
            return slf
        elif isinstance(other, Term):
            slf = copy.deepcopy(self)
            newTerms = []
            for tm in slf.terms:
                tm /= other
                newTerms.append(tm)
            slf.terms = newTerms
            slf.simplify()
            return slf
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
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, other):
        raise ValueError('Division by an expression is not defined')
# =====================================================================================================================
# Overload power (**)
# =====================================================================================================================    
    def __pow__(self, other):
        # Can be valid if integer power
        # Can be approximated if fractional power
        raise NotImplementedError('Power has not been defined as fractional exponents are not readily expandable.  Contact the developers.')
# =====================================================================================================================
# Overload the right power, only applicable for ints and floats
# =====================================================================================================================
    def __rpow__(self,other):
        from corsairlite.optimization.term import Term              # Must be here to prevent circular import on init
        if isinstance(other, (int, float, Q_)):
            if other <=0:
                raise ValueError('Negative number to constant power is not supported, derivative is not defined')
            if len(self.variables_only) != 0:
                raise TypeError('Exponents cannot contain variables, only fixed values and constants.  The following variables are present: %s'%(str(other.variables)))    
            tm = Term()
            tm.constant = 1.0
            tm.hiddenConstants.append(other)
            tm.hiddenConstantExponents.append(copy.deepcopy(self))
            return tm
# =====================================================================================================================
# Overload negative sign
# =====================================================================================================================
    def __neg__(self):
        tms = copy.deepcopy(self.terms)
        newTerms = []
        for tm in tms:
            newTerms.append(-1.0 * tm)
        self.terms = newTerms 
        return self
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
        if isinstance(other, (float,int)):
            tm = Term()
            tm.constant = other * units.dimensionless
            expr = tm.toExpression()
            return expr.isEqual(self)
        elif isinstance(other, Q_):
            tm = Term()
            tm.constant = other
            expr = tm.toExpression()
            return expr.isEqual(self)
        elif isinstance(other, (Variable, Constant)):
            tm = other.toTerm()
            expr = tm.toExpression()
            return expr.isEqual(self)
        elif isinstance(other,Term):
            expr = other.toExpression()
            return expr.isEqual(self)
        elif isinstance(other, Expression):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            slf.simplify()
            othr.simplify()

            if len(slf.terms) != len(othr.terms):
                return False

            for tm in slf.terms:
                tfl = []
                for checkTerm in othr.terms:
                    tfl.append(tm.isEqual(checkTerm))
                if sum(tfl) !=1:
                    return False

            return True
        else:
            raise TypeError('Operation between an Expression and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================    
# Substitutes in a dict of quantities
# =====================================================================================================================    
    def substitute(self, substitutions = {}, preserveType = False, skipConstants=False):
        newExpr = Term()
        newExpr.constant = 0.0 * self.units
        newExpr.toExpression()
        for tm in self.terms:
            subbedTerm = tm.substitute(substitutions = substitutions, skipConstants = skipConstants)
            if isinstance(subbedTerm, Q_):
                s1clean = subbedTerm.magnitude * units.dimensionless
                for ky in list(subbedTerm.units._units.keys()):
                    vl = subbedTerm.units._units[ky]
                    if abs(vl - round(vl)) < 1e-3:
                        s1clean *= units(ky) ** round(vl)
                    else:
                        s1clean *= units(ky) ** vl
                newExpr += s1clean
            else:
                newExpr += subbedTerm

        newExpr.simplify()
        if not preserveType:
            try:
                newExpr = newExpr.toTerm()
                if newExpr.variables == []:
                    rv = newExpr.constant
                else:
                    rv = newExpr
            except:
                rv = newExpr
        else:
            rv = newExpr
        return rv
# =====================================================================================================================    
# Demotes expression to a term if possible
# =====================================================================================================================    
    def toTerm(self):
        self.simplify()

        if len(self.terms) == 1:
            return self.terms[0]
        else:
            raise ValueError('Could not demote the following expression to a term: %s'%(str(self)))
# =====================================================================================================================
# Evaluates the variable to a quantity by taking in a design point
# =====================================================================================================================
    def evaluate(self, x_star):
        res = Q_(0, self.units)
        for term in self.terms:
            res += term.evaluate(x_star)
        return res
# =====================================================================================================================
# Evaluates the sensitivity of a constant
# =====================================================================================================================
    def sensitivity(self, x, subsDict):
        sensList = []
        self.simplify()
        if isinstance(x,Constant):
            for tm in self.terms:
                sensList.append(tm.sensitivity(x,subsDict))
            return sum(sensList)
        else:
            raise ValueError('Can only take sensitivity with respect to a constant')
# =====================================================================================================================
# Take derivative with respect to a variable
# =====================================================================================================================
    def derivative(self, x):
        if isinstance(x,Variable):
            self.simplify()
            newExpr = 0.0 * self.units/x.units
            for i in range(0,len(self.terms)):
                newExpr += self.terms[i].derivative(x)
            if isinstance(newExpr, (Term, Expression)):
                newExpr.simplify()
            return newExpr
        else:
            raise ValueError('Can only take derivative with respect to a Variable')
# =====================================================================================================================
# =====================================================================================================================