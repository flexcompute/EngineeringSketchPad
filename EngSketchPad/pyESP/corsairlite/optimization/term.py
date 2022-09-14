# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
import numpy as np
from itertools import compress
import corsairlite
from corsairlite import units, Q_, _Unit
from corsairlite.core.dataTypes import TypeCheckedList
from corsairlite.optimization import Variable
from corsairlite.optimization.constant import Constant
from corsairlite.core.pint_custom.errors import DimensionalityError
from corsairlite.core.dataTypes.optimizationObject import OptimizationObject
# =====================================================================================================================
# Define the Term Class
# =====================================================================================================================
class Term(OptimizationObject):
    def __init__(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        self.variables = TypeCheckedList((Variable,Constant), [])
        self.constant = 0.0*units.dimensionless
        self.exponents = TypeCheckedList((int, float, Constant, Term, Expression), [])
        self.hiddenConstants = TypeCheckedList((int,float, Q_), [])
        self.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), [])
        self.runtimeExpressions = []
        self.runtimeExpressionExponents = []
# =====================================================================================================================
# The printing function
# =====================================================================================================================
    def __repr__(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        self.simplify()
        pstr = ''
        # Extract constant magnitude
        cst = float(self.constant.magnitude)

        # Only print constant if it is not 1
        if cst != 1.0:
            if cst == -1.0:
                pstr += '-'
            else:
                # Check if integer, print no decimals
                if cst.is_integer():
                    pstr += '%d'%(cst)
                else:  # Print to 3 decimals
                    pstr += '%.3f'%(cst)
                    
        # Print the hidden constants
        for i in range(0,len(self.hiddenConstants)):
            if isinstance(self.hiddenConstants[i],Q_):
                hc = float(self.hiddenConstants[i].magnitude)
            else:
                hc = float(self.hiddenConstants[i])
            if i != 0 or abs(cst) != 1:
                pstr += ' * '
            if hc < 0:
                pstr += '('
            if hc.is_integer():
                pstr += '%d'%(hc)
            else:  # Print to 3 decimals
                pstr += '%.3f'%(hc)
            if hc < 0:
                pstr += ')'
            pstr += '**('
            pstr += str(self.hiddenConstantExponents[i])
            pstr += ')'

        # Print the variables
        for i in range(0,len(self.variables)):
            # Do nothing if exponent is zero (result is 1 with no units, thus it's not relevant)
            if isinstance(self.exponents[i], (int,float)):
                if self.exponents[i] != 0.0:
                    # Appends a multiplication sign if this isn't the first variable or if a constant was present
                    if i != 0 or abs(cst) != 1 or len(self.hiddenConstants)!=0:
                        pstr += ' * '
                    # Append variable name
                    pstr += self.variables[i].name
                    # Append exponent
                    if self.exponents[i] != 1:
                        # Append the ** and exponent value
                        if float(self.exponents[i]).is_integer():
                            # No decimals if integer
                            pstr += '**%d'%(self.exponents[i])
                        else:
                            # Three decimals if float
                            pstr += '**%.3f'%(self.exponents[i])
            elif isinstance(self.exponents[i], (Constant,Term,Expression)):
                if isinstance(self.exponents[i], (Term,Expression)):
                    self.exponents[i].simplify()
                # Appends a multiplication sign if this isn't the first variable or if a constant was present
                if i != 0 or abs(cst) != 1 or len(self.hiddenConstants)!=0:
                    pstr += ' * '
                # Append variable name
                pstr += self.variables[i].name
                pstr += '**(%s)'%(self.exponents[i].__repr__())
            else:
                raise ValueError('Invalid exponent type')
            
        if cst == 1.0 and len(self.variables) == 0 and len(self.exponents) == 0 and len(self.hiddenConstants)==0 and len(self.hiddenConstantExponents)==0:
            pstr = '1'
        if cst == -1.0 and len(self.variables) == 0 and len(self.exponents) == 0 and len(self.hiddenConstants)==0 and len(self.hiddenConstantExponents)==0:
            pstr = '-1'

        return pstr
# =====================================================================================================================
# The variables in the term
# =====================================================================================================================    
    @property
    def variables(self):
        return self._variables
    @variables.setter
    def variables(self, val):
        if isinstance(val, (list, tuple)):
            if all(isinstance(x, (Variable, Constant)) for x in val):
                self._variables = val
            else:
                raise ValueError('Variable list must be a list or tuple.  TypeCheckedList and list types are preferred.')
        else:
            raise ValueError('Variable list must be a list or tuple.  TypeCheckedList and list types are preferred.')
# =====================================================================================================================
# The variables in the term, only variables, not unsubbed constants
# =====================================================================================================================    
    @property
    def variables_only(self):
        vrs = self.variables
        vos = []
        for v in vrs:
            if isinstance(v,Variable):
                vos.append(v)
        return vos
# =====================================================================================================================
# The variables in the term, including the unsubbed constants in the various exponents
# =====================================================================================================================    
    @property
    def variables_all(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        vrs = copy.deepcopy(self.variables)

        exps = [self.hiddenConstantExponents[i] for i in range(0,len(self.hiddenConstantExponents))]
        exps += [self.exponents[i] for i in range(0,len(self.exponents))]
        expVars = []
        for e in exps:
            if isinstance(e, (int,float)):
                pass
            elif isinstance(e, Constant):
                expVars.append(e)
            elif isinstance(e, (Term, Expression)):
                expVars += e.variables_all
            else:
                raise ValueError('Invlid Exponent Type')

        for e in expVars:
            tfl = [vrs[i].isEqual(e) for i in range(0,len(vrs))]
            if not any(tfl):
                vrs.append(e)

        return vrs
# =====================================================================================================================
# The leading constant of the term
# =====================================================================================================================    
    @property
    def constant(self):
        return self._constant
    @constant.setter
    def constant(self, val):
        if type(val)==float or type(val)==int:
            val = float(val) * units.dimensionless
        if isinstance(val, Q_):
            if isinstance(val.magnitude, (float, int)):
                self._constant = val
            else:
                # print(val)
                raise ValueError('Invalid constant.  Must be a scalar pint quantity, float, or int.')
        else:
            raise ValueError('Invalid constant.  Must be a scalar pint quantity, float, or int.')
# =====================================================================================================================
# The exponents of each variable in the term
# =====================================================================================================================    
    @property
    def exponents(self):
        return self._exponents
    @exponents.setter
    def exponents(self, val):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(val, (list, tuple)):
            if all(isinstance(x, (int, float, Constant, Term, Expression)) for x in val):
                for ep in val:
                    if isinstance(val, (Term,Expression)):
                        if len(other.variables_only) != 0:
                            raise TypeError('Exponents cannot contain variables, only fixed values and constants.  The following variables are present: %s'%(str(other.variables)))    
                self._exponents = val
            else:
                raise ValueError('Exponents must be a list or tuple.  TypeCheckedList and list types are preferred.')
        else:
            raise ValueError('Exponents must be a list or tuple.  TypeCheckedList and list types are preferred.')
# =====================================================================================================================
# The units of the term
# =====================================================================================================================    
    @property
    def units(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        self.simplify()

        unt = units.dimensionless
        unt *= self.constant.units
        for i in range(0,len(self.variables)):
            var = self.variables[i]
            if isinstance(self.exponents[i], (int,float)):
                unt *= var.units ** self.exponents[i]
            elif isinstance(self.exponents[i], Constant):
                unt *= var.units ** self.exponents[i].value.magnitude
            elif isinstance(self.exponents[i], (Term, Expression)):
                unt *= var.units ** self.exponents[i].substitute().magnitude
            else:
                raise ValueError('Invlid Exponent Type')

        for i in range(0,len(self.hiddenConstants)):
            if isinstance(self.hiddenConstants[i],Q_):
                c = self.hiddenConstants[i]
                if isinstance(self.hiddenConstantExponents[i], (int,float)):
                    unt *= c.units ** self.hiddenConstantExponents[i]
                elif isinstance(self.hiddenConstantExponents[i], Constant):
                    unt *= c.units ** self.hiddenConstantExponents[i].value.magnitude
                elif isinstance(self.hiddenConstantExponents[i], (Term, Expression)):
                    unt *= c.units ** self.hiddenConstantExponents[i].substitute().magnitude
                else:
                    raise ValueError('Invlid Exponent Type')

        unitNames = list(unt._units.keys())
        unitExps = [unt._units[unitNames[i]] for i in range(0,len(unitNames))]
        changedUnits = False
        for i in range(0,len(unitNames)):
            err = unitExps[i] - float(round(unitExps[i])) 
            if abs(err) <= 1e-3:
                changedUnits = True
                self.constant /= units(unitNames[i])**err
                unt /= units(unitNames[i])**err

        if changedUnits:
            unt = unt.units

        return unt
# =====================================================================================================================
# Pint style unit converter
# =====================================================================================================================
    def to(self, conversionUnits):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(conversionUnits, str):
            if conversionUnits in ['-', '', ' ']:
                conversionUnits = 'dimensionless'
        if isinstance(conversionUnits, str):
            conversionUnits = getattr(units, conversionUnits)
        elif isinstance(conversionUnits, _Unit):
            conversionUnits = conversionUnits
        else:
            # print(self.units)
            # print(conversionUnits)
            raise ValueError('Invalid units.  Must be a string compatible with pint or a unit instance.')

        self.simplify()
        unt = units.dimensionless
        for i in range(0,len(self.variables)):
            var = self.variables[i]
            if isinstance(self.exponents[i], (int,float)):
                unt *= var.units ** self.exponents[i]
            elif isinstance(self.exponents[i], Constant):
                unt *= var.units ** self.exponents[i].value.magnitude
            elif isinstance(self.exponents[i], (Term, Expression)):
                unt *= var.units ** self.exponents[i].substitute().magnitude
            else:
                raise ValueError('Invlid Exponent Type')
        newConstantUnitValue = 1.0 * conversionUnits/unt
        # newConstantUnitValue.ito_reduced_units()
        # print(newConstantUnitValue)
        newConstantUnits = newConstantUnitValue.units
        try:
            self.constant = self.constant.to(newConstantUnits)
        except DimensionalityError as e:
            # This needs to be updated to ignore ANY units that are zero
            # right now it is only valid for if all the units are the same and there are some that are close to zero
            vl1 = copy.deepcopy(self.constant)
            vl2 = 1.0*newConstantUnits
            vl1.ito_base_units()
            vl2.ito_base_units()
            uts1 = vl1.units._units
            uts2 = vl2.units._units
            if set(uts1.keys()) == set(uts2.keys()):
                kys = list(uts1.keys())
                vls1 = np.array([uts1[ky] for ky in kys])
                vls2 = np.array([uts1[ky] for ky in kys])
                err = abs(vls1 - vls2)/vls2
                if all(err <= 1e-5):
                    nwcst = vl1.magnitude * vl2.units
                    nwcst.to(newConstantUnits)
                    self.constant = nwcst
                else:
                    raise ValueError('Could not convert from %s to %s'%(str(vl1.units),str(vl2.units)))
            else:
                # Need to try if actually dimensionless
                vl1 = copy.deepcopy(self.constant)
                vl2 = 1.0*newConstantUnits
                vl1.ito_base_units()
                vl2.ito_base_units()
                uts1 = vl1.units._units
                uts2 = vl2.units._units
                vls1 = np.array([abs(uts1[ky]) for ky in list(uts1.keys())] + [0])
                vls2 = np.array([abs(uts2[ky]) for ky in list(uts2.keys())] + [0])
                if all(vls1 < 1e-5) and all(vls2 < 1e-5): 
                    nwcst = vl1.magnitude * vl2.units
                    nwcst.to(newConstantUnits)
                    self.constant = nwcst
                else:
                    raise ValueError('Could not convert from %s to %s'%(str(vl1.units),str(vl2.units)))
        #     # self.constant = self.constant.to(newConstantUnits)
        # except DimensionalityError as e:
        #     print(self.constant.units.__dict__)
        #     print(newConstantUnits.__dict__)
        #     raise e
# =====================================================================================================================
# Pint style unit converter to base units
# =====================================================================================================================
    def to_base_units(self):
        uq = 1.0 * self.units
        uq = uq.to_base_units()
        self.to(uq.units)
# =====================================================================================================================
# Simplifies the term by combining all like variables 
# =====================================================================================================================    
    def simplify(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init

        simple = True

        # Simplify all the exponents
        for i in range(0,len(self.exponents)):
            ex = copy.deepcopy(self.exponents[i])
            if isinstance(ex,(float, int)):
                if ex == 0:
                    simple = False
            elif isinstance(ex, Constant):
                if ex.value.magnitude == 0.0:
                    simple = False
            elif isinstance(ex, (Term, Expression)):
                ex.simplify()
                if ex.isEqual(0.0):
                    simple = False
            else:
                raise ValueError('Invalid Exponent Type')
            self.exponents[i] = ex
        for i in range(0,len(self.hiddenConstantExponents)):
            ex = copy.deepcopy(self.hiddenConstantExponents[i])
            if isinstance(ex,(float, int)):
                if ex == 0:
                    simple = False
            elif isinstance(ex, Constant):
                if ex.value.magnitude == 0.0:
                    simple = False
            elif isinstance(ex, (Term, Expression)):
                ex.simplify()
                if ex.isEqual(0.0):
                    simple = False
            else:
                raise ValueError('Invalid Exponent Type')
            self.hiddenConstantExponents[i] = ex

         # Multiples of variables or hiddenConstants
        vas = [self.variables, self.hiddenConstants]
        eas = [self.exponents, self.hiddenConstantExponents]

        for iii in [0,1]:
            vbls = vas[iii]
            expnts = eas[iii]
           
            for i in range(0,len(vbls)):
                chks = []
                for j in range(0,len(vbls)):
                    if i != j:
                        if iii == 0:
                            chks.append(vbls[i].isEqual(vbls[j]))
                        else:
                            chks.append(vbls[i] == vbls[j])
                if any(chks):
                    simple = False

            if not simple:
                oldVars = copy.deepcopy(vbls)
                oldExps = copy.deepcopy(expnts)

                if iii == 0:
                    newVars = TypeCheckedList((Variable, Constant), [])
                if iii == 1:
                    newVars = TypeCheckedList((int,float), [])
                newExps = TypeCheckedList((int, float, Constant, Term, Expression), [])

                for i in range(0,len(vbls)):
                    var = oldVars[0]
                    if iii == 0:
                        tfl = [var.isEqual(oldVars[ii]) for ii in range(0,len(oldVars))]
                    else:
                        tfl = [var == oldVars[ii] for ii in range(0,len(oldVars))]
                    exps = list(compress(oldExps, tfl))
                    if all(isinstance(x, (int, float)) for x in exps):
                        ex = sum(exps)
                        if ex != 0:
                            newVars.append(var)
                            newExps.append(ex)
                    elif all(isinstance(x, (int, float, Constant, Term, Expression)) for x in exps):
                        ex = exps[0]
                        for ii in range(1,len(exps)):
                            ex += exps[ii]
                        ex.simplify()
                        if not ex.isEqual(0):
                            newVars.append(var)
                            newExps.append(ex)
                    else:
                        raise ValueError('Invalid Exponent Type')
                    
                    oldVars = list(compress(oldVars, [not ii for ii in tfl]))
                    oldExps = list(compress(oldExps, [not ii for ii in tfl]))

                    if len(oldVars) == 0:
                        break

                if iii == 0:
                    self.variables = newVars
                    self.exponents = newExps
                if iii == 1:
                    self.hiddenConstants = newVars
                    self.hiddenConstantExponents = newExps
            simple = True

        # Check that the hidden constants produce valid numbers
        # Check if zero, check if imaginary, remove any of the magnitude 1s
        eq0  = [0.0  == self.hiddenConstants[i].magnitude for i in range(0,len(self.hiddenConstants)) ]
        eq1  = [1.0  == self.hiddenConstants[i].magnitude for i in range(0,len(self.hiddenConstants)) ]
        # eqn1 = [-1.0 == self.hiddenConstants[i] for i in range(0,len(self.hiddenConstants)) ]
        el0  = [0.0   > self.hiddenConstants[i].magnitude for i in range(0,len(self.hiddenConstants)) ]

        if any(eq0):
            self.variables = TypeCheckedList((Variable,Constant), [])
            self.constant = 0.0*units.dimensionless
            self.exponents = TypeCheckedList((int, float, Constant, Term, Expression), [])
            self.hiddenConstants = TypeCheckedList((int,float,Q_), [])
            self.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), [])
        if any(eq1):
            self.hiddenConstants = TypeCheckedList((int,float,Q_), list(compress(self.hiddenConstants,[not ii for ii in eq1])))
            self.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), list(compress(self.hiddenConstantExponents,[not ii for ii in eq1])))
        if any(el0):
            negExps = TypeCheckedList((int, float, Constant, Term, Expression), list(compress(self.hiddenConstantExponents,el0)))
            for ne in negExps:
                if isinstance(ne, (int,float)):
                    neVal = float(ne)
                elif isinstance(ne, Constant):
                    neVal = float(ne.value.magnitude)
                elif isinstance(ne, (Term, Expression)):
                    neVal = float(ne.substitute().magnitude)
                else:
                    raise ValueError('Invlid Exponent Type')
            if not neVal.is_integer():
                raise ValueError('A negative hidden constant raised to fractional power will produce an imaginary number')
        # if any(eqn1):
        #     neg1Exps = TypeCheckedList((int, float, Constant, Term, Expression), list(compress(self.hiddenConstantExponents,eqn1)))
        #     neValSum = 0.0
        #     for ne in neg1Exps:
        #         if isinstance(ne, (int,float)):
        #             neVal = float(ne)
        #         elif isinstance(ne, Constant):
        #             neVal = float(ne.value.magnitude)
        #         elif isinstance(ne, (Term, Expression)):
        #             neVal = float(ne.substitute().magnitude)
        #         else:
        #             raise ValueError('Invlid Exponent Type')  
        #         neValSum += neVal          
        #     if neValSum%2 == 0:
        #         multiplier = 1.0
        #     else:
        #         multiplier = -1.0
        #     self.hiddenConstants = TypeCheckedList((int,float,Q_), list(compress(self.hiddenConstants,[not ii for ii in eqn1])))
        #     self.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), list(compress(self.hiddenConstantExponents,[not ii for ii in eqn1])))
        #     self.constant *= multiplier

# =====================================================================================================================
# Overload addition
# =====================================================================================================================    
    def __add__(self, other):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint  # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = float(other) * units.dimensionless
        if isinstance(other, Q_):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            term2 = Term()
            term2.constant = other
            res.terms.append(term2)
            res.simplify()
            return res
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on a scalar %s.'%(type(other).__name__.lower()))
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            term2 = Term()
            term2.variables.append(other)
            term2.constant = 1.0*units.dimensionless
            term2.exponents.append(1.0)
            res.terms.append(term2)
            res.simplify()
            return res
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                # print(self)
                # print(other)
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            res.terms.append(copy.deepcopy(other))
            res.simplify()
            return res
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            for tm in other.terms:
                res.terms.append(copy.deepcopy(tm))
            res.simplify()
            return res
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator += self * othr.denominator
            othr.simplify()
            return othr
        elif isinstance(other, Constraint):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            othr.RHS = slf + other.RHS
            othr.RHS.simplify()
            othr.LHS = slf + other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
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
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint  # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = float(other) * units.dimensionless
        if isinstance(other, Q_):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            term2 = Term()
            term2.constant = -1.0 * other
            res.terms.append(term2)
            res.simplify()
            return res
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on a scalar %s.'%(type(other).__name__.lower()))
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            term2 = Term()
            term2.variables.append(other)
            term2.constant = -1.0*units.dimensionless
            term2.exponents.append(1.0)
            res.terms.append(term2)
            res.simplify()
            return res
        elif isinstance(other, Term):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            res = Expression()
            res.terms.append(copy.deepcopy(slf))
            res.terms.append(copy.deepcopy(-1.0*othr))
            res.simplify()
            return res
        elif isinstance(other, Expression):
            try:
                oterm = 1.0*other.units
                oterm.to(self.units)
            except DimensionalityError as e:
                raise e
            res = Expression()
            res.terms.append(copy.deepcopy(self))
            for tm in other.terms:
                res.terms.append(copy.deepcopy(-1.0*tm))
            res.simplify()
            return res
        elif isinstance(other, Fraction):
            othr = copy.deepcopy(other)
            othr.numerator -= self * othr.denominator
            othr.simplify()
            return othr
        elif isinstance(other, Constraint):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            othr.RHS = slf - other.RHS
            othr.RHS.simplify()
            othr.LHS = slf - other.LHS
            othr.LHS.simplify()
            return othr
        else:
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        return other + -1.0 * self
# =====================================================================================================================
# Overload Multiplication
# =====================================================================================================================    
    def __mul__(self, other):
        from corsairlite.optimization.expression import Expression
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.constraint import Constraint
        if isinstance(other, int) or isinstance(other, float):
            other = float(other) * units.dimensionless
        if isinstance(other, (Q_, _Unit)):
            slf = copy.deepcopy(self)
            slf.constant *= other
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            slf = copy.deepcopy(self)
            slf.variables.append(other)
            slf.exponents.append(1.0)
            slf.simplify()
            return slf
        elif isinstance(other, Term):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            if isinstance(slf.variables, tuple):
                slf.variables = list(slf.variables)
            if isinstance(othr.variables, tuple):
                othr.variables = list(othr.variables)
            slf.variables = slf.variables + othr.variables
            if isinstance(slf.exponents, tuple):
                slf.exponents = list(slf.exponents)
            if isinstance(othr.exponents, tuple):
                othr.exponents = list(othr.exponents)
            slf.exponents = slf.exponents + othr.exponents
            if isinstance(slf.hiddenConstants, tuple):
                slf.hiddenConstants = list(slf.hiddenConstants)
            if isinstance(othr.hiddenConstants, tuple):
                othr.hiddenConstants = list(othr.hiddenConstants)
            slf.hiddenConstants = slf.hiddenConstants + othr.hiddenConstants
            if isinstance(slf.hiddenConstantExponents, tuple):
                slf.hiddenConstantExponents = list(slf.hiddenConstantExponents)
            if isinstance(othr.hiddenConstantExponents, tuple):
                othr.hiddenConstantExponents = list(othr.hiddenConstantExponents)
            slf.hiddenConstantExponents = slf.hiddenConstantExponents + othr.hiddenConstantExponents
            slf.constant *= othr.constant
            slf.simplify()
            return slf
        elif isinstance(other, Expression):
            res = Expression()
            for term in other.terms:
                x = self * term
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
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
# --------------------------------------------------------------------------------------------------------------------- 
    def __rmul__(self, other):
        return self.__mul__(other)
# =====================================================================================================================
# Overload Division
# =====================================================================================================================    
    def __truediv__(self, other):
        from corsairlite.optimization.fraction import Fraction      # Must be here to prevent circular import on init
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(other, int) or isinstance(other, float):
            other = float(other) * units.dimensionless
        if isinstance(other, (Q_, _Unit)):
            slf = copy.deepcopy(self)
            slf.constant /= other
            return slf
        elif isinstance(other, (Variable, Constant)):
            if other.size != 0:
                raise ValueError('Operations can only be performed on scalar variables.')
            slf = copy.deepcopy(self)
            slf.variables.append(other)
            slf.exponents.append(-1.0)
            slf.simplify()
            return slf
        elif isinstance(other, Term):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            if isinstance(slf.variables, tuple):
                slf.variables = list(slf.variables)
            if isinstance(othr.variables, tuple):
                othr.variables = list(othr.variables)
            slf.variables = slf.variables + othr.variables
            if isinstance(slf.exponents, tuple):
                slf.exponents = list(slf.exponents)
            othrExps = [-1.0 * othr.exponents[i] for i in range(0,len(othr.exponents))]
            slf.exponents = slf.exponents + othrExps
            if isinstance(slf.hiddenConstants, tuple):
                slf.hiddenConstants = list(slf.hiddenConstants)
            if isinstance(othr.hiddenConstants, tuple):
                othr.hiddenConstants = list(othr.hiddenConstants)
            slf.hiddenConstants = slf.hiddenConstants + othr.hiddenConstants
            if isinstance(slf.hiddenConstantExponents, tuple):
                slf.hiddenConstantExponents = list(slf.hiddenConstantExponents)
            othrHCEs = [-1.0 * othr.hiddenConstantExponents[i] for i in range(0,len(othr.hiddenConstantExponents))]
            slf.hiddenConstantExponents = slf.hiddenConstantExponents + othrHCEs
            slf.constant = self.constant/other.constant
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
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, other):
        slf = copy.deepcopy(self)
        newExps = [-1.0 * slf.exponents[i] for i in range(0,len(slf.exponents))]
        newHCEs = [-1.0 * slf.hiddenConstantExponents[i] for i in range(0,len(slf.hiddenConstantExponents))]
        slf.exponents = newExps
        slf.hiddenConstantExponents = newHCEs
        slf.constant = 1.0/slf.constant
        return other * slf
# =====================================================================================================================
# Overload power (**)
# =====================================================================================================================    
    def __pow__(self, other):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(other, (Term, Expression)):
            if len(other.variables_only) != 0:
                raise TypeError('Exponents cannot contain variables, only fixed values and constants.  The following variables are present: %s'%(str(other.variables))) 

        if isinstance(other, (float, int)):
            slf = copy.deepcopy(self)
            slf.constant = slf.constant ** other
            slf.exponents = TypeCheckedList((int, float, Constant, Term, Expression), [other * ex for ex in slf.exponents])
            slf.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), [other * ex for ex in slf.hiddenConstantExponents])
            return slf
        elif isinstance(other, (Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            slf.exponents = TypeCheckedList((int, float, Constant, Term, Expression), [other * ex for ex in slf.exponents])
            slf.hiddenConstantExponents = TypeCheckedList((int, float, Constant, Term, Expression), [other * ex for ex in slf.hiddenConstantExponents])
            nhc = self.constant.magnitude
            slf.hiddenConstants = TypeCheckedList((int,float,Q_), [nhc] + slf.hiddenConstants)
            slf.hiddenConstantExponents += [other] + slf.hiddenConstantExponents
            slf.constant = 1.0*units.dimensionless 
            return slf
        else:
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
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
        self.constant *= -1.0
        return self
# =====================================================================================================================
# Determines if two terms are similar (share everything but leading constant)
# =====================================================================================================================    
    def similarTerms(self, other):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        if isinstance(other,Term):
            self.simplify()
            other.simplify()

            # Create two dictionaries with keys of variable names and values the variables themselves.  Check if equal
            dct1_vars = {}
            dct2_vars = {}
            for i in range(0,len(self.variables)):
                dct1_vars[self.variables[i].name] = self.variables[i]
            for i in range(0,len(other.variables)):
                dct2_vars[other.variables[i].name] = other.variables[i]
            try:
                varsSame1 = [dct1_vars[ky].isEqual(dct2_vars[ky]) for ky in list(dct1_vars.keys())]
                varsSame2 = [dct2_vars[ky].isEqual(dct1_vars[ky]) for ky in list(dct2_vars.keys())]
            except KeyError as e:
                return False
            
            # Create two dictionaries with keys of variable names and values of exponents.  Check if equal
            dct1_exps = {}
            dct2_exps = {}
            expChecks = []
            for i in range(0,len(self.variables)):
                dct1_exps[self.variables[i].name] = self.exponents[i]
            for i in range(0,len(other.variables)):
                dct2_exps[other.variables[i].name] = other.exponents[i]
            for ky,vl in dct1_exps.items():
                if isinstance(vl,(int,float)):
                    if isinstance(dct2_exps[ky],(int,float)):
                        expChecks.append(vl == dct2_exps[ky])
                    else:
                        return False
                else:
                    expChecks.append(vl.isEqual(dct2_exps[ky]))
                    
            variableChecks = all(expChecks) and all(varsSame1) and all(varsSame2)
            if not variableChecks:
                return variableChecks

            # Check the hidden constants
            hc1 = [self.hiddenConstants[i].to_base_units().magnitude for i in range(0,len(self.hiddenConstants))]
            hc2 = [other.hiddenConstants[i].to_base_units().magnitude for i in range(0,len(other.hiddenConstants))]

            ix1 = np.array(hc1).argsort()
            ix2 = np.array(hc2).argsort()

            e1 = copy.deepcopy(self.hiddenConstantExponents)
            e2 = copy.deepcopy(other.hiddenConstantExponents)
            e1sorted = [e1[ixe1] for ixe1 in ix1]
            e2sorted = [e2[ixe2] for ixe2 in ix2]
            hc1sorted = [hc1[ixe1] for ixe1 in ix1]
            hc2sorted = [hc2[ixe2] for ixe2 in ix2]
            hConstantsSame = hc1sorted==hc2sorted

            if not hConstantsSame:
                return False

            chkList = []
            for i in range(0,len(e1sorted)):
                ev1 = e1sorted[i]
                ev2 = e2sorted[i]
                vls = [0, 1]
                for ii in [0,1]:
                    if ii == 0:
                        vl = ev1
                    else:
                        vl = ev2
                    if isinstance(vl,(int, float)):
                        vls[ii] = float(vl)
                    elif isinstance(vl,Constant):
                        vls[ii] = vl
                    elif isinstance(vl,(Term, Expression)):
                        vl.simplify()
                        vls[ii] = vl
                    else:
                        raise ValueError('Invalid Exponent Type')
                if type(vls[0]) != type(vls[1]):
                    chkList.append(False)
                else:
                    if isinstance(vls[0],(int,float)):
                        chkList.append(vls[0]==vls[1])
                    elif isinstance(vls[0], (Constant,Term,Expression)):
                        chkList.append(vls[0].isEqual(vls[1]))
                    else:
                        raise ValueError('Invalid Exponent Type')
            if sum(chkList)<len(chkList):
                return False
            else:
                return True
        else:
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Determines if two terms are equal
# =====================================================================================================================   
    def isEqual(self, other, tol = 1e-6):
        from corsairlite.optimization.expression import Expression
        from corsairlite.optimization.constraint import Constraint
        self.simplify()
        if isinstance(other, (float,int)):
            tm = Term()
            tm.constant = other * units.dimensionless
            return tm.isEqual(self)
        elif isinstance(other, Q_):
            tm = Term()
            tm.constant = other
            return tm.isEqual(self)
        elif isinstance(other, (Variable, Constant)):
            othr = other.toTerm()
            return self.isEqual(othr)
        elif isinstance(other,Term):
            other.simplify()
            othr = other.substitute(preserveType = True)
            slf = self.substitute(preserveType = True)
            sim = slf.similarTerms(othr)
            oc = copy.deepcopy(othr.constant).to_base_units()
            sc = copy.deepcopy(slf.constant).to_base_units()
            sameMagnitudes = abs(oc.magnitude - sc.magnitude) <= tol
            sameUnits = oc.units == sc.units
            return sim and sameMagnitudes and sameUnits
        elif isinstance(other, Expression):
            other.simplify()
            slf = self.toExpression()
            return other.isEqual(slf)
        else:
            raise TypeError('Operation between a Term and an object of type %s is not defined.'%(type(other).__name__))
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
# Converts the term to a single term expression
# =====================================================================================================================    
    def toExpression(self):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        tmp = Expression()
        tmp.terms.append(self)
        return tmp
# =====================================================================================================================    
# Substitutes in a dict of quantities
# =====================================================================================================================    
    def substitute(self, substitutions = {}, preserveType = False, skipConstants = False):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        res = Term()
        res.constant = 1.0*units.dimensionless
        tm = copy.deepcopy(self)
        cst = tm.constant
        for i in range(0,len(tm.hiddenConstants)):
            hc = copy.deepcopy(tm.hiddenConstants[i])
            hce = copy.deepcopy(tm.hiddenConstantExponents[i])
            if isinstance(hce, (int,float)):
                cst *= hc ** hce
            elif isinstance(hce, Constant):
                if not skipConstants or hce.name in list(substitutions.keys()):
                    try:
                        tmpvar = 1.0 * hce.units
                        tmpvar.to('dimensionless')
                    except:
                        raise ValueError('Exponents must be dimensionless')
                    cst *= hc ** hce.value
                else:
                    res *= hc ** hce
            elif isinstance(hce, (Term, Expression)):
                try:
                    hce.to('dimensionless')
                except:
                    raise ValueError('Exponents must be dimensionless')
                cst *= hc ** hce.substitute(substitutions, skipConstants = skipConstants)
            else:
                raise ValueError('Invalid exponent type')

        for i in range(0,len(tm.variables)):
            v = copy.deepcopy(tm.variables[i])
            expnt = copy.deepcopy(tm.exponents[i])
            if isinstance(expnt, (int, float)):
                pass # do nothing
            elif isinstance(expnt, Constant):
                if not skipConstants or expnt.name in list(substitutions.keys()):
                    expnt = expnt.value
                    try:
                        expnt.to('dimensionless')
                    except:
                        raise ValueError('Exponents must be dimensionless')
                    expnt = expnt.to('dimensionless').magnitude
                else:
                    pass
            elif isinstance(expnt, (Term, Expression)):
                expnt = expnt.substitute(substitutions, skipConstants = skipConstants)
                try:
                    expnt.to('dimensionless')
                except:
                    raise ValueError('Exponents must be dimensionless')
                if isinstance(expnt,(int,float,Q_)):
                    expnt = expnt.to('dimensionless').magnitude
            else:
                raise ValueError('Invalid exponent type')

            if v.name in list(substitutions.keys()):
                sbVal = substitutions[v.name]
                if isinstance(sbVal,(int, float, Q_)):
                    vl = substitutions[v.name].to(v.units)**expnt
                elif isinstance(sbVal,(Constant,Variable)):
                    sbVal = sbVal.toTerm()
                    sbVal.to(v.units)
                    vl = sbVal**expnt
                elif isinstance(sbVal,(Term)):
                    sbVal.to(v.units)
                    vl = sbVal**expnt
                else:
                    raise ValueError('Substitute value must be int, float, quantity, constant, variabe, term')
                if isinstance(vl,(int, float, Q_)):
                    cst *= vl
                else:
                    res *= vl
            elif isinstance(v,Constant):
                if not skipConstants or v.name in list(substitutions.keys()):
                    cst *= v.value**expnt
                else:
                    res *= v**expnt
            else:
                res.variables.append(v)
                res.exponents.append(expnt)
        res.constant *= cst

        res.simplify()
        untDict = res.constant.units.__dict__['_units']._d
        untKeys = untDict.keys()
        delks = []
        for k in untKeys:
            if abs(untDict[k]) < 1e-5:
                delks.append(k)
        for k in delks:
            del untDict[k]

        if not preserveType:
            if res.variables == [] and res.hiddenConstants == []:
                c1clean = res.constant.magnitude * units.dimensionless
                for ky in list(res.constant.units._units.keys()):
                    vl = res.constant.units._units[ky]
                    if abs(vl - round(vl)) < 1e-3:
                        c1clean *= units(ky) ** round(vl)
                    else:
                        c1clean *= units(ky) ** vl
                rv = c1clean
            else:
                rv = res
        else:
            rv = res
        return rv
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
# Find the constant sensitivity
# =====================================================================================================================
    def sensitivity(self, x, subsDict):
        from corsairlite.optimization.expression import Expression  # Must be here to prevent circular import on init
        # print(self)
        indexSelected = False
        self.simplify()
        if isinstance(x,Constant):
            for i in range(0,len(self.hiddenConstants)):
                ex = self.hiddenConstantExponents[i]
                if isinstance(ex, Constant):
                    tf = x.isEqual(ex)
                    if tf:
                        indexSelected = True
                        break
                elif isinstance(ex,(Term,Expression)):
                    tfl = [x.isEqual(v) for v in ex.variables_all]
                    if any(tfl):
                        indexSelected = True
                        break
                else:
                    raise ValueError('Invalid Exponent')
            
            if indexSelected:
                if isinstance(ex,Constant):
                    everythingElse = copy.deepcopy(self)
                    del everythingElse.hiddenConstants[i]
                    del everythingElse.hiddenConstantExponents[i]
                    sens = self.hiddenConstants[i]**(ex.value.magnitude) * np.log(self.hiddenConstants[i])
                    sens *= everythingElse.substitute(subsDict) 
                    sens += self.hiddenConstants[i]**(ex.value.magnitude) * everythingElse.sensitivity(x)
                    return sens
                else: # is Term or Expression
                    everythingElse = copy.deepcopy(self)
                    del everythingElse.hiddenConstants[i]
                    del everythingElse.hiddenConstantExponents[i]

                    sens = self.hiddenConstants[i]**(ex.substitute().magnitude) * np.log(self.hiddenConstants[i]) 
                    sens *= ex.sensitivity(x,subsDict) # sensitivity of exponent
                    sens *= everythingElse.substitute(subsDict) 
                    sens += self.hiddenConstants[i]**(ex.substitute().magnitude) * everythingElse.sensitivity(x, subsDict)
                    return sens
            else:
                for j in range(0,len(self.variables)):
                    v  = self.variables[j]
                    ex = self.exponents[j]
                    if x.isEqual(v):
                        indexSelected = True
                        break
                    if isinstance(ex,(int,float)):
                        pass
                    elif isinstance(ex, Constant):
                        tf = x.isEqual(ex)
                        if tf:
                            indexSelected = True
                            break
                    elif isinstance(ex,(Term,Expression)):
                        tfl = [x.isEqual(v) for v in ex.variables_all]
                        if any(tfl):
                            indexSelected = True
                            break
                    else:
                        raise ValueError('Invalid Exponent')
                        
            if indexSelected:
                if isinstance(ex,(int, float)):
                    everythingElse = copy.deepcopy(self)
                    del everythingElse.variables[j]
                    del everythingElse.exponents[j]
                    
                    if isinstance(v,Variable):
                        vVal = subsDict[v.name]
                    else:
                        vVal = v.value
                    s1 = ex * vVal**(ex-1)*everythingElse.substitute(subsDict)
                    s1clean = s1.magnitude * units.dimensionless
                    for ky in list(s1.units._units.keys()):
                        vl = s1.units._units[ky]
                        if abs(vl - round(vl)) < 1e-3:
                            s1clean *= units(ky) ** round(vl)
                        else:
                            s1clean *= units(ky) ** vl
                    s2 = vVal**ex * everythingElse.sensitivity(x,subsDict)
                    s2clean = s2.magnitude * units.dimensionless
                    for ky in list(s2.units._units.keys()):
                        vl = s2.units._units[ky]
                        if abs(vl - round(vl)) < 1e-3:
                            s2clean *= units(ky) ** round(vl)
                        else:
                            s2clean *= units(ky) ** vl
                    sens = s1clean + s2clean
                    return sens
                elif isinstance(ex,Constant):
                    if x.isEqual(v):
                        everythingElse = copy.deepcopy(self)
                        del everythingElse.variables[j]
                        del everythingElse.exponents[j]
                        
                        vVal = ex.value.magnitude
                        sens = vVal**(ex.value.magnitude) * (np.log(vVal) + 1.0)
                        sens *= everythingElse.substitute(subsDict) 
                        sens += vVal**(ex.value.magnitude) * everythingElse.sensitivity(x,subsDict)      
                        return sens
                    else:
                        everythingElse = copy.deepcopy(self)
                        del everythingElse.variables[j]
                        del everythingElse.exponents[j]
                        
                        if isinstance(v,Variable):
                            vVal = subsDict[v.name]
                        else:
                            vVal = v.value
                        sens = vVal**(ex.value.magnitude) * np.log(vVal.magnitude)
                        sens *= everythingElse.substitute(subsDict) 
                        sens += vVal**(ex.value.magnitude) * everythingElse.sensitivity(x,subsDict)
                        return sens
                else: # is term or expression
                    if x.isEqual(v):
                        everythingElse = copy.deepcopy(self)
                        del everythingElse.variables[j]
                        del everythingElse.exponents[j]

                        if isinstance(v,Variable):
                            vVal = subsDict[v.name]
                        else:
                            vVal = v.value
                        sens = vVal**(ex.substitute().magnitude) 
                        sens *= (np.log(vVal.magnitude) * ex.sensitivity(x,subsDict) + 1/vVal.magnitude * ex.substitute().magnitude) 
                        sens *= everythingElse.substitute(subsDict)
                        sens += vVal**(ex.substitute().magnitude) * everythingElse.sensitivity(x, subsDict)
                        return sens
                    else:
                        everythingElse = copy.deepcopy(self)
                        del everythingElse.variables[j]
                        del everythingElse.exponents[j]

                        if isinstance(v,Variable):
                            vVal = subsDict[v.name]
                        else:
                            vVal = v.value
                        sens = vVal**(ex.substitute().magnitude) * np.log(vVal.magnitude)
                        sens *= ex.sensitivity(x,subsDict) # sensitivity of exponent
                        sens *= everythingElse.substitute(subsDict) 
                        sens += vVal**(ex.substitute().magnitude) * everythingElse.sensitivity(x, subsDict)
                        return sens
            else:
                return 0.0*self.units/x.units
        else:
            raise ValueError('Can only take sensitivity with respect to a constant')
# =====================================================================================================================
# Take derivative with respect to a variable
# =====================================================================================================================
    def derivative(self, x):
        if isinstance(x,Variable):
            self.simplify()
            tfl = [self.variables[i].isEqual(x) for i in range(0,len(self.variables))]
            if any(tfl):
                ix = tfl.index(True)
                ex = self.exponents[ix]
                newTerm = copy.deepcopy(self)
                newTerm.exponents[ix] = newTerm.exponents[ix] - 1
                newTerm *= ex
                newTerm.simplify()
                return newTerm
            else:
                return 0.0 * self.units/x.units
        else:
            raise ValueError('Can only take derivative with respect to a Variable')
# =====================================================================================================================
# =====================================================================================================================    


 