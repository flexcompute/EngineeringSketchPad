# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
from corsairlite import units, Q_, _Unit
from corsairlite.optimization.constant import Constant
from corsairlite.optimization import Variable
from corsairlite.optimization.term import Term
from corsairlite.optimization.expression import Expression
from corsairlite.core.pint_custom.errors import DimensionalityError
# =====================================================================================================================
# Define the Constraint Class
# =====================================================================================================================
class Constraint(object):
    def __init__(self, LHS = None, operator = '==', RHS = None):
        self.LHS = LHS
        self.RHS = RHS
        self.operator = operator
# =====================================================================================================================
# The printing function
# =====================================================================================================================
    def __repr__(self):
        pstr = ''
        pstr += self.LHS.__repr__()
        pstr += '  '
        pstr += self.operator
        pstr += '  '
        pstr += self.RHS.__repr__()
        return pstr
# =====================================================================================================================
# Gets and checks the units of the constraint.  Returns only a valid unit type, not simplified or human intrepretable
# units.  To accomplish that, utilize the converter command: self.to(conversionUnits)
# =====================================================================================================================
    @property
    def units(self):
        # Check that the two sides have the same units (may be different systems)
        try:
            RHScp = copy.deepcopy(self.RHS)
            RHScp.to(self.LHS.units)
        except DimensionalityError as e:
            raise e
        # Defaults to the LHS units
        return self.LHS.units
# =====================================================================================================================
# Pint style unit converter
# =====================================================================================================================
    def to(self, conversionUnits):
        self.RHS.to(conversionUnits)
        self.LHS.to(conversionUnits)
# =====================================================================================================================
# Pint style unit converter to base units
# =====================================================================================================================
    def to_base_units(self):
        self.RHS.to_base_units()
        self.LHS.to_base_units()
# =====================================================================================================================
# The left hand side of the constraint getter and setter
# =====================================================================================================================
    @property
    def LHS(self):
        return self._LHS
    @LHS.setter
    def LHS(self, val):
        if val is None:
            self._LHS = None
        else:
            if isinstance(val, (int, float)):
                val = val * units.dimensionless
            if isinstance(val, Q_):
                term1 = Term()
                term1.constant = val
                val = term1
            if isinstance(val, (Variable, Constant)):
                val = val.toTerm()
            if isinstance(val,Term):
                val = val.toExpression()
            if isinstance(val, Expression):
                self._LHS = val
            else:
                raise ValueError('Invalid type for LHS.  Must be an Expression')
        try:
            self.units
        except DimensionalityError as e:
            if self.LHS.units != self.RHS.units:
                print(self.LHS)
                print(self.RHS)
                print(self.LHS.units)
                print(self.RHS.units)
                print(self.LHS.units == self.RHS.units)
                raise ValueError('Units between LHS and RHS do not match.  Declare a new constraint.')
        except AttributeError as e:
            # No RHS was present
            pass
# =====================================================================================================================
# The right hand side of the constraint getter and setter
# =====================================================================================================================
    @property
    def RHS(self):
        return self._RHS
    @RHS.setter
    def RHS(self, val):
        if val is None:
            self._RHS = None
        else:
            if isinstance(val, (int, float)):
                val = val * units.dimensionless
            if isinstance(val, Q_):
                term1 = Term()
                term1.constant = val
                val = term1
            if isinstance(val, (Variable, Constant)):
                val = val.toTerm()
            if isinstance(val,Term):
                val = val.toExpression()
            if isinstance(val, Expression):
                self._RHS = val
            else:
                raise ValueError('Invalid type for RHS.  Must be an Expression')
        
        try:
            self.units
        except DimensionalityError as e:
            if self.LHS.units != self.RHS.units:
                print(self.LHS)
                print(self.RHS)
                print(self.LHS.units)
                print(self.RHS.units)
                print(self.LHS.units == self.RHS.units)
                raise ValueError('Units between LHS and RHS do not match.  Declare a new constraint.')
        except AttributeError as e:
            # No LHS was present
            pass
# =====================================================================================================================
# The operator getter and setter
# =====================================================================================================================
    @property
    def operator(self):
        return self._operator
    @operator.setter
    def operator(self, val):
        if isinstance(val, str):
            if val in ['==','>=','<=']:
                self._operator = val
            else:
                raise ValueError('Invalid type for operator.  Must be a string which is one of: ==, >=, <=')
        else:
            raise ValueError('Invalid type for operator.  Must be a string which is one of: ==, >=, <=')
# =====================================================================================================================
# Swaps left and right
# =====================================================================================================================
    def swap(self):
        rtmp = copy.deepcopy(self.RHS)
        ltmp = copy.deepcopy(self.LHS)
        otmp = self.operator
        self.LHS = rtmp
        self.RHS = ltmp
        if otmp == '>=':
            self.operator = '<='
        elif otmp == '<=':
            self.operator = '>='
        else:
            self.operator = '=='
# =====================================================================================================================
# Overload addition
# =====================================================================================================================
    def __add__(self, other):
        if isinstance(other, (int, float, Q_, Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            slf.RHS = slf.RHS + other
            slf.LHS = slf.LHS + other
            slf.simplify()
            return slf
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __radd__(self, other):
        if isinstance(other, (int, float, Q_, Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            slf.RHS = other + slf.RHS 
            slf.LHS = other + slf.LHS
            slf.simplify()
            return slf
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Overload subtraction
# =====================================================================================================================
    def __sub__(self, other):
        if isinstance(other, (int, float, Q_, Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            slf.RHS = slf.RHS - other
            slf.LHS = slf.LHS - other
            slf.simplify()
            return slf
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))
# ---------------------------------------------------------------------------------------------------------------------
    def __rsub__(self, other):
        if isinstance(other, (int, float, Q_, Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            slf.RHS = other - slf.RHS
            slf.LHS = other - slf.LHS
            slf.simplify()
            return slf
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Overload multiplication
# =====================================================================================================================
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            slf = copy.deepcopy(self)
            slf.RHS = slf.RHS * other
            slf.LHS = slf.LHS * other
            if other < 0.0:
                if slf.operator == '>=':
                    slf.operator = '<='
                elif slf.operator == '<=':
                    slf.operator = '>='
                else:
                    pass
        elif isinstance(other, Q_):
            slf = copy.deepcopy(self)
            rTmp = slf.RHS * other
            lTmp = slf.LHS * other
            slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
            if other.magnitude < 0.0:
                if slf.operator == '>=':
                    slf.operator = '<='
                elif slf.operator == '<=':
                    slf.operator = '>='
                else:
                    pass
        elif isinstance(other, _Unit):
            slf = copy.deepcopy(self)
            rTmp = slf.RHS * other
            lTmp = slf.LHS * other
            return Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
        elif isinstance(other, (Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            if self.operator == '==':
                rTmp = slf.RHS * other
                lTmp = slf.LHS * other
                slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
            else:
                raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))

        slf.simplify()
        return slf
# ---------------------------------------------------------------------------------------------------------------------
    def __rmul__(self, other):
    # _Unit when used with rmul automatically converts to Q_, so not included
        if isinstance(other, (int, float, Q_)):
            return self.__mul__(other)
        elif isinstance(other, (Variable, Constant, Term, Expression)):
            slf = copy.deepcopy(self)
            if self.operator == '==':
                rTmp = slf.RHS * other
                lTmp = slf.LHS * other
                slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
                slf.simplify()
                return slf
            else:
                raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
        else:
            raise TypeError('Operation between a Constraint and an object of type %s is not defined.'%(type(other).__name__))
# =====================================================================================================================
# Overload division
# =====================================================================================================================
    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            slf = copy.deepcopy(self)
            slf.RHS = slf.RHS / other
            slf.LHS = slf.LHS / other
            if other < 0.0:
                if slf.operator == '>=':
                    slf.operator = '<='
                elif slf.operator == '<=':
                    slf.operator = '>='
                else:
                    pass
        elif isinstance(other, Q_):
            slf = copy.deepcopy(self)
            rTmp = slf.RHS / other
            lTmp = slf.LHS / other
            slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
            if other.magnitude < 0.0:
                if slf.operator == '>=':
                    slf.operator = '<='
                elif slf.operator == '<=':
                    slf.operator = '>='
                else:
                    pass
        elif isinstance(other, _Unit):
            slf = copy.deepcopy(self)
            rTmp = slf.RHS / other
            lTmp = slf.LHS / other
            slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
        elif isinstance(other, (Variable, Constant)):
            slf = copy.deepcopy(self)
            if self.operator == '==':
                rTmp = slf.RHS / other
                lTmp = slf.LHS / other
                slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
            else:
                raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
        elif isinstance(other, Term):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            othr.simplify()
            if othr.variables == []:
                slf = slf / othr.constant
            elif self.operator == '==':
                rTmp = slf.RHS / other
                lTmp = slf.LHS / other
                slf = Constraint(LHS = lTmp, operator = self.operator, RHS = rTmp)
            else:
                raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
        elif isinstance(other, Expression):
            slf = copy.deepcopy(self)
            othr = copy.deepcopy(other)
            othr.simplify()
            if len(othr.terms) == 1:
                if othr.terms[0].variables == []:
                    slf = slf / othr.terms[0].constant
                else:
                    raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
            else:
                raise TypeError('Operation between an inequality constraint and an object of type %s is ambiguous and therefore not defined.'%(type(other).__name__))
        else:
            raise TypeError('Operation between a constraint and an object of type %s is not defined.'%(type(other).__name__))

        slf.simplify()
        return slf
# ---------------------------------------------------------------------------------------------------------------------
    def __rtruediv__(self, other):
        raise TypeError('Division by a constraint is not defined')
# =====================================================================================================================
# Overload negative sign
# =====================================================================================================================
    def __neg__(self):
        slf = copy.deepcopy(self)
        return -1.0 * slf
# =====================================================================================================================
# Determines if two constraints are equal
# =====================================================================================================================   
    def isEqual(self, other):
        slf = self.toZeroed()
        othr = other.toZeroed()
        if slf.operator == '==':
            if othr.operator == '>=' or othr.operator == '<=':
                return False
            elif othr.operator == '==':
                if slf.LHS.isEqual(othr.LHS) and slf.RHS.isEqual(othr.RHS):
                    return True
                elif slf.LHS.isEqual(othr.RHS) and slf.RHS.isEqual(othr.LHS):
                    return True
                else:
                    return False
            else:
                raise AttributeError('Invalid Operator Type')
        elif slf.operator == '>=':
            if othr.operator == '==':
                return False
            elif othr.operator == '>=':
                if slf.LHS.isEqual(othr.LHS) and slf.RHS.isEqual(othr.RHS):
                    return True
                else:
                    return False
            elif othr.operator == '<=':
                if slf.LHS.isEqual(othr.RHS) and slf.RHS.isEqual(othr.LHS):
                    return True
                else:
                    return False
            else:
                raise AttributeError('Invalid Operator Type')
        elif slf.operator == '<=':
            if othr.operator == '==':
                return False
            elif othr.operator == '>=':
                if slf.LHS.isEqual(othr.RHS) and slf.RHS.isEqual(othr.LHS):
                    return True
                else:
                    return False
            elif othr.operator == '<=':
                if slf.LHS.isEqual(othr.LHS) and slf.RHS.isEqual(othr.RHS):
                    return True
                else:
                    return False
            else:
                raise AttributeError('Invalid Operator Type')
        else:
            raise AttributeError('Invalid Operator Type')
# =====================================================================================================================
# Simplifies both sides independently
# =====================================================================================================================  
    def simplify(self):
        self.RHS.simplify()
        self.LHS.simplify()
# =====================================================================================================================    
# Substitutes in a dict of quantities
# =====================================================================================================================    
    def substitute(self, substitutions = {}, skipConstants=False):
        newLHS = self.LHS.substitute(substitutions = substitutions, skipConstants = skipConstants)
        newRHS = self.RHS.substitute(substitutions = substitutions, skipConstants = skipConstants)
        newCon = Constraint(LHS = newLHS, operator = self.operator, RHS = newRHS)
        newCon.simplify()
        ck1 = len(newCon.LHS.terms) == 1
        ck2 = len(newCon.RHS.terms) == 1
        if ck1 and ck2:
            ick1 = len(newCon.LHS.terms[0].variables) == 0
            ick2 = len(newCon.RHS.terms[0].variables) == 0
            if ick1 and ick2:
                if newCon.operator == '==':
                    if newCon.LHS.terms[0].constant != newCon.RHS.terms[0].constant:
                        raise ValueError('A substitution resulted in an invalid equality constraint')
                if newCon.operator == '>=':
                    if newCon.LHS.terms[0].constant < newCon.RHS.terms[0].constant:
                        raise ValueError('A substitution resulted in an invalid inequality constraint')
                if newCon.operator == '<=':
                    if newCon.LHS.terms[0].constant > newCon.RHS.terms[0].constant:
                        raise ValueError('A substitution resulted in an invalid inequality constraint')
        return newCon
# =====================================================================================================================
# Evaluates the variable to a quantity by taking in a design point
# =====================================================================================================================
    def evaluate(self, x_star):
        newLHS = self.LHS.evaluate(x_star)
        newRHS = self.RHS.evaluate(x_star)
        if self.operator == '==':
            rv = newLHS == newRHS
        if self.operator == '>=':
            rv = newLHS >= newRHS
        if self.operator == '<=':
            rv = newLHS <= newRHS
        return rv
# =====================================================================================================================    
# Sets a constraint so that the RHS is zero, returns copy
# ===================================================================================================================== 
    def toZeroed(self):
        con = copy.deepcopy(self)
        if con.operator == '<=':
            con.swap()
        for tm in con.RHS.terms:
            con -= tm
        return con
# =====================================================================================================================    
# Sets a constraint so that all terms have a positive leading constant, returns copy
# ===================================================================================================================== 
    def toAllPositive(self):
        con = self.toZeroed()
        for tm in con.LHS.terms:
            if tm.constant < 0 :
                con -= tm
        return con
# =====================================================================================================================    
# Sets a constraint so that the RHS is 1, with the LHS being a valid posynomial. Returns copy, does not enforce
# GP compatibility (posynomial equalities will pass through by default).  Setting the GPcompatible flag to true
# will return a list of GP compatible constraint(s).
# ===================================================================================================================== 
    def toPosynomial(self, GPCompatible = False, bypassChecks = False, leaveMonomialEqualities = False):
        con = self.toZeroed()
        csts = [con.LHS.terms[i].substitute(preserveType=True).constant for i in range(0,len(con.LHS.terms))]
        posTF = [csts[i] > 0 for i in range(0,len(csts))]
        numberPositive = sum(posTF)
        negTF = [csts[i] < 0 for i in range(0,len(csts))]
        numberNegative = sum(negTF)
        if con.operator == '==' and numberNegative == 1 and numberPositive != 1 :
            # Flips a posynomial equality if terms are on wrong side
            con *= -1
            csts = [con.LHS.terms[i].substitute(preserveType=True).constant for i in range(0,len(con.LHS.terms))]
            posTF = [csts[i] > 0 for i in range(0,len(csts))]
            numberPositive = sum(posTF)
            negTF = [csts[i] < 0 for i in range(0,len(csts))]
            numberNegative = sum(negTF)
        if numberPositive == 1:
            pix = posTF.index(True)
            newLHS = con.LHS / con.LHS.terms[pix]
            newRHS = con.RHS / con.LHS.terms[pix]
            con = Constraint(LHS = newLHS, operator = con.operator, RHS = newRHS)
            con = -1.0 * con
            con = con + 1.0 * con.units
            if con.operator == '==':
                con /= con.RHS
        else:
            if not bypassChecks:
                raise ValueError('This constraint cannot be written as a posynomial: %s'%(str(self)))
        if bypassChecks and numberNegative==1 and numberPositive!=1:
            nix = negTF.index(True)
            newLHS = con.LHS / con.LHS.terms[nix]
            newRHS = con.RHS / con.LHS.terms[nix]
            if con.operator == '>=':
                newOperator =  '<='
            if con.operator == '<=':
                newOperator =  '>='
            if con.operator == '==':
                newOperator =  '=='
            con = Constraint(LHS = newLHS, operator = newOperator, RHS = newRHS)
            con = -1.0 * con
            con = con + 1.0 * con.units
        # if bypassChecks and numberNegative>=1 and numberPositive>=1:
        #     con = con + 1.0 * con.units


        if GPCompatible:
            if con.operator == '==':
                if len(con.LHS.terms) > 1:
                    raise ValueError('This posynomial equality constraint is not GP compatible: %s'%(str(self)))
                else:
                    if leaveMonomialEqualities:
                        constraintSet = [con]
                    else:
                        con1 = copy.deepcopy(con)
                        con2 = copy.deepcopy(con)
                        con1.operator = '>='
                        con1.swap()
                        newLHS = con1.LHS / con1.RHS
                        newRHS = con1.RHS / con1.RHS
                        con1 = Constraint(LHS = newLHS, operator = con1.operator, RHS = newRHS)
                        con2.operator = '<='
                        constraintSet = [con1, con2]
            else:
                constraintSet = [con]
            return constraintSet
        else:
            return con
# =====================================================================================================================    
# Sets constraints to con == 1 or similar
# ===================================================================================================================== 
    def toLSQP_StandardForm(self):
        con = self.toAllPositive()
        if con.operator == '>=':
            con.swap()
            
        if len(con.RHS.terms) == 1:
            newLHS = con.LHS/con.RHS.terms[0]
            newRHS = con.RHS/con.RHS.terms[0]
            newCon = Constraint(LHS = newLHS, operator = con.operator, RHS = newRHS)
            return newCon

        for tm in con.RHS.terms:
            con -= tm
        con += 1*con.units
        con /= 1*con.units
        return con
# =====================================================================================================================    
# Sets a constraint so that the RHS is 0, with the LHS being a valid signomial.  Due to the way this has been written,
# this is identical to the toZero function, but may change over time.
# ===================================================================================================================== 
    def toSignomial(self):
        return self.toZeroed()
# =====================================================================================================================    
# Sets a constraint so that the RHS is 0, with the LHS being a quadratic function. 
# ===================================================================================================================== 
    def toQuadratic(self):
        con = self.toZeroed()
        subbedCon = con.substitute()

        VEL = [ [],
                [0.0],
                [1.0],
                [2.0],
                [1.0, 1.0],
                ]
        for tm in subbedCon.LHS.terms:
            if tm.exponents not in VEL:
                raise ValueError('This constraint cannot be written as quadratic: %s'%(str(self)))

        return con
# =====================================================================================================================    
# Sets a constraint so that the RHS is 0, with the LHS being a linear function. 
# ===================================================================================================================== 
    def toLinear(self, bypassLinearCheck = False):
        con = self.toZeroed()
        subbedCon = con.substitute()

        if not bypassLinearCheck:
            VEL = [ [],
                    [0.0],
                    [1.0],
                    ]
            for tm in subbedCon.LHS.terms:
                if tm.exponents not in VEL:
                    raise ValueError('This constraint cannot be written as linear: %s'%(str(self)))

        return con
# =====================================================================================================================    
# Approximates a signomial as GP compatible by replacing the RHS with an affine approximation
# ===================================================================================================================== 
    def approximateConvexifiedSignomial(self, point):
        def apprxMonomial(expr, point):
            expr = copy.deepcopy(expr)
            # Evaluate expression at point
            fval = expr.substitute(point)
            if isinstance(fval, Q_):
                apprx = fval * units.dimensionless
            else:
                raise ValueError('Insufficient variables provided in the reference point, variables remain in the function evaluation')

            # Get Variables on the expression
            exprv = []
            for term in expr.substitute(preserveType = True).terms:
                for vr in term.variables:
                    exprv.append(vr)

            for vr in exprv:
                # Compute partial derivative
                partial = Term()
                partial.constant = 0.0 * expr.units/vr.units
                for tm in expr.terms:
                    tfl1 = [tm.variables[i].isEqual(vr) for i in range(0,len(tm.variables))]
                    if any(tfl1):
                        newTerm = copy.deepcopy(tm)
                        tfl = [vr.isEqual(newTerm.variables[i]) for i in range(0,len(newTerm.variables))]
                        if True in tfl:
                            ix = tfl.index(True)
                            newTerm.constant *= float(newTerm.exponents[ix]) 
                            newTerm.exponents[ix] = newTerm.exponents[ix] - 1.0
                            # if partial.units == units.dimensionless:
                                # partial *= newTerm.units
                            partial += newTerm
                            
                partial.simplify()

                a = point[vr.name] / fval * partial.substitute(point)
                a = a.to('dimensionless').magnitude
                apprx *= (vr/point[vr.name])**a
            return apprx

        con = self.toZeroed()
        if con.operator == '>=':
            con *= -1
        for tm in con.LHS.terms:
            if tm.constant.to_base_units().magnitude < 0:
                con -= tm
                
        con.simplify()

        if con.operator == '==':
            # https://convex.mit.edu/publications/SignomialEquality.pdf
            # # Method A
            # con1 = copy.deepcopy(con)
            # con2 = copy.deepcopy(con)
            # con1.operator = '>='
            # con1.swap()
            # con1.RHS = apprxMonomial(con1.RHS, point)
            # con2.operator = '<='
            # con2.RHS = apprxMonomial(con2.RHS, point)
            # CS = [con1, con2]

            # # Method B
            # con1 = copy.deepcopy(con)
            # con2 = copy.deepcopy(con)
            # con1.operator = '>='
            # con1.swap()
            # con1.RHS = apprxMonomial(con1.RHS, point)
            # con2.operator = '<='
            # con2.LHS *= 0.95
            # con2.RHS = apprxMonomial(con2.RHS, point)
            # CS = [con1, con2]

            # Method C
            CS = []
            LA = apprxMonomial(con.LHS, point)
            RA = apprxMonomial(con.RHS, point)
            c2 = Constraint(LHS = LA, operator = con.operator, RHS = RA)
            CS.append(c2)

        else:
            CS = []
            apprx = apprxMonomial(con.RHS, point)
            # try:
            c2 = Constraint(LHS = con.LHS, operator = con.operator, RHS = apprx)
            # except:
            #     pass
                # print(con.LHS)
                # print(apprx)
                # print(con.LHS.units)
                # print(apprx.units)
                # print(con.LHS.units == apprx.units)
                # cody
            CS.append(c2)

        return CS
# =====================================================================================================================
# =====================================================================================================================

# toQuadraticSet--> Every signomial can be written as a set of quadratic constraints with substituted variables 
# (ie, x^3 -> y = x^2, x^3 = x*y), and thus we can write a tool that breaks signomials into their quadratic sets


# toSlack... (constrained=False) --> Adds a slack variable to each constraint.   May be different for GP constraints

# approximateLinear(point)
# approximateQuadratic(point)
# approximateTaylor(point, order) # special case is linear, quadratic
# approximateMonomial(point)


