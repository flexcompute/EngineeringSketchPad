# =====================================================================================================================
# Import Statements
# =====================================================================================================================
import copy
import warnings
import os
import shutil
from corsairlite import units
from corsairlite.core.dataTypes import TypeCheckedList, Container
from corsairlite.optimization.constant import Constant
from corsairlite.optimization import Variable
from corsairlite.optimization.term import Term
from corsairlite.optimization.expression import Expression
from corsairlite.optimization.runtimeConstraint import RuntimeConstraint
from corsairlite.optimization.constraint import Constraint
# =====================================================================================================================
# Define the Formulation Class
# =====================================================================================================================
class Formulation(object):
    def __init__(self, objective = None, constraints = [], solverOptions=None):
        self.skipEnforceBounds = False

        self.objective = objective
        self.constraints = constraints
        self.solverOptions = solverOptions
        
# =====================================================================================================================
# The printing function
# =====================================================================================================================
    def __repr__(self):
        pstr = ''
        pstr += 'Objective:\n'
        pstr += '\t'
        pstr += self.objective.__repr__()
        pstr += '\n'
        pstr += 'Constraints:\n'
        rcstr = 'Runtime Constraints:\n'
        for eqn in self.constraints:
            if isinstance(eqn, Constraint):
                pstr += '\t'
                pstr += eqn.__repr__()
                pstr += '\n'
            elif isinstance(eqn, RuntimeConstraint):
                for i in range(0,len(eqn.outputs)):
                    rcstr += '\t'
                    rcstr += eqn.outputs[i].name
                    rcstr += ' %s '%(eqn.operators[i])
                    rcstr += 'f('
                    for ipt in eqn.inputs:
                        rcstr += ipt.name
                        rcstr += ', '
                    rcstr = rcstr[0:-2]
                    rcstr += ')'
                    rcstr += '\n'
        return pstr + rcstr
# =====================================================================================================================
# The objective function getter and setter
# =====================================================================================================================
    @property
    def objective(self):
        return self._objective
    @objective.setter
    def objective(self, val):
        if val is None:
            self._objective = None
        else:
            if isinstance(val, Variable):
                try:
                    val = val.toTerm()
                except:
                    pass
            if isinstance(val,Term):
                tmp = Expression()
                tmp.terms.append(val)
                val = tmp
            if isinstance(val, Expression):
                self._objective = val
                self.objective.simplify()
                self.enforceBounds()
            else:
                raise ValueError('Invalid type for the objective.  Must be an Expression')
# =====================================================================================================================
# The constraint set getter and setter
# =====================================================================================================================
    @property
    def constraints(self):
        return self._constraints
    @constraints.setter
    def constraints(self, val):
        if isinstance(val, (list, tuple)):
            if all([isinstance(x, (Constraint,RuntimeConstraint)) for x in val]):
                self._constraints = val
                if not self.skipEnforceBounds:
                    self.enforceBounds()
            else:
                raise ValueError('One of the constraints is not a valid type.  Must be a Constraint or a RuntimeConstraint')
        else:
            raise ValueError('Constraint set must be a list or tuple.  TypeCheckedList and list types are preferred.')
# =====================================================================================================================
# The solverOptions
# =====================================================================================================================
    @property
    def solverOptions(self):
        return self._solverOptions
    @solverOptions.setter
    def solverOptions(self, val):
        if val is None:
            self._solverOptions = Container()
        elif isinstance(val, Container):
            self._solverOptions = val
        else:
            raise ValueError('The solverOptions attribute must be a Container')
# =====================================================================================================================
# Creates the bounds defined in the variable
# =====================================================================================================================
    def enforceBounds(self):
        if not hasattr(self, '_constraints'):
            self.constraints = []
        if self.objective is not None:
            for var in self.variables_only:
                bds = var.bounds
                if bds is not None:
                    con1 = var >= bds[0]
                    con2 = var <= bds[1]
                    c1a = True
                    c2a = True
                    for c in self.constraints:
                        if c.isEqual(con1):
                            c1a = False
                        if c.isEqual(con2):
                            c2a = False
                    if c1a:
                        self.constraints.append(con1)
                    if c2a:
                        self.constraints.append(con2)
# =====================================================================================================================
# Gets the equality constraints
# =====================================================================================================================
    @property
    def constraints_equality(self):
        cons = self.constraints
        equalityConstraints = TypeCheckedList(Constraint, [])
        for con in cons:
            if isinstance(con,Constraint):
                if con.operator == '==':
                    equalityConstraints.append(con)
        return equalityConstraints

# Alternate
    @property
    def equalityConstraints(self):
        return self.constraints_equality
# =====================================================================================================================
# Gets the inequality constraints
# =====================================================================================================================
    @property
    def constraints_inequality(self):
        cons = self.constraints
        inequalityConstraints = TypeCheckedList(Constraint, [])
        for con in cons:
            if isinstance(con,Constraint):
                if con.operator != '==':
                    inequalityConstraints.append(con)
        return inequalityConstraints

# Alternate
    @property
    def inequalityConstraints(self):
        return self.constraints_inequality
# =====================================================================================================================
# Function for sorting out bounding constraints
# =====================================================================================================================
    def _boundingConstraintFunction(self, mode):
        cons = copy.deepcopy(self.constraints)
        bounds = TypeCheckedList(Constraint, [])
        otherConstraints = TypeCheckedList((Constraint, RuntimeConstraint), [])
        for con in cons:
            if isinstance(con,Constraint):
                isBound = False
                cpy = con.toZeroed().substitute()
                cpyNS = con.toZeroed()
                if isinstance(cpy.LHS, Expression):
                    if len(cpy.LHS.terms) == 1:
                        if len(cpy.LHS.terms[0].hiddenConstants) == 0:
                            if len(cpy.LHS.terms[0].variables) == 1:
                                if isinstance(cpy.LHS.terms[0].variables[0],Variable):
                                    if cpy.LHS.terms[0].exponents == [1]:
                                        v = cpy.LHS.terms[0].variables[0]
                                        cst = cpy.LHS.terms[0].constant
                                        cpyNS /= cst
                                        cpyNS *= -1
                                        cpyNS += v
                                        cpyNS.swap()
                                        bounds.append(cpyNS)
                                        isBound = True


                    if len(cpy.LHS.terms) == 2:
                        if len(cpy.LHS.terms[0].hiddenConstants) == 0:
                            if len(cpy.LHS.terms[0].variables) == 1:
                                if len(cpy.LHS.terms[1].variables) == 0:
                                    if isinstance(cpy.LHS.terms[0].variables[0],Variable):
                                        if cpy.LHS.terms[0].exponents == [1]:
                                            v = cpy.LHS.terms[0].variables[0]
                                            cst = cpy.LHS.terms[0].constant
                                            cpyNS /= cst
                                            cpyNS *= -1
                                            cpyNS += v
                                            cpyNS.swap()
                                            bounds.append(cpyNS)
                                            isBound = True
                        if len(cpy.LHS.terms[1].hiddenConstants) == 0:
                            if len(cpy.LHS.terms[0].variables) == 0:
                                if len(cpy.LHS.terms[1].variables) == 1:
                                    if isinstance(cpy.LHS.terms[1].variables[0],Variable):
                                        if cpy.LHS.terms[1].exponents == [1]:
                                            v = cpy.LHS.terms[1].variables[0]
                                            cst = cpy.LHS.terms[1].constant
                                            cpyNS /= cst
                                            cpyNS *= -1
                                            cpyNS += v
                                            cpyNS.swap()
                                            bounds.append(cpyNS)
                                            isBound = True

                    if not isBound:
                        otherConstraints.append(con)
            elif isinstance(con,RuntimeConstraint):
                otherConstraints.append(con)
            else:
                raise ValueError('Invalid constraint type')
                

        if mode == 5:
            return otherConstraints

        names = [vbl.name for vbl in self.variables_only]
        lbs = [[] for i in range(0,len(names))]
        ubs = [[] for i in range(0,len(names))]
        ebs = [[] for i in range(0,len(names))]

        lbs_us = [[] for i in range(0,len(names))]
        ubs_us = [[] for i in range(0,len(names))]
        ebs_us = [[] for i in range(0,len(names))]
        
        for constraint in bounds:
            vl = constraint.RHS.substitute()
            vl_us = constraint.RHS
            vr = constraint.LHS.terms[0].variables[0]
            ix = names.index(vr.name)
            if constraint.operator == '==':
                ebs[ix].append(vl)
                ebs_us[ix].append(vl_us)
            if constraint.operator == '>=':
                lbs[ix].append(vl)
                lbs_us[ix].append(vl_us)
            if constraint.operator == '<=':
                ubs[ix].append(vl)
                ubs_us[ix].append(vl_us)


        if mode == 1:
            finalBounds = TypeCheckedList(Constraint, [])
            for i in range(0,len(self.variables_only)):
                v = self.variables_only[i]
                lbvs = lbs[i]
                ubvs = ubs[i]
                ebvs = ebs[i]
                if len(lbvs) > 0:
                    lb = max(lbvs)
                    finalBounds.append(Constraint(v,'>=',lb))
                if len(ubvs) > 0:
                    ub = min(ubvs)
                    finalBounds.append(Constraint(v,'<=',ub))
                if len(ebvs) > 0:
                    if any([ebvs[ii] != ebvs[0] for ii in range(0,len(ebvs))]):
                        raise ValueError('Incompatible equality bounding constraints for variable %s'%(v.name))
                    else:
                        finalBounds.append(Constraint(v,'==',ebvs[0]))
            return finalBounds

        if mode == 2:
            outputMessages = []
            finalBounds = TypeCheckedList(Constraint, [])
            for i in range(0,len(self.variables_only)):
                v = self.variables_only[i]
                lbvs = lbs[i]
                ubvs = ubs[i]
                ebvs = ebs[i]
                lbvs_us = lbs_us[i]
                ubvs_us = ubs_us[i]
                ebvs_us = ebs_us[i]

                if len(lbvs) > 0:
                    lb = max(lbvs)
                    tfl = [lbvs[ii] == lb for ii in range(0,len(lbvs))]
                    if sum(tfl) > 1:
                        raise ValueError('Two constraints are imposing the same lower bound on variable %s'%(v.name))
                    ix = tfl.index(True)
                    finalBounds.append(Constraint(v,'>=',lbvs_us[ix]))
                if len(ubvs) > 0:
                    ub = min(ubvs)
                    tfl = [ubvs[ii] == ub for ii in range(0,len(ubvs))]
                    if sum(tfl) > 1:
                        raise ValueError('Two constraints are imposing the same upper bound on variable %s'%(v.name))
                    ix = tfl.index(True)
                    finalBounds.append(Constraint(v,'<=',ubvs_us[ix]))
                if len(ebvs) > 0:
                    if any([ebvs[ii] != ebvs[0] for ii in range(0,len(ebvs))]):
                        raise ValueError('Incompatible equality bounding constraints for variable %s'%(v.name))
                    if len(ebvs) > 1:
                        raise ValueError('Two constraints are imposing the same upper bound on variable %s'%(v.name))
                    finalBounds.append(Constraint(v,'==',ebvs_us[0]))

            return finalBounds

        if mode == 4:
            lowerBounds = []
            upperBounds = []
            for i in range(0,len(self.variables_only)):
                v = self.variables_only[i]
                lbvs = lbs[i]
                ubvs = ubs[i]
                ebvs = ebs[i]

                if len(lbvs) > 0:
                    lb = max(lbvs)
                else:
                    lb = -float('inf') * v.units

                if len(ubvs) > 0:
                    ub = min(ubvs)
                else:
                    ub = float('inf') * v.units

                if len(ebvs) > 0:
                    if any([ebvs[ii] != ebvs[0] for ii in range(0,len(ebvs))]):
                        raise ValueError('Incompatible equality bounding constraints for variable %s'%(v.name))
                    eb = ebvs[0]
                    if eb < lb or eb > ub:
                        raise ValueError('Equality constraint does not fall between bounding constraints for variable %s'%(v.name))

                    lowerBounds.append(eb)
                    upperBounds.append(eb)
                else:
                    lowerBounds.append(lb)
                    upperBounds.append(ub)

            return lowerBounds, upperBounds


        if mode == 3:
            outputMessages = []
            finalBoundsMatrix = []
            for i in range(0,len(self.variables_only)):
                v = self.variables_only[i]
                lbvs = lbs[i]
                ubvs = ubs[i]
                ebvs = ebs[i]
                lbvs_us = lbs_us[i]
                ubvs_us = ubs_us[i]
                ebvs_us = ebs_us[i]

                if len(lbvs) > 0:
                    lb = max(lbvs)
                    tfl = [lbvs[ii] == lb for ii in range(0,len(lbvs))]
                    if sum(tfl) > 1:
                        outputMessages.append('Two constraints are imposing the same lower bound on variable %s.  Sensitivities may be innacurate.'%(v.name))
                    ix = tfl.index(True)
                    fbr = [Constraint(v,'>=',lbvs_us[ix]), v , 'lower']
                    finalBoundsMatrix.append(fbr)
                if len(ubvs) > 0:
                    ub = min(ubvs)
                    tfl = [ubvs[ii] == ub for ii in range(0,len(ubvs))]
                    if sum(tfl) > 1:
                        outputMessages.append('Two constraints are imposing the same upper bound on variable %s.  Sensitivities may be innacurate.'%(v.name))
                    ix = tfl.index(True)
                    fbr = [Constraint(v,'<=',ubvs_us[ix]), v , 'upper']
                    finalBoundsMatrix.append(fbr)
                if len(ebvs) > 0:
                    if any([ebvs[ii] != ebvs[0] for ii in range(0,len(ebvs))]):
                        raise ValueError('Incompatible equality bounding constraints for variable %s'%(v.name))
                    if len(ebvs) > 1:
                        outputMessages.append('Two constraints are imposing the same upper bound on variable %s.  Sensitivities may be innacurate.'%(v.name))
                    fbr = [Constraint(v,'==',ebvs_us[0]), v , 'equality']
                    finalBoundsMatrix.append(fbr)

            return finalBoundsMatrix, outputMessages
# =====================================================================================================================
# Gets the variable bounding constraints
# =====================================================================================================================
    @property
    def constraints_bounding(self):
        return self._boundingConstraintFunction(mode=1)

    # Alternate
    @property
    def boundingConstraints(self):
        return self.constraints_bounding
# =====================================================================================================================
# Gets the constraints which are not variable bounding constraints
# =====================================================================================================================
    @property
    def constraints_notBounding(self):
        return self._boundingConstraintFunction(mode=5)

    # Alternate
    @property
    def notBoundingConstraints(self):
        return self.constraints_notBounding
# =====================================================================================================================
# Returns a list of all the terms in the formulation in order of appearance
# =====================================================================================================================
    @property
    def terms(self):
        termList = TypeCheckedList(Term, [])

        # Want to include terms in the objective
        cons = copy.deepcopy(self.constraints)
        obj = copy.deepcopy(self.objective)
        cons = [obj == 0*obj.units] + cons

        for constraint in cons:
            if isinstance(constraint, Constraint):
                LHS = constraint.LHS
                RHS = constraint.RHS
                expressions = [LHS, RHS]
                for expression in expressions:
                    if isinstance(expression, Expression):
                        for term in expression.terms:
                            if term.variables != []:
                                termList.append(term)
                    else:
                        raise ValueError('Invalid expression type')
            elif isinstance(constraint, RuntimeConstraint):
                pass
            else:
                raise ValueError('Invalid constraint type')
        return termList
# =====================================================================================================================
# Returns a list of all the variables in the formulation in alphabetical order
# =====================================================================================================================
    @property
    def variables(self):
        variableList = TypeCheckedList((Variable,Constant), [])

        # Want to include variables in the objective
        # cons = copy.deepcopy(self.constraints)
        cons = self.constraints
        obj = copy.deepcopy(self.objective)
        cons = [obj == 0*obj.units] + cons

        for constraint in cons:
            if isinstance(constraint, Constraint):
                LHS = constraint.LHS
                RHS = constraint.RHS
                expressions = [LHS, RHS]
                for expression in expressions:
                    if isinstance(expression, Expression):
                        for term in expression.terms:
                            for var in term.variables_all:
                                variableList.append(var)
                    else:
                        raise ValueError('Invalid expression type')
            elif isinstance(constraint, RuntimeConstraint):
                variableList += constraint.inputs
                variableList += constraint.outputs
            else:
                raise ValueError ('Invalid constraint type')
                        
        names = [vbl.name for vbl in variableList]
        for i in range(0,len(names)):
            for j in range(0,len(names)):
                if i != j:
                    if names[i] == names[j] and not variableList[i].isEqual(variableList[j]):
                        raise ValueError('The name %s is duplicated.  Use unique variable and constant names.'%(names[i]))

        nms = [variableList[i].name for i in range(0,len(variableList))]
        dct = dict(zip(nms, variableList))
        kys = sorted(dct.keys())
        vls = [dct[kys[i]] for i in range(0,len(kys))]
        variableList = vls
        return variableList
# =====================================================================================================================
# Returns a list of all the variables in the formulation in alphabetical order, without constants
# =====================================================================================================================
    @property
    def variables_only(self):
        vrs =  self.variables
        csts = self.constants
        cl = []
        for i in range(0,len(vrs)):
            v = vrs[i]
            cklist = []
            for c in csts:
                if c.isEqual(v):
                    cklist.append(True)
                else:
                    cklist.append(False)
            if any(cklist):
                cl.append(False)
            else:
                cl.append(True)
        rv = [vrs[i] for i in range(0,len(vrs)) if cl[i]]
        return rv
# =====================================================================================================================
# Returns a list of all the constants in the formulation in alphabetical order
# =====================================================================================================================
    @property
    def constants(self):
        vrs = self.variables
        csts = TypeCheckedList(Constant, [])
        for v in vrs:
            if isinstance(v, Constant):
                csts.append(v)
        return csts
# =====================================================================================================================
# Substitutes all the Constants with their values
# =====================================================================================================================
    def substitute(self, substitutions = {}, skipConstants=False):
        # Set up a substitutions dict which contains variable names and values
        subs = copy.deepcopy(substitutions)
        kys = [vbl.name for vbl in self.variables]
        vls = [None for i in range(0,len(kys))]
        allSubsDict = dict(zip(kys, vls))
        varDict = dict(zip(kys, copy.deepcopy(self.variables)))
        for key, val in allSubsDict.items():
            vr = varDict[key]
            if isinstance(vr, Constant) and not skipConstants:
                allSubsDict[key] = vr.value
            if key in subs.keys():
                sv = subs[key]
                del subs[key]
                if isinstance(sv, (int, float)):
                    sv *= units.dimensionless
                if vr.units == sv.units:
                    allSubsDict[key] = sv
                else:
                    raise ValueError('Invalid units provided on substitution to %s.  Expected %s'%(key, str(vr.units)))
        if len(subs.keys()) > 0:
            warnings.warn('The following substitution keys were unused: %s'%(str(list(subs.keys()))))

        subsDict = {}
        for key, val in allSubsDict.items():
            if val is not None:
                subsDict[key] = val

        # Substitute using the subsDict
        newFormulation = Formulation()
        newObjective = self.objective.substitute(substitutions = subsDict, skipConstants=skipConstants)
        newFormulation.objective = newObjective

        newConstraints = []
        for con in self.constraints:
            if isinstance(con, Constraint):
                newCon = con.substitute(substitutions = subsDict, skipConstants=skipConstants)
                newConstraints.append(newCon)
            elif isinstance(con, RuntimeConstraint):
                newCon = con.substitute(substitutions = subsDict)
                newConstraints.append(newCon)
            else:
                raise ValueError('Invalid constraint type')
                
        for con in newConstraints:
            if isinstance(con,Constraint):
                con.to_base_units()

        newFormulation.constraints = newConstraints
        newFormulation.solverOptions = self.solverOptions

        return newFormulation
# =====================================================================================================================
# Sets all equations to the form: f(x) == 0, g(x) >= 0 
# =====================================================================================================================
    def toZeroed(self):
        newFormulation = Formulation()
        newConstraints = []
        for con in self.constraints:
            if isinstance(con,Constraint):
                con = con.toZeroed()
                con.to_base_units()
                newConstraints.append(con)
            elif isinstance(con,RuntimeConstraint):
                newConstraints.append(con)
            else:
                raise ValueError('Invalid constraint type')
            
        newFormulation.constraints = newConstraints
        newFormulation.objective = copy.deepcopy(self.objective)
        newFormulation.solverOptions = copy.deepcopy(self.solverOptions)
        return newFormulation
# =====================================================================================================================
# Sets all equations to the form f(x) >= 0 if possible, returns a valid linear program
# =====================================================================================================================
    def toLinearProgram(self, substituted = True):
        cons = self.constraints
        if any([not isinstance(con, Constraint) for con in cons]):
            raise ValueError('Cannot construct a linear program with a black boxed RuntimeConstraint')
        
        zeroed = self.toZeroed()
        if substituted:
            formulation = zeroed.substitute()
        else:
            formulation = zeroed
        newFormulation = Formulation()
        newConstraints = []
        for con in formulation.constraints:
            con = con.toLinear(bypassLinearCheck = not substituted)
            con *= -1
            con.to_base_units()
            for tm in con.LHS.terms:
                if tm.exponents == []:
                    con -= tm
                    con.simplify()
                    break
            newConstraints.append(con)

        newFormulation.constraints = newConstraints
        # Need to ensure objective is linear
        obj = copy.deepcopy(formulation.objective)
        if substituted:
            VEL = [ [],
                    [0.0],
                    [1.0],
                    ]
            for tm in obj.terms:
                if tm.exponents not in VEL:
                    raise ValueError('Objective is not a valid linear expression: %s'%(str(obj)))
        newFormulation.objective = obj
        newFormulation.solverOptions = copy.deepcopy(self.solverOptions)
        return newFormulation
# =====================================================================================================================
# Sets all equations to the form f(x) >= 1 if possible (ie, makes the formulation GP compatible)
# =====================================================================================================================
    def toGeometricProgram(self, substituted=True):
        cons = self.constraints
        if any([not isinstance(con, Constraint) for con in cons]):
            raise ValueError('Cannot construct a geometric program with a black boxed RuntimeConstraint')
        
        zeroed = self.toZeroed()
        if substituted:
            formulation = zeroed.substitute()
        else:
            formulation = zeroed
        newFormulation = Formulation()
        newConstraints = []
        for con in formulation.constraints:
            conSet = con.toPosynomial(GPCompatible = True, bypassChecks = not substituted)
            for sc in conSet:
                sc.to_base_units()
                newConstraints.append(sc)

        newFormulation.constraints = newConstraints

        # Need to ensure objective is posynomial
        obj = copy.deepcopy(formulation.objective)
        csts = [obj.terms[i].constant for i in range(0,len(obj.terms))]
        negTF = [csts[i] < 0 for i in range(0,len(csts))]
        numberNegative = sum(negTF)
        if substituted:
            if numberNegative == 0:
                newFormulation.objective = obj
            else:
                raise ValueError('Ojective is not a valid posynomial: %s'%(str(obj)))
        else:
            newFormulation.objective = obj
        newFormulation.solverOptions = copy.deepcopy(self.solverOptions)
        return newFormulation
# =====================================================================================================================
# Sets all equations to the form f(x) >= 0 if possible, returns a valid quadratic program
# =====================================================================================================================
    def toQuadraticProgram(self, substituted=True):
        cons = self.constraints
        if any([not isinstance(con, Constraint) for con in cons]):
            raise ValueError('Cannot construct a quadratic program with a black boxed RuntimeConstraint')
        
        zeroed = self.toZeroed()
        if substituted:
            formulation = zeroed.substitute()
        else:
            formulation = zeroed
        newFormulation = Formulation()
        newConstraints = []
        for con in formulation.constraints:
            con = con.toLinear(bypassLinearCheck = not substituted)
            con *= -1
            con.to_base_units()
            for tm in con.LHS.terms:
                if tm.exponents == []:
                    con -= tm
                    con.simplify()
                    break
            newConstraints.append(con)

        newFormulation.constraints = newConstraints
        # Need to ensure objective is quadratic
        obj = copy.deepcopy(formulation.objective)
        if substituted:
            VEL = [ [],
                    [0.0],
                    [1.0],
                    [2.0],
                    [1.0, 1.0],
                    ]
            for tm in obj.terms:
                if tm.exponents not in VEL:
                    raise ValueError('Objective is not a valid quadratic expression: %s'%(str(obj)))
        newFormulation.objective = obj
        newFormulation.solverOptions = copy.deepcopy(self.solverOptions)
        return newFormulation
# =====================================================================================================================
# Sets all equations to the form f(x) >= 0 if possible, returns a signomial progarm.  At the moment, this does nothing
# =====================================================================================================================
    def toSignomialProgram(self, substituted=True):
        cons = self.constraints
        if any([not isinstance(con, Constraint) for con in cons]):
            raise ValueError('Cannot construct a signomial program with a black boxed RuntimeConstraint')
        
        zeroed = self.toZeroed()
        if substituted:
            formulation = zeroed.substitute()
        else:
            formulation = zeroed
        newFormulation = Formulation()
        newConstraints = []
        for con in formulation.constraints:
            try:
                conSet = con.toPosynomial(GPCompatible = True, bypassChecks = not substituted)
                for sc in conSet:
                    sc.to_base_units()
                    newConstraints.append(sc)
            except:
                newConstraints.append(con)

        newFormulation.constraints = newConstraints

        # Need to ensure objective is posynomial
        obj = copy.deepcopy(self.objective)
        csts = [obj.terms[i].constant for i in range(0,len(obj.terms))]
        negTF = [csts[i] < 0 for i in range(0,len(csts))]
        numberNegative = sum(negTF)
        if numberNegative > 0:
            raise ValueError('Ojective is not a valid posynomial: %s'%(str(obj)))

        newFormulation.objective = copy.deepcopy(formulation.objective)
        newFormulation.solverOptions = copy.deepcopy(self.solverOptions)
        return newFormulation
# =====================================================================================================================
# Sets formulation to standard form for LSQP
# =====================================================================================================================
    def toLSQP_StandardForm(self, substituted=True):
        zeroed = self.toZeroed()
        # if substituted:
        formulation = zeroed.substitute()
        # else:
            # formulation = zeroed

        for i in range(0,len(formulation.constraints)):
            if isinstance(formulation.constraints[i],Constraint):
                formulation.constraints[i] = formulation.constraints[i].toLSQP_StandardForm()
        return formulation
# =====================================================================================================================
# Returns a string of the formulation matrix
# =====================================================================================================================
    def table(self, Ndecimal = 1):
        fstr = '%.' + str(Ndecimal) + 'f'
        
        cons = self.constraints
        if any([not isinstance(con, Constraint) for con in cons]):
            raise ValueError('Table cannot be constructed for a problem with runtime constraints')

        zeroed = self.toZeroed()
        subbed = zeroed.substitute()
        vbls = subbed.variables
        vblNames = [vb.name for vb in vbls]

        # constraint | term | constant | var1exp | var2exp | var3exp |
        terms = []
        tpl = []
        for con in subbed.constraints:
            # if not con.RHS.isEqual(0 * con.units):
            #     print(zeroed)
            #     print('\n\n\n\n')
            #     print(subbed)
            #     print('\n\n\n\n')
            #     print(con)
            if not con.RHS.isEqual(0 * con.units):
                raise ValueError('Constraint is not zeroed')
            terms.append(con.LHS.terms)
            if con.operator == '==':
                tpl.append('e')
            else:
                tpl.append('i')

        constraintList = []
        typeList = []
        termList = []
        constantList = []
        unitList = []
        expMat = []
        for i in range(0,len(terms)):
            for j in range(0,len(terms[i])):
                constraintList.append(i)
                typeList.append(tpl[i])
                termList.append(j)
                constantList.append(terms[i][j].constant)
                unitList.append(terms[i][j].units)
                rw = [0 for ii in range(0,len(vblNames))]
                for k in range(0,len(terms[i][j].variables)):
                    vr = terms[i][j].variables[k]
                    exp = terms[i][j].exponents[k]
                    varIndex = vblNames.index(vr.name)
                    rw[varIndex] = exp
                expMat.append(rw)

        data = {}
        data['constraintNumber'] = constraintList
        data['constraintType'] = typeList
        data['termNumber'] = termList 
        data['constant'] = constantList
        data['units'] = unitList
        data['variableNames'] = vblNames
        data['exponents'] = expMat

        cols = [[] for i in range(0,len(data['variableNames'])+3)]
        topRow = ['Index', 'Type', 'Term', 'Constant', 'Units'] + data['variableNames']
        cols = [[] for i in range(0,len(topRow))]
        for i in range(0,len(topRow)):
            cols[i].append(topRow[i])
        
        subbed.objective.to_base_units()
        for i in range(0,len(subbed.objective.terms)):
            rw = []
            tm = subbed.objective.terms[i]
            rw.append(0)
            rw.append('o')
            rw.append(i)
            rw.append(tm.constant.magnitude)
            rw.append(tm.units)
            for j in range(0,len(vblNames)):
                rw.append(0)
            for j in range(0,len(tm.variables)):
                vr = tm.variables[j]
                exp = tm.exponents[j]
                varIndex = vblNames.index(vr.name)
                rw[varIndex+5] = exp
            for k in range(0,len(rw)):
                cols[k].append(rw[k])

        for i in range(0,len(data['constraintNumber'])):
            rw = []
            rw.append(data['constraintNumber'][i]+1)
            rw.append(data['constraintType'][i])
            rw.append(data['termNumber'][i])
            rw.append(data['constant'][i].magnitude)
            rw.append(data['units'][i])
            rw += data['exponents'][i]
            for k in range(0,len(rw)):
                cols[k].append(rw[k])

        # Column 0 is all integers, no action needed other than convert to string
        for i in range(1,len(cols[0])):
            cols[0][i] = str(cols[0][i])
        # Column 1 is all strings already
        # Column 2 is all integers, no action needed other than convert to string
        for i in range(1,len(cols[2])):
            cols[2][i] = str(cols[2][i])
        # Column 3 is all floats in general, but can be all ints
        if any([not float(cols[3][i]).is_integer() for i in range(1, len(cols[3]))]):
            for i in range(1,len(cols[3])):
                    if float(cols[3][i]).is_integer():
                        cols[3][i] = '%d'%(cols[3][i])+ ' '*(Ndecimal+1)
                    else:
                        cols[3][i] = fstr%(cols[3][i])
        else:
            for i in range(1,len(cols[3])):
                cols[3][i] = '%d'%(cols[3][i])
        # Column 4 is all units and therefore needs a string conversion
        for i in range(1,len(cols[4])):
            cols[4][i] = str(cols[4][i])
            if cols[4][i] == 'dimensionless':
                cols[4][i] = '-'
        # Remaining columns are similar to column 3
        for ii in range(5,len(cols)):
            if any([not float(cols[ii][i]).is_integer() for i in range(1, len(cols[ii]))]):
                for i in range(1,len(cols[ii])):
                        if float(cols[ii][i]).is_integer():
                            cols[ii][i] = '%d'%(cols[ii][i])+ ' '*(Ndecimal+1)
                        else:
                            cols[ii][i] = fstr%(cols[ii][i])
            else:
                for i in range(1,len(cols[ii])):
                    cols[ii][i] = '%d'%(cols[ii][i])


        columnWidths = []
        for col in cols:
            mw = max([len(itm) for itm in col])
            columnWidths.append(mw)
        # -----
        # Sets variable columns to be the same width 
        varWidth = max(columnWidths[5:])
        w = [varWidth for i in range(0,len(columnWidths))]
        for i in [0,1,2,3,4]:
            w[i] = columnWidths[i]
        columnWidths = w
        # -----

        line = ''
        for i in range(0,len(cols)):
            line += cols[i][0].center(columnWidths[i]+2)
            line += '|'
        tableWidth = len(line)
        line += '\n'
        
        el = '='*tableWidth + '\n'
        sl = '-'*tableWidth + '\n'

        pstr = el
        pstr += line

        for i in range(1,len(cols[0])):
            line = ''
            for j in range(0,len(cols)):
                if j in [1,4]:
                    line += cols[j][i].center(columnWidths[j]+2)
                else:
                    s = cols[j][i].rjust(columnWidths[j])
                    line += s.center(len(s)+2)

                line += '|'
            line += '\n'
            pstr += line

        pstr += el
        return pstr
# =====================================================================================================================
# Generates the solver options container
# =====================================================================================================================
    def generateSolverOptions(self):
        if hasattr(self.solverOptions, 'solveType'):
            if not hasattr(self.solverOptions,'processID'):
                self.solverOptions.processID = 0
            if self.solverOptions.solveType.lower() in ['lp','qp','gp','sp','pccp','sqp','lsqp','slcp']:
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                if not hasattr(self.solverOptions,'returnRaw'):
                    self.solverOptions.returnRaw = False
                if not hasattr(self.solverOptions,'unknownTolerance'):
                    self.solverOptions.unknownTolerance = 1e-4
            if self.solverOptions.solveType.lower() in ['sp', 'pccp']:
                if not hasattr(self.solverOptions,'x0'):
                    self.solverOptions.x0 = None
                if not hasattr(self.solverOptions,'maxIter'):
                    self.solverOptions.maxIter = 100
                if not hasattr(self.solverOptions,'trustRegions'):
                    self.solverOptions.trustRegions = None
                if not hasattr(self.solverOptions,'tol'):
                    self.solverOptions.tol = 1e-3
                if not hasattr(self.solverOptions,'convergenceCheckIterations'):
                    self.solverOptions.convergenceCheckIterations = 50
                if not hasattr(self.solverOptions,'convergenceCheckSlope'):
                    self.solverOptions.convergenceCheckSlope = -0.01
            if self.solverOptions.solveType.lower() == 'pccp':
                if not hasattr(self.solverOptions,'penaltyExponent'):
                    self.solverOptions.penaltyExponent = 5
                if not hasattr(self.solverOptions,'slackTol'):
                    self.solverOptions.slackTol = 1e-3
            if self.solverOptions.solveType.lower() == 'ga':
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                if not hasattr(self.solverOptions,'tol'):
                    self.solverOptions.tol = 1e-6
                if not hasattr(self.solverOptions,'constraintTol'): 
                    self.solverOptions.constraintTol = 1e-4
                if not hasattr(self.solverOptions,'populationSize'):
                    self.solverOptions.populationSize = 20*len(self.variables_only)
                if not hasattr(self.solverOptions,'processType'):
                    self.solverOptions.processType = 'parallel'
                if not hasattr(self.solverOptions,'penaltyConstant'):
                    self.solverOptions.penaltyConstant = 1e6
                if not hasattr(self.solverOptions,'nbe'):
                    self.solverOptions.nbe = 6
                if not hasattr(self.solverOptions,'nbm'):
                    self.solverOptions.nbm = 8
                if not hasattr(self.solverOptions,'probabilityOfMutation'):
                    self.solverOptions.probabilityOfMutation = 0.3
                if not hasattr(self.solverOptions,'Ncrossovers'):
                    self.solverOptions.Ncrossovers = 3
                if not hasattr(self.solverOptions,'maxGenerations'):
                    self.solverOptions.maxGenerations = 100
                if not hasattr(self.solverOptions,'convergedGenerations'):
                    self.solverOptions.convergedGenerations = 20
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'GA_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'GA_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() == 'pso':
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                if not hasattr(self.solverOptions,'tol'):
                    self.solverOptions.tol = 1e-6
                if not hasattr(self.solverOptions,'constraintTol'): 
                    self.solverOptions.constraintTol = 1e-4
                if not hasattr(self.solverOptions,'swarmSize'):
                    self.solverOptions.swarmSize = 20
                if not hasattr(self.solverOptions,'processType'):
                    self.solverOptions.processType = 'parallel'
                if not hasattr(self.solverOptions,'penaltyConstant'):
                    self.solverOptions.penaltyConstant = 1e6
                if not hasattr(self.solverOptions,'dt'):
                    self.solverOptions.dt = 0.2
                if not hasattr(self.solverOptions,'inertiaFactor'):
                    self.solverOptions.inertiaFactor = 0.9
                if not hasattr(self.solverOptions,'selfConfidence'):
                    self.solverOptions.selfConfidence = 1.75
                if not hasattr(self.solverOptions,'swarmConfidence'):
                    self.solverOptions.swarmConfidence = 2.25
                if not hasattr(self.solverOptions,'maxTimeSteps'):
                    self.solverOptions.maxTimeSteps = 100
                if not hasattr(self.solverOptions,'convergedTimeSteps'):
                    self.solverOptions.convergedTimeSteps = 5
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'PSO_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'PSO_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() == 'sa':
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                # if not hasattr(self.solverOptions,'tol'): # Not needed for SA which terminates after fixed number of steps
                #     self.solverOptions.tol = 1e-6
                if not hasattr(self.solverOptions,'constraintTol'): 
                    self.solverOptions.constraintTol = 1e-4
                if not hasattr(self.solverOptions,'maximumExponent'):
                    self.solverOptions.maximumExponent = 100
                if not hasattr(self.solverOptions,'Nsteps'):
                    self.solverOptions.Nsteps = 1000
                if not hasattr(self.solverOptions,'scheduleBase'):
                    self.solverOptions.scheduleBase = 0.85
                if not hasattr(self.solverOptions,'initialTemperature'):
                    self.solverOptions.initialTemperature = 100
                if not hasattr(self.solverOptions,'penaltyConstant'):
                    self.solverOptions.penaltyConstant = 1e6
            if self.solverOptions.solveType.lower() == 'nms':
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                if not hasattr(self.solverOptions,'tol'): 
                    self.solverOptions.tol = 1e-4
                if not hasattr(self.solverOptions,'constraintTol'): 
                    self.solverOptions.constraintTol = 1e-4
                if not hasattr(self.solverOptions,'maxIterations'): 
                    self.solverOptions.maxIterations = 1000
                if not hasattr(self.solverOptions,'processType'):
                    self.solverOptions.processType = 'parallel'
                if not hasattr(self.solverOptions,'h'):
                    self.solverOptions.h = 0.05
                if not hasattr(self.solverOptions,'alpha'):
                    self.solverOptions.alpha = 1.0
                if not hasattr(self.solverOptions,'beta'):
                    self.solverOptions.beta = 0.5
                if not hasattr(self.solverOptions,'gamma'):
                    self.solverOptions.gamma = 2.0
                if not hasattr(self.solverOptions,'delta'):
                    self.solverOptions.delta = 0.5
                if not hasattr(self.solverOptions,'greedy'):
                    self.solverOptions.greedy = 'expansion' # 'expansion' or 'minimization'
                if not hasattr(self.solverOptions,'penaltyConstant'):
                    self.solverOptions.penaltyConstant = 1e6
                if not hasattr(self.solverOptions,'convergenceCheckIterations'):
                    self.solverOptions.convergenceCheckIterations = 10
                if not hasattr(self.solverOptions,'convergenceCheckSlope'):
                    self.solverOptions.convergenceCheckSlope = -0.01
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'NMS_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'NMS_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() in ['sqp','lsqp','slcp']:
                if not hasattr(self.solverOptions,'x0'):
                    self.solverOptions.x0 = None
                if not hasattr(self.solverOptions,'solver'):
                    self.solverOptions.solver = None
                if not hasattr(self.solverOptions,'maxIterations'): 
                    self.solverOptions.maxIterations = 500
                if not hasattr(self.solverOptions,'penaltyConstant'):
                    self.solverOptions.penaltyConstant = 1e15
                if not hasattr(self.solverOptions,'constraintTolerance'):
                    self.solverOptions.constraintTolerance = 1e-3
                if not hasattr(self.solverOptions,'stepMagnitudeTolerance'):
                    self.solverOptions.stepMagnitudeTolerance = 1e-4
                if not hasattr(self.solverOptions,'lagrangianGradientTolerance'):
                    self.solverOptions.lagrangianGradientTolerance = 1e-6
                if not hasattr(self.solverOptions,'eta'):
                    self.solverOptions.eta = 1e-4 #(0,0.5)
                if not hasattr(self.solverOptions,'zeta'):
                    self.solverOptions.zeta = 0.9 #(eta,1)
                if not hasattr(self.solverOptions,'tau'):
                    self.solverOptions.tau = 0.9 #(0,1.0)
                if not hasattr(self.solverOptions,'rho'):
                    self.solverOptions.rho = 0.8 #(0,1)
                if not hasattr(self.solverOptions,'mu_margin'):
                    self.solverOptions.mu_margin = 1.2 #(1,inf)
                if not hasattr(self.solverOptions,'maxStepSizeTries'):
                    self.solverOptions.maxStepSizeTries = 30
                if not hasattr(self.solverOptions,'watchdogIterations'):
                    self.solverOptions.watchdogIterations = 5
                if not hasattr(self.solverOptions,'stepSchedule'):
                    self.solverOptions.stepSchedule = [] #[0.2,0.5,1.0]
                if not hasattr(self.solverOptions,'baseStepSchedule'):
                    self.solverOptions.baseStepSchedule = [] #[0.2,0.5,1.0]
                if not hasattr(self.solverOptions,'scaleObjective'):
                    self.solverOptions.scaleObjective = False
                if not hasattr(self.solverOptions,'scaleConstraints'):
                    self.solverOptions.scaleConstraints = False
                if not hasattr(self.solverOptions,'scaleVariables'):
                    self.solverOptions.scaleVariables = False
                if not hasattr(self.solverOptions,'relativeTolerance'):
                    self.solverOptions.relativeTolerance = None
            if self.solverOptions.solveType.lower() in ['sp','pccp']:
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'SP_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'SP_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() == 'sqp':
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'SQP_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'SQP_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() == 'lsqp':
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'LSQP_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'LSQP_progress'
                    self.solverOptions.progressFilename = fnm
            if self.solverOptions.solveType.lower() == 'slcp':
                if not hasattr(self.solverOptions,'progressFilename'):
                    cwd = os.getcwd() + os.sep
                    rdName = 'SLCP_runDirectory'
                    if os.path.isdir(cwd + rdName):
                        shutil.rmtree(cwd + rdName)
                    os.mkdir(cwd + rdName, mode=0o777)
                    runDirectory = cwd + rdName + os.sep 
                    fnm = runDirectory + 'SLCP_progress'
                    self.solverOptions.progressFilename = fnm
        else:
            self.solverOptions.solveType = None
# =====================================================================================================================
# =====================================================================================================================