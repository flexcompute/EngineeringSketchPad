import copy
from corsairlite import units, Q_
from corsairlite.optimization import Variable, Constraint, Formulation, Term, Expression, Constant, RuntimeConstraint
from corsairlite.optimization.solve import solve
import numpy as np
from corsairlite.core.dataTypes import OptimizationOutput
import dill
from corsairlite.core.pickleFunctions import setPickle, unsetPickle

from corsairlite.optimization.solverWrappers.cvxopt.sequentialQuadraticProgram import stripUnitsMat, substituteMatrix

import cvxopt

def substituteMatrixTransformed(M, xStar, skipConstants = True):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7

    try:
        len(M[0])
        isMatrix = True
    except:
        isMatrix = False

    if isMatrix:
        mat = []
        for i in range(0,len(M)):
            matRow = []
            for j in range(0,len(M[0])):
                vl = M[i][j]
                if isinstance(vl, (int, float, Q_)):
                    vl *= units.dimensionless
                    matRow.append(vl.to('').magnitude)
                else:
                    vla = vl.substitute(xStar,skipConstants=skipConstants)
                    matRow.append(vla.to('').magnitude)
            mat.append(matRow)
        return mat
    else: # is a Vector
        vec = []
        for i in range(0,len(M)):
            vl = M[i]
            if isinstance(vl, (int, float, Q_)):
                vec.append(vl.to('').magnitude)
            else:
                vla = vl.substitute(xStar,skipConstants=skipConstants)
                vec.append(vla.to('').magnitude)
        return vec

def meritFunction(formulation, stepSize, fullStep_y, ydict, subbedOG, mu_k, lagrangeMultipliers, Ncons, Ncons_explicit, constraint_operators, constraintFunctionEvaluations = None, mode = 1):
    vrs = formulation.variables_only
    Nvars = len(vrs)
    Ncons_form = len(formulation.constraints)
    names = ['y_{' + vrs[i].name + '}' for i in range(0,len(vrs))]
    
    x_k_dict = {}
    i = 0
    for ky,vl in ydict.items():
        x_k_dict[ky[3:-1]] = np.exp(vl + stepSize * fullStep_y[names.index(ky)]) * vrs[i].units
        i += 1
        
    c_1Norm = []

    runConstraintEvaluations = True
    if mode == 1:
        if constraintFunctionEvaluations is not None:
            constraintIndex = -1
            for i in range(0,Ncons_form):
                con = formulation.constraints[i]
                if isinstance(con,Constraint):
                    constraintIndex += 1
                    cval = np.log(constraintFunctionEvaluations[i].magnitude)
                    if constraint_operators[constraintIndex] == '==':
                        c_1Norm.append(abs(cval))
                    elif constraint_operators[constraintIndex] == '>=':
                        c_1Norm.append(abs(min([0,cval])))
                    else:
                        c_1Norm.append(abs(max([0,cval])))
                else: # RuntimeConstraint
                    rs0 = constraintFunctionEvaluations[i][0]
                    for ii in range(0,len(rs0)):
                        constraintIndex += 1
                        cval = np.log(rs0[ii].magnitude)
                        if constraint_operators[constraintIndex] == '==':
                            c_1Norm.append(abs(cval))
                        elif constraint_operators[constraintIndex] == '>=':
                            c_1Norm.append(abs(min([0,cval])))
                        else:
                            c_1Norm.append(abs(max([0,cval])))
                
                runConstraintEvaluations = False
    
    if runConstraintEvaluations:
        constraintFunctionEvaluations = []
        constraintIndex = -1
        for i in range(0,Ncons_form):
            con = formulation.constraints[i]
            if isinstance(con,Constraint):
                constraintIndex += 1
                cval_nolog = con.LHS.substitute(x_k_dict) 
                cval = np.log(cval_nolog.magnitude)
                constraintFunctionEvaluations.append(cval_nolog)
                if constraint_operators[constraintIndex] == '==':
                    c_1Norm.append(abs(cval))
                elif constraint_operators[constraintIndex] == '>=':
                    c_1Norm.append(abs(min([0,cval])))
                else:
                    c_1Norm.append(abs(max([0,cval])))
            else: # RuntimeConstraint
                rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
                RC_xk = [x_k_dict[rcnm] for rcnm in rc_varnames]
                rs = con.returnConstraintSet(RC_xk, 'lsqp')
                constraintFunctionEvaluations.append(rs)
                for ii in range(0,len(rs[0])):
                    constraintIndex += 1
                    cval = np.log(rs[0][ii].magnitude)
                    if constraint_operators[constraintIndex] == '==':
                        c_1Norm.append(abs(cval))
                    elif constraint_operators[constraintIndex] == '>=':
                        c_1Norm.append(abs(min([0,cval])))
                    else:
                        c_1Norm.append(abs(max([0,cval])))

    objectiveGradient = []
    for ogv in subbedOG:
        if hasattr(ogv,'substitute'):
            objectiveGradient.append(ogv.substitute())
        else:
            objectiveGradient.append(ogv)
            
    if mode == 1:
        delfp = sum([objectiveGradient[i] * fullStep_y[i] for i in range(0,Nvars)])
        for i in range(0,Ncons):
            mu_k[i] = max([abs(lagrangeMultipliers[i]), 0.5*(mu_k[i] + abs(lagrangeMultipliers[i]))])
        D = delfp - sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        phi_k = np.log(formulation.objective.substitute(x_k_dict).magnitude) + sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        return phi_k, D, mu_k, constraintFunctionEvaluations
    else:
        delfp = sum([objectiveGradient[i] * fullStep_y[i] for i in range(0,Nvars)])
        D = delfp - sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        phi_k = np.log(formulation.objective.substitute(x_k_dict).magnitude) + sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        return phi_k, D, constraintFunctionEvaluations

def processUnknownResult(sol,newFormulation):
    if newFormulation.solverOptions.solver == 'cvxopt':
        cdns = []
        cdns.append( sol['status']                    == 'unknown' )
        cdns.append( abs(sol['gap'])                  <= 1e-3 )
        cdns.append( abs(sol['relative gap'])         <= 1 )
        cdns.append( abs(sol['primal infeasibility']) <= 1e-3 )
        cdns.append( abs(sol['dual infeasibility'])   <= 1e-3 )
        cdns.append( abs(sol['primal slack'])         <= 1e-3 )
        cdns.append( abs(sol['dual slack'])           <= 1 )

        if all(cdns):
            isValid = True
        else:
            isValid = False
            return sol, isValid

        names = [vbl.name for vbl in newFormulation.variables_only]

        res = OptimizationOutput()
        res.rawResult = copy.deepcopy(sol)
        vrs = list(sol['x'])
        for i in range(0,len(vrs)):
            unt = newFormulation.variables_only[i].units
            vrs[i] = vrs[i] * unt
        obj = sol['primal objective'] * newFormulation.objective.units
        vDict = dict(zip(names, vrs))
        res.variables = vDict
        res.objective = obj
        for tm in newFormulation.objective.terms:
            if isinstance(tm.substitute(), (int,float,Q_)):
                res.objective += tm.substitute()
                
        unsubbedQP = newFormulation.toQuadraticProgram(substituted = False)
        cons   = copy.deepcopy(unsubbedQP.inequalityConstraints) + copy.deepcopy(unsubbedQP.equalityConstraints)
        objUQP = copy.deepcopy(unsubbedQP.objective)
        cons = [objUQP == 0*objUQP.units] + cons
        expressions = [c.LHS for c in cons]
        la = np.array( list(sol['z']) + list(sol['y']) )
        idv = list(sol['z'])
        edv = list(sol['y'])
        icc = 0
        ecc = 0
        dv = []
        for con in unsubbedQP.constraints:
            if con.operator == '==':
                dv.append(edv[ecc])
                ecc += 1
            else:
                dv.append(idv[icc])
                icc +=1
        res.dualVariables = dv

        return res, isValid

def LogspaceSequentialQuadraticProgramSolver(formulation):
    # Setup output
    res = OptimizationOutput()
    res.initialFormulation = copy.deepcopy(formulation)
    res.solveType = formulation.solverOptions.solveType
    res.solver = formulation.solverOptions.solver.lower()
    res.stepNumber = []
    res.designPoints = []
    res.lagrangeMultipliers = []
    res.objectiveFunctionValues = []
    res.intermediateFormulations = []
    res.intermediateSolves = []
    res.error = []
    res.hessianApproximations = []
    res.steps = []
    res.magnitudeOfStep = []
    res.gradientOfLagrangian = []
    res.messages = []
    res.terminationStatus = 'Did not converge, maximum number of iterations reached (%d)'%(formulation.solverOptions.maxIterations)
    res.FOC = []

    # Setup key values used throughout run
    formulation = formulation.toLSQP_StandardForm()
    explicitConstraints = [con for con in formulation.constraints if isinstance(con,Constraint) ]
    runtimeConstraints  = [con for con in formulation.constraints if isinstance(con,RuntimeConstraint)]
    rtcVariables = [ [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs]   for rc in runtimeConstraints]
    allConstraints = explicitConstraints + runtimeConstraints
    formulation.skipEnforceBounds = True
    formulation.constraints = explicitConstraints + runtimeConstraints
    formulation.skipEnforceBounds = False
    formulation.solverOptions = copy.deepcopy(res.initialFormulation.solverOptions)
    vrs = formulation.variables_only
    names = [vrs[i].name for i in range(0,len(vrs))]
    explicitConstraintsOriginal = copy.deepcopy(explicitConstraints)
    runtimeConstraintsOriginal = copy.deepcopy(runtimeConstraints)
    vrsOriginal = copy.deepcopy(vrs)
    namesOriginal = copy.deepcopy(names)
    Nvars = len(vrs)
    Ncons_explicit = len(explicitConstraints)

    # Set the initial guess
    if formulation.solverOptions.x0 is not None:
        guessDict = copy.deepcopy(formulation.solverOptions.x0)
        guessVec = []
        for i in range(0,Nvars):
            v = vrs[i]
            guessVec.append(guessDict[v.name].to(v.units).magnitude * v.units)
        x_k = copy.deepcopy(guessVec)
    else:
        guessVec = []
        for i in range(0,Nvars):
            v = vrs[i]
            guessVec.append(v.guess)
        x_k = copy.deepcopy(guessVec)
        guessDict = dict(zip(names,copy.deepcopy(guessVec)))
        
    # Normalize as desired
    if formulation.solverOptions.scaleObjective:
        objectiveNormalizer = abs(formulation.objective.substitute(guessDict))
        divisor = objectiveNormalizer
        if divisor.magnitude != 0.0:
            formulation.objective = formulation.objective/divisor
        res.messages.append('Objective was scaled by the value at the initial guess')
        
    # Scaling constraints not enabled because the constraints of form g(x) <= 1 and h(x) == 1 are naturally scaled.  This code was not tested
    # if formulation.solverOptions.scaleConstraints:
    #     for i in range(0,len(allConstraints)):
    #         if isinstance(allConstraints[i], Constraint):
    #             divisor = abs(allConstraints[i].LHS.substitute(guessDict))
    #             if divisor.magnitude != 0.0:
    #                 allConstraints[i] = allConstraints[i]/divisor
    #         else: #RuntimeConstraint
    #             allConstraints[i].scaledConstraintVector = [abs(allConstraints[i].outputs[ii].guess) for ii in range(0,len(allConstraints[i].outputs))]
    #     if len(explicitConstraints) > 0:
    #         res.messages.append('Constraints were scaled by the value at the initial guess')
    #     if len(runtimeConstraints) > 0:
    #         res.messages.append('RuntimeConstraints were scaled by the guess provided for the corresponding output')
        
    if formulation.solverOptions.scaleVariables:
        guessList = [formulation.variables_only[i].guess for i in range(0,Nvars)]
        for i in range(0,len(guessList)):
            if guessList[i].magnitude == 0.0:
                guessList[i] = 1.0 * guessList[i].units
        vrsSub = [Variable('v_%d'%(i),1,'') for i in range(0,Nvars)]
        subsList = [guessList[i]*vrsSub[i] for i in range(0,Nvars)]
        subsDict = dict(zip([formulation.variables_only[i].name for i in range(0,Nvars)],[se for se in subsList]))
        formulation.objective = formulation.objective.substitute(subsDict,skipConstants=True)
        for i in range(0,len(allConstraints)):
            if isinstance(allConstraints[i],Constraint):
                allConstraints[i] = allConstraints[i].substitute(subsDict,skipConstants=True)
            else: # RuntimeConstraint
                rc_vars = rtcVariables[i-Ncons_explicit]
                variableIndicies = []
                for nm in rc_vars:
                    variableIndicies.append(names.index(nm))
                allConstraints[i].scaledVariablesVector = [vrs[ii].guess for ii in variableIndicies]
                rcSubsList = [ vrsSub[ii] for ii in variableIndicies]
                rc_subsDict = dict(zip(rc_vars,rcSubsList))
                allConstraints[i] = allConstraints[i].substitute(rc_subsDict)
        x_k = [1.0 * units.dimensionless for i in range(0,Nvars)]
        names = [vrsSub[i].name for i in range(0,len(vrsSub))]
        vrs = vrsSub
        explicitConstraints = [con for con in allConstraints if isinstance(con,Constraint) ]
        runtimeConstraints  = [con for con in allConstraints if isinstance(con,RuntimeConstraint)]
        res.messages.append('Variables were scaled by their provided initial guesses')

    constraint_units = []
    constraint_operators = []
    for i in range(0,len(allConstraints)):
        if isinstance(allConstraints[i], Constraint):
            if allConstraints[i].operator == '>=':
                allConstraints[i] *= -1
            constraint_units.append(allConstraints[i].units)
            constraint_operators.append(allConstraints[i].operator)
        elif isinstance(allConstraints[i], RuntimeConstraint):
            for ii in range(0,len(allConstraints[i].operators)):
                if formulation.solverOptions.scaleConstraints:
                    constraint_units.append(units.dimensionless)
                else:
                    if formulation.solverOptions.scaleVariables:
                        constraint_units.append(runtimeConstraintsOriginal[i-Ncons_explicit].outputs[ii].units)
                    else:
                        constraint_units.append(allConstraints[i].outputs[ii].units)
                if allConstraints[i].operators[ii] == '==':
                    constraint_operators.append('==')
                else:
                    constraint_operators.append('<=')
        else:
            raise ValueError('Invalid constraint type')
    
    Ncons = len(constraint_operators)
        
    # Set up the transformed variables
    yVariableList = []
    yNames = []
    for i in range(0,Nvars):
        v = vrs[i]
        nm = 'y_{%s}'%(v.name)
        yNames.append(nm)
        y = Variable(name = nm, 
                     guess = 1.0, 
                     units = '',#v.units,   
                     description = "Transformed variable for %s"%(v.name))
        yVariableList.append(y)

    # Set up the slack variables
    slackVariableList = []
    slackNames = []
    slackCounter = 0
    for i in range(0,Ncons):
        slackNames.append("__sqpSlack_%d"%(slackCounter))
        s = Variable(name = "__sqpSlack_%d"%(slackCounter), 
                     guess = 1.0, 
                     units = "-",   
                     description = "Slack number %d for LSQP"%(slackCounter))
        slackCounter += 1
        slackVariableList.append(s)

    # Set up the step variables
    stepVariableList = []
    stepNames = []
    for i in range(0,Nvars):
        v = vrs[i]
        nm = 'd_{%s}'%(v.name)
        stepNames.append(nm)
        s = Variable(name = nm, 
                     guess = 1.0, 
                     units = '',#v.units,   
                     description = "Step variable for log of %s"%(v.name))
        stepVariableList.append(s)

    # Set up the lagrange multiplier variables
    lamv = []
    lamNames = []
    lamCounter = 0
    for i in range(0,Ncons):
        lamNames.append("__lambda_%d"%(lamCounter))
        lu = formulation.objective.units/constraint_units[i]
        lm = Variable(name = "__lambda_%d"%(lamCounter), 
                     guess = 1.0, 
                     units = "%s"%(str(lu)),   
                     description = "Dual Variable %d for LSQP"%(lamCounter))
        lamCounter += 1
        lamv.append(lm)

    # Initalize the lagrange mutipliers
    lam = [1.0 for i in range(0,Ncons)]

    # Set up hesian approximation (initalizes to identity)
    Bmat = np.eye(Nvars).tolist()
    res.hessianApproximations.append(Bmat)

    # Determine Gradients
    objectiveGradient = [formulation.objective.derivative(v) for v in vrs]
    constraintGradients = [[allConstraints[i].LHS.derivative(v) for v in vrs] for i in range(0,Ncons_explicit)] # only holds values for explicit constraints

    # Set up the lagrangian
    Lagrangian_us = copy.deepcopy(formulation.objective)
    expr = Expression()
    prodList = np.array(lamv[0:Ncons_explicit]) * np.array([allConstraints[i].LHS - 1.0 * allConstraints[i].units for i in range(0,Ncons_explicit)])
    for i in range(0,len(prodList)):
        expr.terms.extend(prodList[i].terms)
    expr.simplify()
    Lagrangian_us += expr

    Lagrangian_us_gradient = [None] * len(vrs)
    for i in range(0,len(vrs)):
        Lagrangian_us_gradient[i] = objectiveGradient[i]
        expr = Expression()
        prodList = [lamv[iii] * [constraintGradients[ii][i] for ii in range(0,Ncons_explicit)][iii] for iii in range(0,Ncons_explicit)]
        for ii in range(0,len(prodList)):
            if isinstance(prodList[ii],Expression):
                expr.terms.extend(prodList[ii].terms)
            elif isinstance(prodList[ii],Term):
                expr.terms.extend([prodList[ii]])
        expr.simplify()
        Lagrangian_us_gradient[i] += expr


    # Compute initial values
    x_k_dict = dict(zip(names,x_k))
    objectiveValue = formulation.objective.substitute(x_k_dict)
    res.objectiveFunctionValues.append(objectiveValue)

    # Main Iterative Loop
    watchdogCounter = 0
    stepNumberIndex = -1
    constraintFunctionEvaluations = None
    for iterate in range(0,formulation.solverOptions.maxIterations):
        stepNumberIndex += 1
        res.stepNumber.append(stepNumberIndex)
        lamSubs = [lam[i] * lamv[i].units for i in range(0,Ncons)]
        lamDict = dict(zip(lamNames,lamSubs))
        
        x_k_dict = dict(zip(names,x_k))
        subbedOH = Bmat
        if iterate != 0:
            objectiveValue = objectiveValue_kp1
        objectiveGradient_y = [ x_k[i]/objectiveValue * objectiveGradient[i] for i in range(0,Nvars) ]
        subbedOG = substituteMatrixTransformed(objectiveGradient_y, x_k_dict)
        y_k = []
        for i in range(0,len(x_k)):
            y_k.append(np.log(x_k[i].to(vrs[i].units).magnitude))
        y_k_dict = dict(zip(yNames,y_k))

        # =====================================================
        # =====================================================
        # =====================================================
        Nvars_qp = Ncons+Nvars
        Praw = np.zeros([Nvars_qp,Nvars_qp])
        qraw = np.zeros([Nvars_qp,1])
        
        np.fill_diagonal(Praw, formulation.solverOptions.penaltyConstant)
        Praw[-Nvars:,-Nvars:] = stripUnitsMat(subbedOH)
        qraw[-Nvars:,:] = np.array(subbedOG).reshape([Nvars,1])
        P = cvxopt.matrix(Praw)
        q = cvxopt.matrix(qraw)

        crv = []
        subbedCG = []
        if constraintFunctionEvaluations is None: # should only activate on the first iteration
            constraintGradients_subbed = substituteMatrix(constraintGradients,x_k_dict,False)
            constraintFunctionEvaluations = []
            for i in range(0,Ncons_explicit):
                vl = allConstraints[i].LHS.substitute(x_k_dict)
                crv.append(-1*np.log(vl.magnitude))
                constraintFunctionEvaluations.append(vl)
                subbedCG.append( [ x_k[ii] / vl * constraintGradients_subbed[i][ii] for ii in range(0,len(constraintGradients_subbed[i])) ] )
            for i in range(0,len(runtimeConstraints)):
                rc = runtimeConstraints[i]
                rc_varnames = [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs] 
                RC_xk = [x_k_dict[rcnm] for rcnm in rc_varnames]
                rc_rs = rc.returnConstraintSet(RC_xk, 'lsqp')
                constraintFunctionEvaluations.append(rc_rs)
                rc_eval = rc_rs[0]
                rc_grad = rc_rs[1]
                for sii in range(0,len(rc_eval)):
                    vl = rc_eval[sii]
                    crv.append(-1*np.log(vl.magnitude))
                for sii in range(0,len(rc_grad)):
                    vl = rc_eval[sii]
                    gd = rc_grad[sii]
                    rc_grad_corrected = []
                    obj_uts = formulation.objective.units
                    for ix in range(0,len(names)):
                        nm = names[ix]
                        uts = vrs[ix].units
                        if nm in rc_varnames:
                            rc_grad_corrected.append(gd[rc_varnames.index(nm)])
                        else:
                            rc_grad_corrected.append(0.0*obj_uts/uts)
                    subbedCG.append( [ x_k[ii] / vl * rc_grad_corrected[ii] for ii in range(0,len(rc_grad_corrected)) ] )
        else: # Function evaluations are made in the compution of step in the previous iteration
            cfe_old = copy.deepcopy(constraintFunctionEvaluations)
            constraintFunctionEvaluations = []
            for i in range(0,Ncons_explicit):
                vl = cfe_old[i]
                crv.append(-1*np.log(vl.magnitude))
                constraintFunctionEvaluations.append(vl)
            for i in range(Ncons_explicit,Ncons_explicit+len(runtimeConstraints)):
                constraintFunctionEvaluations.append(cfe_old[i])
                rs0 = cfe_old[i][0]
                rs1 = cfe_old[i][1]
                for sii in range(0,len(rs0)):
                    vl = rs0[sii]
                    crv.append(-1*np.log(vl.magnitude))
            # subbedCG comes from pervious iteration directly
            if iterate != 0:
                subbedCG = subbedCG_kp1

        crv = np.reshape(np.array(crv),[Ncons,1])
        
        subbedCG = stripUnitsMat(subbedCG)
        
        cgm = np.diag(-1*np.ones(Ncons))
        cgm = np.hstack((cgm,subbedCG))
        
        Graw = None
        Araw = None
        hraw = None
        braw = None
        
        for conNum in range(0,len(constraint_operators)):
            conop = constraint_operators[conNum]
            if conop == '==':
                if Araw is None:
                    Araw = []
                    braw = []
                Araw.append(cgm[conNum])
                braw.append(crv[conNum])
            else:
                if Graw is None:
                    Graw = []
                    hraw = []
                Graw.append(cgm[conNum])
                hraw.append(crv[conNum])   

        Graw = (np.array(Graw).T).tolist()
        Araw = (np.array(Araw).T).tolist()
        hraw = (np.array(hraw).T).tolist()
        braw = (np.array(braw).T).tolist()

        if Graw is not None:
            G = cvxopt.matrix(Graw)
        else:
            G = None
        if hraw is not None:
            h = cvxopt.matrix(hraw)
        else:
            h = None
        if Araw is not None:
            A = cvxopt.matrix(Araw)
        else:
            A = None
        if braw is not None:
            b = cvxopt.matrix(braw)
        else:
            b = None
            
        if A is not None and G is not None:
            sol = cvxopt.solvers.qp( P, q, G, h, A, b)
        elif A is not None:
            sol = cvxopt.solvers.qp( P, q, A = A, b = b)
        elif G is not None:
            sol = cvxopt.solvers.qp( P, q, G, h)
        else:
            sol = cvxopt.solvers.qp( P, q)

        idv = list(sol['z'])
        edv = list(sol['y'])
        icc = 0
        ecc = 0
        dv = []
        for conop in constraint_operators:
            if conop == '==':
                dv.append(edv[ecc])
                ecc += 1
            else:
                dv.append(idv[icc])
                icc +=1
        QP_dualVariables = dv
        QP_slacks = abs(sol['x'][0:Ncons])
        QP_stepVars = sol['x'][-Nvars:]
#        # print(QP_stepVars)
#        # print(QP_dualVariables)
        # cody

        # =====================================================
        # =====================================================
        # =====================================================
        for watchdog_Iterate in [0,1]:
            slacks = QP_slacks
            flps = []
            for conop in constraint_operators:
                if conop == '<=' or conop == '==':
                    flps.append(1.0)
                else:
                    flps.append(-1.0)

            newLagrangeMultipliers = [flps[i]*QP_dualVariables[i] for i in range(0,Ncons)]
            fullStep_l = [newLagrangeMultipliers[i] - lam[i] for i in range(0,Ncons)]
            fullStep_y = QP_stepVars
            y_kp1_full = [y_k[i] + fullStep_y[i] for i in range(0,Nvars)]
            x_kp1_full = [np.exp(y_kp1_full[i]) * vrs[i].units for i in range(0,Nvars)] 
            fullStep_x = [x_kp1_full[i] - x_k[i] for i in range(0,Nvars)]
            
            if stepNumberIndex in range(0,len(formulation.solverOptions.stepSchedule)):
                stpLm = formulation.solverOptions.stepSchedule[stepNumberIndex]
                maxAllowedStep = min(stpLm * abs(np.array([(x_k[idx]/fullStep_x[idx]).to('').magnitude for idx in range(0,Nvars)])))
            else:
                maxAllowedStep = 1.0

            if iterate == 0:
                mu_k = np.ones(Ncons).tolist()

            eta = formulation.solverOptions.eta
            zeta = formulation.solverOptions.zeta
            tau = formulation.solverOptions.tau
            stepSize = 1.0

            # nlmu = [newLagrangeMultipliers[i] * lamv[i].units for i in range(0,Ncons)]
            cfe_saveForCache = copy.deepcopy(constraintFunctionEvaluations)
            
            MF_output = meritFunction(formulation, 0.0, fullStep_y, y_k_dict, subbedOG, mu_k, newLagrangeMultipliers, Ncons, Ncons_explicit, constraint_operators, constraintFunctionEvaluations, 1)
            phi_k = MF_output[0]
            D     = MF_output[1]
            mu_k  = MF_output[2]
            constraintFunctionEvaluations_previousIteration = MF_output[3]

            for i in range(0,formulation.solverOptions.maxStepSizeTries):
                doWolfeChecking = True
                try:
                    phi_kp1, D_kp1, constraintFunctionEvaluations = meritFunction(formulation, stepSize, fullStep_y, y_k_dict, subbedOG, mu_k, newLagrangeMultipliers, Ncons, Ncons_explicit, constraint_operators, None, 2)
                except:
                    # Step was too big, one of the constraints became imaginary
                    doWolfeChecking = False
                    res.messages.append('Attempted step of size %.2e at iteration %d is too large and is causing a constraint to return as imaginary.  Decreasing according to the tau schedule.'%(
                                        stepSize, stepNumberIndex))
                    stepSize *= tau
                    noRerun = True

                if doWolfeChecking:
                    wolfe_1 = phi_kp1 <= phi_k + eta * stepSize * D
                    if wolfe_1:
                        watchdogCounter = 0
                        checkpointStep = stepNumberIndex
                        noRerun = True
                        break
                    else:
                        watchdogCounter += 1
                        if watchdogCounter == 1:
                            checkpointStep = stepNumberIndex
                            y_k_cache = copy.deepcopy(y_k)
                            lam_cache = copy.deepcopy(lam)
                            Bmat_cache = copy.deepcopy(Bmat)
                            QP_dualVariables_cache = copy.deepcopy(QP_dualVariables)
                            QP_slacks_cache = copy.deepcopy(QP_slacks)
                            QP_stepVars_cache = copy.deepcopy(QP_stepVars)
                            subbedOG_cache = copy.deepcopy(subbedOG)
                            constraintFunctionEvaluations_cache = copy.deepcopy(cfe_saveForCache)
                            noRerun = True
                            res.messages.append('Wolfe Condition 1 failed at iteration %d (Step %d).  Watchdog counter is %d'%(iterate,stepNumberIndex,watchdogCounter))
                            break
                        elif watchdogCounter > formulation.solverOptions.watchdogIterations:
                            if watchdog_Iterate == 0:
                                res.messages.append('Wolfe Condition 1 failed at iteration %d (Step %d).  '%(iterate,stepNumberIndex) + 
                                                     'Resetting to Step %d.'%(checkpointStep))
                                noRerun = False
                                stepNumberIndex = checkpointStep
                                y_k = y_k_cache
                                lam = lam_cache
                                Bmat = Bmat_cache
                                QP_dualVariables = QP_dualVariables_cache
                                QP_slacks = QP_slacks_cache
                                QP_stepVars = QP_stepVars_cache
                                constraintFunctionEvaluations = constraintFunctionEvaluations_cache
                                y_k_dict = dict(zip(yNames,y_k))
                                subbedOG = subbedOG_cache
                                break
                            else:
                                stepSize *= tau
                        else:
                            #watchdog requires no action yet, accept full step
                            res.messages.append('Wolfe Condition 1 failed at iteration %d (Step %d).  Watchdog counter is %d'%(iterate,stepNumberIndex,watchdogCounter))
                            noRerun = True
                            break

            if noRerun:
                break
                
        if stepSize > maxAllowedStep:
            res.messages.append('Step %d (%.2e) is too large given the uncertianty in the curvature matrix.  Decreased to a smaller step (%.2e)'%(
                                       stepNumberIndex, stepSize, maxAllowedStep))
            stepSize = maxAllowedStep
            # Not taking the full step called for by the algorithm, must rerun the constraint functions at call
            constraintFunctionEvaluations = []
            tmp_x_k_dict = dict(zip(names,[x_k[i]  + stepSize * fullStep_x[i] for i in range(0,Nvars)]))
            for i in range(0,len(allConstraints)):
                con = allConstraints[i]
                if isinstance(con,Constraint):
                    cval = con.LHS.substitute(tmp_x_k_dict) 
                    constraintFunctionEvaluations.append(cval)
                else: # RuntimeConstraint
                    rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
                    RC_xk = [tmp_x_k_dict[rcnm] for rcnm in rc_varnames]
                    rs = con.returnConstraintSet(RC_xk, 'lsqp')
                    constraintFunctionEvaluations.append(rs)
        
        if watchdogCounter > formulation.solverOptions.watchdogIterations:
            watchdogCounter = 0

        res.intermediateSolves.append(sol)
        res.steps.append(stepSize)
        y_kp1   = [y_k[i]  + stepSize * fullStep_y[i] for i in range(0,Nvars)]
        x_kp1   = [x_k[i]  + stepSize * fullStep_x[i] for i in range(0,Nvars)]
        lam_kp1 = [lam[i]  + stepSize * fullStep_l[i] for i in range(0,Ncons)]
        x_kp1_dict = dict(zip(names,x_kp1))
        
        newObj = formulation.objective.substitute(x_kp1_dict)
        res.objectiveFunctionValues.append(newObj)
        if iterate == 0:
            res.designPoints.append(x_k)     
        res.designPoints.append(x_kp1)
        
        objectiveValue_kp1 = formulation.objective.substitute(x_kp1_dict)
        objectiveGradient_y_kp1 = [ x_kp1[i]/objectiveValue_kp1 * objectiveGradient[i] for i in range(0,Nvars) ]
        subbedOG_kp1 = substituteMatrixTransformed(objectiveGradient_y_kp1, x_kp1_dict)
        
        lagrangianGradient_k = np.array([subbedOG[j] +  sum([lam_kp1[i] * subbedCG[i][j] for i in range(0,Ncons)]) for j in range(0,Nvars)])
        constraintValues_kp1 = []
        for i in range(0,Ncons_explicit):
            constraintValues_kp1.append(constraintFunctionEvaluations[i])
        for i in range(0,len(runtimeConstraints)):
            ix = i + Ncons_explicit
            for ii in range(0,len(runtimeConstraints[i].outputs)):
                constraintValues_kp1.append(constraintFunctionEvaluations[ix][0][ii])
                
        lam_x_kp1 = [objectiveValue_kp1 / constraintValues_kp1[i] * lam_kp1[i] for i in range(0,Ncons)]
        lamDict = dict(zip(lamNames,lam_x_kp1))
        
        subbedCG_kp1 = []
        constraintGradients_subbed = substituteMatrix(constraintGradients,x_kp1_dict,False)
        for i in range(0,Ncons_explicit):
            vl = constraintFunctionEvaluations[i]
            subbedCG_kp1.append( [ x_kp1[ii] / vl * constraintGradients_subbed[i][ii] for ii in range(0,len(constraintGradients_subbed[i])) ] )
        for i in range(Ncons_explicit,Ncons_explicit+len(runtimeConstraints)):
            rs0 = constraintFunctionEvaluations[i][0]
            rs1 = constraintFunctionEvaluations[i][1]
            rtcIdx = i-Ncons_explicit
            rc_varnames = rtcVariables[i-Ncons_explicit]
            for sii in range(0,len(rs1)):
                rc_varnames = [ip.name for ip in runtimeConstraints[rtcIdx].inputs] + [runtimeConstraints[rtcIdx].outputs[sii].name] 
                vl = rs0[sii]
                gd = rs1[sii]
                rc_grad_corrected = []
                for ix in range(0,len(names)):
                    nm = names[ix]
                    uts = vrs[ix].units
                    if nm in rc_varnames:
                        rc_grad_corrected.append(gd[rc_varnames.index(nm)])
                    else:
                        rc_grad_corrected.append(0.0*vl.units/uts)
                subbedCG_kp1.append( [ x_kp1[ii] / vl * rc_grad_corrected[ii] for ii in range(0,len(rc_grad_corrected)) ] )

        lagrangianGradient_kp1_withUnits = [subbedOG_kp1[j] + sum([lam_kp1[i] * subbedCG_kp1[i][j] for i in range(0,Ncons)]) for j in range(0,Nvars)]
        lagrangianGradient_kp1 = np.array([lagrangianGradient_kp1_withUnits[i].to('').magnitude for i in range(0,len(lagrangianGradient_kp1_withUnits))])

        FOC = [] # First order conditions (grad of lagrangian in untransformed space).  Computation allows for common comparison with SQP using the same termination condition
                 # Generally should not use in practice, as the ordinary logspace termination condition is faster and of sufficient tolerance
        lagrangianGradient_untransformed = [Lagrangian_us_gradient[vr_ix].substitute(lamDict) for vr_ix in range(0,len(vrs))] # The version with only explicit constraints, runtime added on below
        rtc_LGk = [0.0*formulation.objective.units / v.units for v in vrs]
        rtc_LGkp1 = [0.0*formulation.objective.units / v.units for v in vrs]
        rtcEvaluations = constraintFunctionEvaluations[Ncons_explicit:]
        rtc_vrs = [ [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs]   for rc in runtimeConstraints]
        for i in range(0,Nvars):
            nm = names[i]
            constraintNum = Ncons_explicit - 1
            for ii in range(0,len(rtcVariables)):
                constraintNum += 1
                rcNames = rtc_vrs[ii]
                if nm in rcNames:
                    # I think there is an indexing issue here for multi output
                    dc_dxkp1 = constraintFunctionEvaluations[constraintNum][1][0][rcNames.index(nm)]
                    current_lam = lam_x_kp1[constraintNum]
                    rtc_LGkp1[i] += current_lam * dc_dxkp1
        for i in range(0,len(lagrangianGradient_untransformed)):
            lg = lagrangianGradient_untransformed[i]
            if hasattr(lg,'substitute'):
                FOC.append(lg.substitute(dict(zip(names,x_kp1))) + rtc_LGkp1[i])
            else:
                FOC.append(lg+rtc_LGkp1[i])     
        FOC = np.array([FOC[i].magnitude for i in range(0,len(FOC))])
        
        s_k    = [stepSize * fullStep_y[i] for i in range(0,Nvars)]
        gl_k   = lagrangianGradient_k
        gl_kp1 = lagrangianGradient_kp1
        magStep = sum([(x_kp1[i].magnitude - x_k[i].magnitude)**2 for i in range(0,Nvars)])**0.5 
        gradLag = gl_kp1
        res.magnitudeOfStep.append(magStep)
        res.gradientOfLagrangian.append(gradLag)
        res.FOC.append(FOC)
        
        tc1 = magStep < formulation.solverOptions.stepMagnitudeTolerance
        tc2 = all(abs(gl_kp1) < formulation.solverOptions.lagrangianGradientTolerance)
        # tc2 = all(abs(FOC) < formulation.solverOptions.lagrangianGradientTolerance)
        tc3 = all(np.array(res.magnitudeOfStep[-5:]) <= formulation.solverOptions.stepMagnitudeTolerance)

        if tc2:
            res.terminationStatus = ('Converged, gradient of Lagrangian (%e) within tolerance (%e) '%( max(abs(gl_kp1)),
                                               formulation.solverOptions.lagrangianGradientTolerance) )
            break
        elif tc3:
            res.terminationStatus = ('Converged, magnitude of step (%.2e)'%(res.magnitudeOfStep[-1]) +
                                     ' is witin specified tolerance (%.2e)'%(formulation.solverOptions.stepMagnitudeTolerance))
            break
        else:
            #Update
            y_k = [gl_kp1[i] - gl_k[i] for i in range(0,len(gl_k))]
            sTy = sum([s_k[i] * y_k[i] for i in range(0,len(s_k))])
            Bs = [sum([Bmat[i][j] * s_k[j] for j in range(0,len(s_k))]) for i in range(0,len(Bmat))]
            sTBs = sum([s_k[i] * Bs[i] for i in range(0,len(s_k))])
            thetaCon = sTy >= 0.2*sTBs
            if thetaCon:
                theta = 1.0 * units.dimensionless
            else:
                theta = 0.8*sTBs / (sTBs - sTy)

            r_k = [theta * y_k[i] + (1-theta)*Bs[i] for i in range(0,len(y_k))]
            sTr = sum([s_k[i] * r_k[i] for i in range(0,len(s_k))])
            rrT_sTr = [[r_k[i] * r_k[j]/sTr for i in range(0,len(r_k))] for j in range(0,len(r_k))]

            sTB = [sum([s_k[i] * Bmat[i][j] for i in range(0,len(s_k))]) for j in range(0,len(Bmat))]
            BssTB_sTBs = [[Bs[i] * sTB[j]/sTBs for i in range(0,len(Bs))] for j in range(0,len(sTB))]

            Bmat_new = [[Bmat[i][j] - BssTB_sTBs[i][j] + rrT_sTr[i][j] for i in range(0,len(Bmat)) ] for j in range(0,len(Bmat[0]))]
            Bmat = Bmat_new
            res.hessianApproximations.append(Bmat)

            # for rw in Bmat_new:
            #     print(stripUnitsMat(rw))

            x_k = x_kp1
            lam = lam_kp1

            if formulation.solverOptions.progressFilename is not None:
                setPickle()
                with open(formulation.solverOptions.progressFilename + "_iteration-%d.dl"%(iterate+1), "wb") as dill_file:
                    dill.dump(res, dill_file)
                with open(formulation.solverOptions.progressFilename, "wb") as dill_file:
                    dill.dump(res, dill_file)
                unsetPickle()

    vDict = dict(zip(names, x_kp1))
    res.solvedFormulation = copy.deepcopy(formulation)
    if formulation.solverOptions.scaleVariables:
        finalVariables = {}
        ix_ctr = 0
        vrsOriginal = res.initialFormulation.variables_only
        namesOriginal = [vrsOriginal[i].name for i in range(0,len(vrsOriginal))]
        for nm in namesOriginal:
            mulFac = vDict['v_%d'%(ix_ctr)]
            finalVariables[nm] = mulFac * guessDict[nm]
            ix_ctr += 1
        res.variables = finalVariables
    else:
        res.variables = vDict
    if formulation.solverOptions.scaleObjective:
        res.objective = newObj*objectiveNormalizer
    else:
        res.objective = newObj
    res.numberOfIterations = iterate + 1
    #     res.sensitivities = res.intermediateSolves[-1].sensitivities

    for i in range(0,len(slacks)):
        if abs(slacks[i]) > formulation.solverOptions.constraintTolerance:
            res.messages.append('WARNING: One or more constraint tolerances not met at the final solution (tol=%.2e)'%(formulation.solverOptions.constraintTolerance))

    return res