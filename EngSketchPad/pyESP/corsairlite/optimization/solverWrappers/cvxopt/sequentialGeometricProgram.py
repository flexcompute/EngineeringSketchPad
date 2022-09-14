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

def meritFunction(  vrs,
                    rho,
                    mu_margin,
                    objective,
                    constraints, 
                    stepSize, 
                    fullStep_x, 
                    xdict, 
                    subbedOG, 
                    mu_k, 
                    lagrangeMultipliersRaw, 
                    constraint_operators, 
                    constraintFunctionEvaluations = None, 
                    mode = 1):
    sigma = 1.0
    
    Nvars = len(vrs)
    names = [vrs[i].name for i in range(0,len(vrs))]
    
    x_k_dict = {}
    for ky,vl in xdict.items():
        x_k_dict[ky] = vl + stepSize * fullStep_x[names.index(ky)]
    
    c_1Norm = []
    constraintValues = []
    
    def returncval(vl_fcn,op):
        vl_fcn -= 1.0*vl_fcn.units
        if op == '==':
            return abs(vl_fcn)
        elif op == '>=':
            return abs(min([0*vl_fcn.units,vl_fcn]))
        else:
            return abs(max([0*vl_fcn.units,vl_fcn]))
    
    runConstraintEvaluations = True
    if mode == 1:
        constraintIndex = -1
        if constraintFunctionEvaluations is not None:
            for jj in [0,1]:
                for i in range(0,len(constraintFunctionEvaluations[jj])):
                    constraintIndex += 1
                    cval = constraintFunctionEvaluations[jj][i]
                    constraintValues.append(cval)
                    c_1Norm.append(returncval(cval,constraint_operators[constraintIndex]))
            for i in range(0,len(constraintFunctionEvaluations[2])):
                rs0 = constraintFunctionEvaluations[2][i]['con'][0]
                for ii in range(0,len(rs0)):
                    constraintIndex += 1
                    cval = rs0[ii]
                    constraintValues.append(cval)
                    c_1Norm.append(returncval(cval,constraint_operators[constraintIndex]))
                
        runConstraintEvaluations = False
    
    if runConstraintEvaluations:
        constraintFunctionEvaluations = [[],[],[]]
        constraintIndex = -1
        for con in constraints[0]:
            constraintIndex += 1
            vl = con.LHS.substitute(x_k_dict)
            constraintFunctionEvaluations[0].append(vl)
            constraintValues.append(vl)
            c_1Norm.append(returncval(vl,constraint_operators[constraintIndex]))
        for fr in constraints[1]:
            constraintIndex += 1
            vl = fr.substitute(x_k_dict)
            constraintFunctionEvaluations[1].append(vl)
            constraintValues.append(vl)
            c_1Norm.append(returncval(vl,constraint_operators[constraintIndex]))
        for con in constraints[2]:
            rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
            RC_xk = [x_k_dict[rcnm] for rcnm in rc_varnames]
            rs = con.returnConstraintSet(RC_xk, 'sgp') 
            constraintFunctionEvaluations[2].append(rs)
            for ii in range(0,len(rs['con'][0])):
                constraintIndex += 1
                cval = rs['con'][0][ii]
                constraintValues.append(cval)
                c_1Norm.append(returncval(cval,constraint_operators[constraintIndex]))

    objectiveValue = objective.substitute(x_k_dict)
    lagrangeMultipliers = [lagrangeMultipliersRaw[i]*objectiveValue/constraintValues[i] for i in range(0,len(lagrangeMultipliersRaw))] 
#     lagrangeMultipliers = [lagrangeMultipliersRaw[i]/constraintValues[i] for i in range(0,len(lagrangeMultipliersRaw))] 
        
    objectiveGradient = []
    for ogv in subbedOG:
        if hasattr(ogv,'substitute'):
            objectiveGradient.append(ogv.substitute())
        else:
            objectiveGradient.append(ogv)
            
    Ncons = len(lagrangeMultipliers)
                
    if mode == 1:
        delfp = sum([objectiveGradient[i] * fullStep_x[i] for i in range(0,Nvars)])
        for i in range(0,Ncons):
            mu_k[i] = max([abs(lagrangeMultipliers[i]), 0.5*(mu_k[i] + abs(lagrangeMultipliers[i]))])
        D = delfp - sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        phi_k = objective.substitute(x_k_dict) + sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        return phi_k, D, mu_k, constraintFunctionEvaluations, lagrangeMultipliers
    else:
        delfp = sum([objectiveGradient[i] * fullStep_x[i] for i in range(0,Nvars)])
        D = delfp - sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        phi_k = objective.substitute(x_k_dict) + sum([mu_k[i] * c_1Norm[i] for i in range(0,Ncons)])
        return phi_k, D, constraintFunctionEvaluations, lagrangeMultipliers




def SequentialGeometricProgramSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7
    
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
    res.hessianApproximations = []
    res.steps = []
    res.magnitudeOfStep = []
    res.gradientOfLagrangian = []
    res.messages = []
    res.terminationStatus = 'Did not converge, maximum number of iterations reached (%d)'%(formulation.solverOptions.maxIterations)
    res.terminationFlag = False

    # Setup key values used throughout run
    gpcConstraints = []
    ngpcConstraints = []
    runtimeConstraints = []
    for con in formulation.constraints:
        if isinstance(con,Constraint):
            try: # is gp compatible
                gpConstraintSet = con.toPosynomial(GPCompatible = True)
                gpcConstraints.extend(gpConstraintSet)
            except: # is explicit but not gp compatible
                zeroedCon = con.toZeroed()
                if zeroedCon.operator == '==':
                    con1 = Constraint(zeroedCon.LHS,'<=',0.0*zeroedCon.units)
                    con2 = Constraint(zeroedCon.LHS,'>=',0.0*zeroedCon.units)
                    con2 *= -1
                    ngpcConstraints.append(con1)
                    ngpcConstraints.append(con2)
                elif zeroedCon.operator == '>=':
                    zeroedCon *= -1
                    ngpcConstraints.append(zeroedCon)
                else:
                    ngpcConstraints.append(zeroedCon)
        else: #runtime constraint
            runtimeConstraints.append(copy.deepcopy(con))
    
    # This is temporary until the fraction class is fully integrated
    for i in range(0,len(gpcConstraints)):
        gpcConstraints[i] = gpcConstraints[i].substitute()
    for i in range(0,len(ngpcConstraints)):
        ngpcConstraints[i] = ngpcConstraints[i].substitute()
    
    allConstraints = gpcConstraints + ngpcConstraints + runtimeConstraints
    formulation.constraints = copy.deepcopy(allConstraints)
    ngpcConstraintsLHS = []
    for i in range(0,len(ngpcConstraints)):
        con = ngpcConstraints[i]
        posCon = con.toAllPositive()
        if posCon.operator == '>=':
            posCon.swap()
        frac = posCon.LHS/posCon.RHS
        ngpcConstraintsLHS.append(frac)
    
    allConstraints = gpcConstraints + ngpcConstraints + runtimeConstraints
    explicitConstraints = gpcConstraints + ngpcConstraints

    rtcVariables = [ [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs]   for rc in runtimeConstraints]
    formulation.solverOptions = copy.deepcopy(res.initialFormulation.solverOptions)
    vrs = formulation.variables_only
    names = [vrs[i].name for i in range(0,len(vrs))]
    Nvars = len(vrs)
    Ncons_explicit = len(gpcConstraints + ngpcConstraints)
    
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
                        constraint_units.append(units.dimensionless)
                if allConstraints[i].operators[ii] == '==':
                    constraint_operators.append('==')
                else:
                    constraint_operators.append('<=')
        else:
            raise ValueError('Invalid constraint type')
    
    Ncons = len(constraint_operators)

    # Set up the slack variables
    slackVariableList = []
    slackNames = []
    slackCounter = 0
    for i in range(0,Ncons):
        slackNames.append("__sgpSlack_%d"%(slackCounter))
        s = Variable(name = "__sgpSlack_%d"%(slackCounter), 
                     guess = 1.0, 
                     units = "-",   
                     description = "Slack number %d for SGP"%(slackCounter))
        slackCounter += 1
        slackVariableList.append(s)

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
                         description = "Dual Variable %d for SGP"%(lamCounter))
            lamCounter += 1
            lamv.append(lm)

    # Initalize the lagrange mutipliers
    lam = [1.0 for i in range(0,Ncons)]

    # Set the hessian matricies for the runtime constraints
    runtimeConstraintHessians = []
    for rc in runtimeConstraints:
        hessianStack = []
        for op in rc.outputs:
            Nvars_rtc = len(rc.inputs)
            Bmat = np.eye(Nvars_rtc).tolist()
            # Bmat = [[0,1,1],[1,0,0],[1,0,0]]
#             Bmat = [[1,1,1],[1,1,1],[1,1,1]]
            Bmat_withUnits = copy.deepcopy(Bmat)
            Nvars_rtc = len(rc.inputs)
            for i in range(0,Nvars_rtc):
                for j in range(0,Nvars_rtc):
                    nm1 = rc.inputs[i].name
                    nm2 = rc.inputs[j].name
                    Bmat_withUnits[i][j] *= op.units / guessDict[nm1].units / guessDict[nm2].units
            hessianStack.append(Bmat_withUnits)
        runtimeConstraintHessians.append(hessianStack)
    
    # Determine Gradients
    objectiveGradient = [formulation.objective.derivative(v) for v in vrs]
    #This is temporary until the fraction class gets fully integrated
    constraintGradients = [[gpcConstraints[i].LHS.derivative(v) for v in vrs] for i in range(0,len(gpcConstraints))] # only the gp compatible constraints
    constraintGradients.extend([[ngpcConstraintsLHS[i].derivative(v) for v in vrs] for i in range(0,len(ngpcConstraintsLHS))]) # does the rest of the explicit constraints
        
#     # Set up the lagrangian, only uses explicit constraints
#     Lagrangian_us = copy.deepcopy(formulation.objective)
#     expr = Expression()
#     prodList = np.array(lamv[0:Ncons_explicit]) * np.array([formulation.constraints[i].LHS for i in range(0,Ncons_explicit)])
#     for i in range(0,len(prodList)):
#         expr.terms.extend(prodList[i].terms)
#     expr.simplify()
#     Lagrangian_us += expr
    
#     Lagrangian_us_gradient = [None] * len(vrs)
#     for i in range(0,len(vrs)):
#         Lagrangian_us_gradient[i] = objectiveGradient[i]
#         expr = Expression()
#         prodList = [lamv[iii] * [constraintGradients[ii][i] for ii in range(0,Ncons_explicit)][iii] for iii in range(0,Ncons_explicit)]
#         print(prodList)
#         for ii in range(0,len(prodList)):
#             if isinstance(prodList[ii],Expression):
#                 expr.terms.extend(prodList[ii].terms)
#             elif isinstance(prodList[ii],Term):
#                 expr.terms.extend([prodList[ii]])
#         expr.simplify()
#         Lagrangian_us_gradient[i] += expr

    # Main Iterative Loop
    x_k_dict = dict(zip(names,x_k))
    objectiveValue = formulation.objective.substitute(x_k_dict)
    res.objectiveFunctionValues.append(objectiveValue)
    watchdogCounter = 0
    stepNumberIndex = -1
    constraintFunctionEvaluations = None
    for iterate in range(0,formulation.solverOptions.maxIterations):
        stepNumberIndex += 1
        res.stepNumber.append(stepNumberIndex)

        x_k_dict = dict(zip(names,x_k))
        subbedOG = substituteMatrix(objectiveGradient,x_k_dict,False)
    
        # =====================================================
        # =====================================================
        # =====================================================
        
        gpVarNames = names + [sv.name for sv in slackVariableList]
        Fmat = [[] for i in range(0,len(gpVarNames))]
        gvec = []
        K = []
        
        relaxedObjective = copy.deepcopy(formulation.objective)
        for vr in slackVariableList:
            relaxedObjective *= vr**formulation.solverOptions.penaltyExponent

        Nterms = 0
        for tm in relaxedObjective.terms:
            Nterms += 1
            gvec.append(tm.constant.magnitude)
            exps = [0.0 for i in range(0,len(gpVarNames))]
            for i in range(0,len(tm.variables)):
                vr = tm.variables[i]
                exps[gpVarNames.index(vr.name)] = tm.exponents[i]
            for i in range(0,len(gpVarNames)):
                Fmat[i].append(exps[i]) 
        K.append(Nterms)

        if constraintFunctionEvaluations is None: # should only be first iteration
            constraintFunctionEvaluations = [[],[],[]]
            for con in gpcConstraints:
                vl = con.LHS.substitute(x_k_dict)
                constraintFunctionEvaluations[0].append(vl)
            for fr in ngpcConstraintsLHS:
                vl = fr.substitute(x_k_dict)
                constraintFunctionEvaluations[1].append(vl)
            for con in runtimeConstraints:
                rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
                RC_xk = [x_k_dict[rcnm] for rcnm in rc_varnames]
                vl = con.returnConstraintSet(RC_xk, 'sgp') 
                constraintFunctionEvaluations[2].append(vl)
                
        constraintSlackIndex = -1
        
        for con in gpcConstraints:
            Nterms = 0
            constraintSlackIndex += 1
            for tm in con.LHS.terms:
                Nterms += 1
                gvec.append(tm.constant.magnitude)
                exps = [0.0 for i in range(0,len(gpVarNames))]
                for i in range(0,len(tm.variables)):
                    vr = tm.variables[i]
                    exps[gpVarNames.index(vr.name)] = tm.exponents[i]
                exps[gpVarNames.index(slackVariableList[constraintSlackIndex].name)] = -1.0
                for i in range(0,len(gpVarNames)):
                    Fmat[i].append(exps[i]) 
            K.append(Nterms)
            
        for con in ngpcConstraints:
            approxCon = con.approximateConvexifiedSignomial(x_k_dict)
            approxCon = approxCon[0]
            ncApprox = Constraint(approxCon.LHS/approxCon.RHS, '<=', 1.0*units.dimensionless)
            Nterms = 0
            constraintSlackIndex += 1
            for tm in ncApprox.LHS.terms:
                Nterms += 1
                gvec.append(tm.constant.magnitude)
                exps = [0.0 for i in range(0,len(gpVarNames))]
                for i in range(0,len(tm.variables)):
                    vr = tm.variables[i]
                    exps[gpVarNames.index(vr.name)] = tm.exponents[i]
                exps[gpVarNames.index(slackVariableList[constraintSlackIndex].name)] = -1.0
                for i in range(0,len(gpVarNames)):
                    Fmat[i].append(exps[i]) 
            K.append(Nterms)  
            
        for ii in range(0,len(runtimeConstraints)):
            con = runtimeConstraints[ii]  
            for j in range(0,len(con.outputs)):
                # Do the taylor expansion
                op = con.outputs[j]
                hes = runtimeConstraintHessians[ii][j]
                jac = constraintFunctionEvaluations[2][ii]['raw'][1][j]
                fvl = constraintFunctionEvaluations[2][ii]['raw'][0][j]
                teRHS = fvl
                for jj in range(0,len(con.inputs)):
                    teRHS += jac[jj]*(con.inputs[jj]-x_k_dict[con.inputs[jj].name])
                for jj in range(0,len(con.inputs)):
                    for jjj in range(0,len(con.inputs)):
                        teRHS += 1/2*(con.inputs[jj]-x_k_dict[con.inputs[jj].name])*hes[jj][jjj]*(con.inputs[jjj]-x_k_dict[con.inputs[jjj].name])
                teCon = Constraint(op, con.operators[j], teRHS)
                
                # Now treat like any other constraint
                approxConSet = teCon.approximateConvexifiedSignomial(x_k_dict)
                if teCon.operator == '==':
                    approxConSet = [Constraint(approxCon.LHS, '<=', approxCon.RHS),
                                    Constraint(approxCon.RHS, '<=', approxCon.LHS)]
                for approxCon in approxConSet:
                    ncApprox = Constraint(approxCon.LHS/approxCon.RHS, approxCon.operator, 1.0*units.dimensionless)
                    Nterms = 0
                    constraintSlackIndex += 1
                    for tm in ncApprox.LHS.terms:
                        Nterms += 1
                        gvec.append(tm.constant.magnitude)
                        exps = [0.0 for i in range(0,len(gpVarNames))]
                        for i in range(0,len(tm.variables)):
                            vr = tm.variables[i]
                            exps[gpVarNames.index(vr.name)] = tm.exponents[i]
                        exps[gpVarNames.index(slackVariableList[constraintSlackIndex].name)] = -1.0
                        for i in range(0,len(gpVarNames)):
                            Fmat[i].append(exps[i]) 
                    K.append(Nterms)  

        for sv in slackVariableList:
            nm = sv.name
            gvec.append(1.0)
            exps = [0.0 for i in range(0,len(gpVarNames))]
            exps[gpVarNames.index(nm)] = -1.0 
            for i in range(0,len(gpVarNames)):
                Fmat[i].append(exps[i]) 
            K.append(1)

        F = cvxopt.matrix(Fmat)
        g = cvxopt.log(cvxopt.matrix(gvec))
        sol = cvxopt.solvers.gp(K, F, g)

        gpConsLen = len(constraint_operators) + len(slackVariableList) 
        if len(sol['znl']) == gpConsLen: # only constraints included
            dv = sol['znl'][0:len(constraint_operators)]
        elif len(sol['znl']) == gpConsLen+1: # constraints and objective included, shouldn't happen but sometimes does
            dv = sol['znl'][1:len(constraint_operators)+1]
        else:
            raise ValueError('Returned an unexpected number of dual variables')
    
        
        GP_dualVariablesRaw = copy.deepcopy(dv)
        
        GP_vars = sol['x'][0:Nvars]
        GP_slacks = abs(sol['x'][Nvars:])
        GP_stepVars = [np.exp(GP_vars[i])*x_k[i].units-x_k[i] for i in range(0,Nvars)]

        # =====================================================
        # =====================================================
        # =====================================================
        for watchdog_Iterate in [0,1]:
            newObj = sol['primal objective'] * formulation.objective.units + formulation.objective.substitute(x_k_dict)

            slacks = GP_slacks
            fullStep_x = GP_stepVars #[rs.variables[nm] for nm in stepNames]
            
            if stepNumberIndex in range(0,len(formulation.solverOptions.stepSchedule)):
                stpLm = formulation.solverOptions.stepSchedule[stepNumberIndex]
                maxAllowedStep = min(stpLm * abs(np.array([(x_k[idx]/fullStep_x[idx]).to('').magnitude for idx in range(0,Nvars)])))
            else:
                maxAllowedStep = 1.0

            if iterate == 0:
                mu_k = [1.0 * formulation.objective.units/constraint_units[i] for i in range(0,Ncons)]
                
            eta = formulation.solverOptions.eta
            zeta = formulation.solverOptions.zeta
            tau = formulation.solverOptions.tau
            stepSize = 1.0

            cfe_saveForCache = copy.deepcopy(constraintFunctionEvaluations)
            
            MF_output = meritFunction(  vrs,
                                        formulation.solverOptions.rho,
                                        formulation.solverOptions.mu_margin,
                                        formulation.objective,
                                        [gpcConstraints, ngpcConstraintsLHS, runtimeConstraints], 
                                        0.0, 
                                        fullStep_x, 
                                        x_k_dict, 
                                        subbedOG, 
                                        mu_k, 
                                        GP_dualVariablesRaw, 
                                        constraint_operators, 
                                        constraintFunctionEvaluations, 
                                        mode = 1)
            
            phi_k = MF_output[0]
            D     = MF_output[1]
            mu_k  = MF_output[2]
            constraintFunctionEvaluations_previousIteration = MF_output[3]
            newLagrangeMultipliers = MF_output[4]          
            for i in range(0,formulation.solverOptions.maxStepSizeTries):
                doWolfeChecking = True
                
                MF_output = meritFunction(  vrs,
                                            formulation.solverOptions.rho,
                                            formulation.solverOptions.mu_margin,
                                            formulation.objective,
                                            [gpcConstraints, ngpcConstraintsLHS, runtimeConstraints], 
                                            stepSize, 
                                            fullStep_x, 
                                            x_k_dict, 
                                            subbedOG, 
                                            mu_k, 
                                            GP_dualVariablesRaw, 
                                            constraint_operators, 
                                            None, 
                                            mode = 2)

                phi_kp1 = MF_output[0]
                D_kp1   = MF_output[1]
                constraintFunctionEvaluations = MF_output[2]
                newLagrangeMultipliers = MF_output[3]
                
                try:
                    MF_output = meritFunction(  vrs,
                                                formulation.solverOptions.rho,
                                                formulation.solverOptions.mu_margin,
                                                formulation.objective,
                                                [gpcConstraints, ngpcConstraintsLHS, runtimeConstraints], 
                                                stepSize, 
                                                fullStep_x, 
                                                x_k_dict, 
                                                subbedOG, 
                                                mu_k, 
                                                GP_dualVariablesRaw, 
                                                constraint_operators, 
                                                None, 
                                                mode = 2)

                    phi_kp1 = MF_output[0]
                    D_kp1   = MF_output[1]
                    constraintFunctionEvaluations = MF_output[2]
                    newLagrangeMultipliers = MF_output[3]
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
                            x_k_cache = copy.deepcopy(x_k)
                            lam_cache = copy.deepcopy(lam)
                            GP_dualVariablesRaw_cache = copy.deepcopy(GP_dualVariablesRaw)
                            GP_slacks_cache = copy.deepcopy(GP_slacks)
                            GP_stepVars_cache = copy.deepcopy(GP_stepVars)
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
                                x_k = x_k_cache
                                lam = lam_cache
                                GP_dualVariablesRaw = GP_dualVariablesRaw_cache
                                GP_slacks = GP_slacks_cache
                                GP_stepVars = GP_stepVars_cache
                                constraintFunctionEvaluations = constraintFunctionEvaluations_cache
                                x_k_dict = dict(zip(names,x_k))
                                subbedOG = substituteMatrix(objectiveGradient,x_k_dict)
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
            constraintFunctionEvaluations = [[],[],[]]
            for con in gpcConstraints:
                vl = con.LHS.substitute(tmp_x_k_dict)
                constraintFunctionEvaluations[0].append(vl)
            for fr in ngpcConstraintsLHS:
                vl = fr.substitute(tmp_x_k_dict)
                constraintFunctionEvaluations[1].append(vl)
            for con in runtimeConstraints:
                rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
                RC_xk = [tmp_x_k_dict[rcnm] for rcnm in rc_varnames]
                vl = con.returnConstraintSet(RC_xk, 'sgp') 
                constraintFunctionEvaluations[2].append(vl)
                
        if watchdogCounter > formulation.solverOptions.watchdogIterations:
            watchdogCounter = 0

        res.intermediateSolves.append(sol)
        res.steps.append(stepSize)
        
        x_kp1   = [x_k[i]  + stepSize * fullStep_x[i] for i in range(0,Nvars)]
        fullStep_l = [newLagrangeMultipliers[i].magnitude - lam[i] for i in range(0,Ncons)]  
        lam_kp1 = [lam[i]  + stepSize * fullStep_l[i] for i in range(0,Ncons)]
        lamSubs = [lam_kp1[i] * lamv[i].units for i in range(0,Ncons)]
        res.designPoints.append(x_kp1)
        res.lagrangeMultipliers.append(lamSubs)
        x_kp1_dict = dict(zip(names,x_kp1))
        newObj = formulation.objective.substitute(x_kp1_dict)
        res.objectiveFunctionValues.append(newObj)
        
        lamDict = dict(zip(lamNames,copy.deepcopy(lamSubs)))
        s_k    = [stepSize * fullStep_x[i] for i in range(0,Nvars)]
        gl_k = []
        gl_kp1 = []

        
        rtc_LGk = [0.0*formulation.objective.units / v.units for v in vrs]
        rtc_LGkp1 = [0.0*formulation.objective.units / v.units for v in vrs]
        rtc_vrs = [ [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs]   for rc in runtimeConstraints]
        for i in range(0,Nvars):
            nm = names[i]
            constraintIndex = Ncons_explicit - 1
            for ii in range(0,len(runtimeConstraints)):
                rc = runtimeConstraints[ii]
                for jj in range(0,len(rc.outputs)):
                    op = rc.outputs[jj]
                    constraintIndex += 1
                    rcNames = [ip.name for ip in rc.inputs] + [op.name] 
                    current_lam = lamSubs[constraintIndex]
                    if nm in rcNames:
                        dc_dxkp1 = constraintFunctionEvaluations[2][ii]['con'][1][jj][rcNames.index(nm)]
                        rtc_LGkp1[i] += current_lam * dc_dxkp1
        
        combinedDict = copy.deepcopy(x_k_dict)
        
        # In regular space
        # for ky, vl in lamDict.items():
        #     combinedDict[ky] = vl
        # for i in range(0,len(vrs)):
        #     if hasattr(objectiveGradient[i],'substitute'):
        #         lg = objectiveGradient[i].substitute(combinedDict)
        #     else: # is a quantity
        #         lg = objectiveGradient[i]
        #     for ii in range(0,len(constraintGradients)):
        #         cg = constraintGradients[ii]
        #         if hasattr(cg[i],'substitute'):
        #             lg += lamSubs[ii] * cg[i].substitute(combinedDict)
        #         else: # is quantity
        #             lg += lamSubs[ii] * cg[i]
        #     lg += rtc_LGkp1[i]
        #     gl_kp1.append(lg)   

        # In logspace
        for ky, vl in lamDict.items():
            combinedDict[ky] = vl
        for i in range(0,len(vrs)):
            if hasattr(objectiveGradient[i],'substitute'):
                ogv = objectiveGradient[i].substitute(combinedDict)
            else: # is a quantity
                ogv = objectiveGradient[i]
            lg = 1/newObj * ogv * x_kp1[i]
            for ii in range(0,len(constraintGradients)):
                if ii < len(gpcConstraints):
                    conval = constraintFunctionEvaluations[0][ii]
                else:
                    conval = constraintFunctionEvaluations[1][ii-len(gpcConstraints)]
                cg = constraintGradients[ii]
                if hasattr(cg[i],'substitute'):
                    lg += GP_dualVariablesRaw[ii] * cg[i].substitute(combinedDict) * x_kp1[i] / conval
                else: # is quantity
                    lg += GP_dualVariablesRaw[ii] * cg[i] * x_kp1[i] / conval
            # lg += rtc_LGkp1[i]
            gl_kp1.append(lg)  


        magStep = sum([svl.magnitude**2 for svl in s_k])**0.5 
        gradLag = stripUnitsMat(gl_kp1)
        res.magnitudeOfStep.append(magStep)
        res.gradientOfLagrangian.append(gradLag)
        tc1 = all(abs(gradLag) < formulation.solverOptions.lagrangianGradientTolerance)
        tc2 = res.magnitudeOfStep[-1] <= formulation.solverOptions.stepMagnitudeTolerance

        if tc1:
            res.terminationStatus = ('Converged, gradient of Lagrangian (%e) within tolerance (%e) '%( max(abs(gradLag)),
                                               formulation.solverOptions.lagrangianGradientTolerance) )
            res.terminationFlag = True
            break
        elif tc2:
            res.terminationStatus = ('Converged, magnitude of step (%.2e)'%(res.magnitudeOfStep[-1]) +
                                     ' is witin specified tolerance (%.2e)'%(formulation.solverOptions.stepMagnitudeTolerance))
            res.terminationFlag = True
            break
        else:
            #Update Hessians
            for ii in range(0,len(runtimeConstraints)):
                con = runtimeConstraints[ii]
                for j in range(0,len(con.outputs)):
                    Bmat = runtimeConstraintHessians[ii][j]
                    g_k = constraintFunctionEvaluations_previousIteration[2][ii]['raw'][1][j]
                    g_kp1 = constraintFunctionEvaluations[2][ii]['raw'][1][j]
            
                    names_local = [ ip.name for ip in con.inputs]
                    step_indicies = [names.index(nm) for nm in names_local]
                    s_local = [s_k[i] for i in step_indicies] # need to extract only relevant step values
                    
                    # ===================================
                    # BFGS
                    # ===================================
                    y_k = [g_kp1[i] - g_k[i] for i in range(0,len(g_k))]
                    sTy = sum([s_local[i] * y_k[i] for i in range(0,len(s_local))])
                    Bs = [sum([Bmat[i][j] * s_local[j] for j in range(0,len(s_local))]) for i in range(0,len(Bmat))]
                    sTBs = sum([s_local[i] * Bs[i] for i in range(0,len(s_local))])
                    #------------------------------------
                    # undamped
                    #------------------------------------
#                     sTy = sum([s_local[i] * y_k[i] for i in range(0,len(s_local))])
#                     yyT_sTy = [[y_k[i] * y_k[j]/sTy for i in range(0,len(y_k))] for j in range(0,len(y_k))]
#                     sTB = [sum([s_local[i] * Bmat[i][j] for i in range(0,len(s_local))]) for j in range(0,len(Bmat))]
#                     BssTB_sTBs = [[Bs[i] * sTB[j]/sTBs for i in range(0,len(Bs))] for j in range(0,len(sTB))]
#                     Bmat_new = [[Bmat[i][j] - BssTB_sTBs[i][j] + yyT_sTy[i][j] for i in range(0,len(Bmat)) ] for j in range(0,len(Bmat[0]))]
                    #------------------------------------
                    # damped
                    #------------------------------------
#                     thetaCon = sTy >= 0.2*sTBs
#                     if thetaCon:
#                         theta = 1.0 * units.dimensionless
#                     else:
#                         theta = 0.8*sTBs / (sTBs - sTy)
#                     r_k = [theta * y_k[i] + (1-theta)*Bs[i] for i in range(0,len(y_k))]
#                     sTr = sum([s_local[i] * r_k[i] for i in range(0,len(s_local))])
#                     rrT_sTr = [[r_k[i] * r_k[j]/sTr for i in range(0,len(r_k))] for j in range(0,len(r_k))]
#                     sTB = [sum([s_local[i] * Bmat[i][j] for i in range(0,len(s_local))]) for j in range(0,len(Bmat))]
#                     BssTB_sTBs = [[Bs[i] * sTB[j]/sTBs for i in range(0,len(Bs))] for j in range(0,len(sTB))]
#                     Bmat_new = [[Bmat[i][j] - BssTB_sTBs[i][j] + rrT_sTr[i][j] for i in range(0,len(Bmat)) ] for j in range(0,len(Bmat[0]))]
                    #------------------------------------
                    # ===================================
                    # SR1
                    # ===================================
                    r = 1e-8 # from nocedal and wright
                    y_k = [g_kp1[i] - g_k[i] for i in range(0,len(g_k))]
                    Bs = [sum([Bmat[i][j] * s_local[j] for j in range(0,len(s_local))]) for i in range(0,len(Bmat))]
                    ymBs = [y_k[i] - Bs[i] for i in range(0,len(y_k))]
                    nmr = [[ymBs[i]*ymBs[j] for i in range(0,len(ymBs))] for j in range(0,len(ymBs))]
                    ymBsT_s = sum([ymBs[i] * s_local[i] for i in range(0,len(s_local))])
                    
                    normS = sum([sv**2 for sv in s_local])**(0.5) 
                    normO = sum([vl**2 for vl in ymBs])**(0.5)
                    
                    if abs(ymBsT_s) >= r*normS*normO:
                        upm = [[ymBs[i]*ymBs[j]/ymBsT_s for i in range(0,len(ymBs))] for j in range(0,len(ymBs))]
                        Bmat_new = [[Bmat[i][j] + upm[i][j] for i in range(0,len(Bmat)) ] for j in range(0,len(Bmat[0]))]
                    else:
                        # else no update required, hessian is already a good approximation
                        Bmat_new = Bmat
                    # ===================================
                    runtimeConstraintHessians[ii][j] = Bmat_new
                    # ===================================
                    
            # update x_k, lam_k
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
    # res.sensitivities = res.intermediateSolves[-1].sensitivities

    for i in range(0,len(slacks)):
        if abs(slacks[i]) > formulation.solverOptions.constraintTolerance:
            res.messages.append('WARNING: One or more constraint tolerances not met at the final solution (tol=%.2e)'%(formulation.solverOptions.constraintTolerance))

            
    return res