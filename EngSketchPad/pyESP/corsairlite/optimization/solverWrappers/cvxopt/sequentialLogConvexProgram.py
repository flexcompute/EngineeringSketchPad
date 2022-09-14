import copy
import math
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
        # gpc
        for con in constraints[0]:
            constraintIndex += 1
            vl = con.LHS.substitute(x_k_dict)
            constraintFunctionEvaluations[0].append(vl)
            constraintValues.append(vl)
            c_1Norm.append(returncval(vl,constraint_operators[constraintIndex]))
        # ngpc
        for fr in constraints[1]:
            constraintIndex += 1
            vl = fr.substitute(x_k_dict)
            constraintFunctionEvaluations[1].append(vl)
            constraintValues.append(vl)
            c_1Norm.append(returncval(vl,constraint_operators[constraintIndex]))
        # runtime
        for con in constraints[2]:
            rc_varnames = [ip.name for ip in con.inputs] + [op.name for op in con.outputs] 
            RC_xk = [x_k_dict[rcnm] for rcnm in rc_varnames]
            rs = con.returnConstraintSet(RC_xk, 'sgp')
            testValues = [(vl*units.dimensionless).magnitude for vl in rs['raw'][0]]
            if any(np.isnan(testValues)):
                raise RuntimeError('A RuntimeConstraint threw a Nan value')
            if any([tv <=0 for tv in testValues]):
                raise RuntimeError('A RuntimeConstraint threw a negative value')
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

def SequentialLogConvexProgramSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7
    
    # Setup output
    res = OptimizationOutput()
    # res.initialFormulation = copy.deepcopy(formulation)
    res.initialFormulation = formulation
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
    res.magnitudeOfStepScaled = []
    res.gradientOfLagrangian = []
    res.messages = []
    res.terminationStatus = 'Did not converge, maximum number of iterations reached (%d)'%(formulation.solverOptions.maxIterations)
    res.terminationFlag = False
    # debuggingCounter = formulation.solverOptions.counter

    # Setup key values used throughout run
    gpcConstraints = []
    ngpcConstraints = []
    runtimeConstraints = []
    for con in formulation.constraints:
        if isinstance(con,Constraint):
            try: # is gp compatible
                gpConstraintSet = con.toPosynomial(GPCompatible = True, leaveMonomialEqualities=True)
                for i in range(0,len(gpConstraintSet)):
                    gpConstraintSet[i].to_base_units()
                gpcConstraints.extend(gpConstraintSet)
            except: # is explicit but not gp compatible
                zeroedCon = con.toZeroed()
                zeroedCon.to_base_units()
                # stdcon = con.toLSQP_StandardForm()
                # stdcon.to_base_units()
                # stdcon = stdcon.substitute()
                # ngpcConstraints.append(stdcon.substitute())
                if zeroedCon.operator == '>=':
                    zeroedCon *= -1
                    ngpcConstraints.append(zeroedCon)
                else:
                    ngpcConstraints.append(zeroedCon)
        else: #runtime constraint
            # runtimeConstraints.append(copy.deepcopy(con))
            runtimeConstraints.append(con)
    
    # This is temporary until the fraction class is fully integrated
    for i in range(0,len(gpcConstraints)):
        gpcConstraints[i] = gpcConstraints[i].substitute()
    for i in range(0,len(ngpcConstraints)):
        ngpcConstraints[i] = ngpcConstraints[i].substitute()
    
    allConstraints = gpcConstraints + ngpcConstraints + runtimeConstraints
    formulation.skipEnforceBounds = True
    # formulation.constraints = copy.deepcopy(allConstraints)
    formulation.constraints = allConstraints
    formulation.skipEnforceBounds = False
    ngpcConstraintsLHS = []
    for i in range(0,len(ngpcConstraints)):
        con = ngpcConstraints[i]
        posCon = con.toAllPositive()
        if posCon.operator == '>=':
            posCon.swap()
        frac = posCon.LHS/posCon.RHS
        ngpcConstraintsLHS.append(frac)
    # for i in range(0,len(ngpcConstraints)):
    #     ngpcConstraintsLHS.append(ngpcConstraints[i].LHS)
    
    allConstraints = gpcConstraints + ngpcConstraints + runtimeConstraints
    explicitConstraints = gpcConstraints + ngpcConstraints

    rtcVariables = [ [ip.name for ip in rc.inputs] + [op.name for op in rc.outputs]   for rc in runtimeConstraints]
    formulation.solverOptions = copy.deepcopy(formulation.solverOptions)
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
            constraint_units.append(units.dimensionless)
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
                     description = "Slack number %d for SLCP"%(slackCounter))
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
                         description = "Dual Variable %d for SLCP"%(lamCounter))
            lamCounter += 1
            lamv.append(lm)

    # Initalize the lagrange mutipliers
    lam = [1.0 for i in range(0,Ncons)]

    # Set up hesian approximation (initalizes to identity)
    Bmat = np.eye(Nvars).tolist()
    res.hessianApproximations.append(Bmat)
        
    # Determine Gradients
    objectiveGradient = [formulation.objective.derivative(v) for v in vrs]
    #This is temporary until the fraction class gets fully integrated
    constraintGradients = [[gpcConstraints[i].LHS.derivative(v) for v in vrs] for i in range(0,len(gpcConstraints))] # only the gp compatible constraints
    constraintGradients_ngp = [[ngpcConstraintsLHS[i].derivative(v) for v in vrs] for i in range(0,len(ngpcConstraintsLHS))] # does the rest of the explicit constraints
    constraintGradients.extend(constraintGradients_ngp)

    # Main Iterative Loop
    x_k_dict = dict(zip(names,x_k))
    objectiveValue = formulation.objective.substitute(x_k_dict)
    res.objectiveFunctionValues.append(objectiveValue)
    watchdogCounter = 0
    stepNumberIndex = -1
    constraintFunctionEvaluations = None
    Astatic = None
    bstatic = None
    Gstatic = None
    hstatic = None
    posyCons = None
    posyInd = None
    for iterate in range(0,formulation.solverOptions.maxIterations):
        stepNumberIndex += 1
        res.stepNumber.append(stepNumberIndex)

        x_k_dict = dict(zip(names,x_k))
        subbedOG = substituteMatrix(objectiveGradient,x_k_dict,False)
    
        # =====================================================
        # =====================================================
        # =====================================================
        
        gpVarNames = names + [sv.name for sv in slackVariableList]
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

        A = []
        b = []
        G = []
        h = []
        posyCons = []
        moniInd = []
        moneInd = []
        posyInd = []
        conIndex = -1
        for i in range(0,len(gpcConstraints)):
            conIndex += 1
            con = copy.deepcopy(gpcConstraints[i])
            if con.operator == '==':
                #monomial equality, goes to A,b
                const = con.LHS.terms[0].constant * units.dimensionless
                termVars = con.LHS.terms[0].variables
                termNames = [v.name for v in termVars]
                expnts = con.LHS.terms[0].exponents
                vridxs = [names.index(nm) for nm in termNames]
                Arow = [0.0]*Nvars
                for ii in range(0,len(expnts)):
                    Arow[vridxs[ii]] = expnts[ii]
                brow = -np.log(const.magnitude) - sum([Arow[ii] * np.log(x_k[ii].magnitude) for ii in range(0,len(x_k))])
                A.append(Arow)
                b.append(brow)
                moneInd.append(conIndex)
            else:
                if len(con.LHS.terms) == 1:
                    #monomial
                    const = con.LHS.terms[0].constant * units.dimensionless
                    termVars = con.LHS.terms[0].variables
                    termNames = [v.name for v in termVars]
                    expnts = con.LHS.terms[0].exponents
                    vridxs = [names.index(nm) for nm in termNames]
                    Grow = [0.0]*Nvars
                    for ii in range(0,len(expnts)):
                        Grow[vridxs[ii]] = expnts[ii]
                    hrow = -np.log(const.magnitude) - sum([Grow[ii] * np.log(x_k[ii].magnitude) for ii in range(0,len(x_k))])
                    G.append(Grow)
                    h.append(hrow)
                    moniInd.append(conIndex)
                else:
                    #posynomial: save for function F
                    posyCons.append(con)
                    posyInd.append(conIndex)
                        
                        
        # Generate Function F
        Fmat = [[] for i in range(0,len(gpVarNames))]
        gvec = []
        K = []
        Nterms = 0

        Nterms += 1
        gvec.append(1)
        exps = [0.0 for i in range(0,len(gpVarNames))]
        for ii in range(0,len(gpVarNames)):
            Fmat[ii].append(exps[ii]) 
        K.append(Nterms)
        
        for ii in range(0,len(posyCons)):
            con = posyCons[ii]
            Nterms = 0
            lhs = copy.deepcopy(con.LHS)/slackVariableList[Ncons-len(posyCons)+ii]
            for tm in lhs.terms:
                Nterms += 1
                gvec.append(tm.constant.magnitude)
                exps = [0.0 for iii in range(0,len(gpVarNames))]
                for j in range(0,len(tm.variables)):
                    vr = tm.variables[j]
                    exps[gpVarNames.index(vr.name)] = tm.exponents[j]
                for jj in range(0,len(gpVarNames)):
                    Fmat[jj].append(exps[jj]) 
            K.append(Nterms)

        F = cvxopt.matrix(Fmat)
        g = cvxopt.log(cvxopt.matrix(gvec))
        
        if type(K) is not list or [ k for k in K if type(k) is not int 
            or k <= 0 ]:
            raise TypeError("'K' must be a list of positive integers")
        mnl = len(K)-1
        l = sum(K)
        if type(F) not in (cvxopt.base.matrix, cvxopt.base.spmatrix) or F.typecode != 'd' or \
            F.size[0] != l:
            raise TypeError("'F' must be a dense or sparse 'd' matrix "\
                "with %d rows" %l)
        if type(g) is not cvxopt.base.matrix or g.typecode != 'd' or g.size != (l,1): 
            raise TypeError("'g' must be a dene 'd' matrix of "\
                "size (%d,1)" %l)
        n = F.size[1]
        y = cvxopt.base.matrix(0.0, (l,1))
        Fsc = cvxopt.base.matrix(0.0, (max(K),n))
        cs1 = [ sum(K[:i]) for i in range(mnl+1) ] 
        cs2 = [ cs1[i] + K[i] for i in range(mnl+1) ]
        ind = list(zip(range(mnl+1), cs1, cs2))
            
        objVal = formulation.objective.substitute(x_k_dict)
        TS0 = np.log(objVal.magnitude)
        TS1 = [(subbedOG[i] * x_k[i] / objVal).magnitude for i in range(0,Nvars)]
        TS2 = stripUnitsMat(Bmat).tolist()

        def Fnonlin(x = None, z = None):
            # ===================================================================
            # This stuff is the GP solver from cvxopt, used to do the constraints
            # ===================================================================
            if x is None: 
                return mnl, cvxopt.base.matrix(0.0, (n,1))

            # need to correct so that posynomials use variables d, not x
            current_xk = cvxopt.matrix([np.log(vl.magnitude) for vl in x_k] + [0.0]*(len(x)-Nvars))
            f = cvxopt.base.matrix(0.0, (mnl+1,1))
            Df = cvxopt.base.matrix(0.0, (mnl+1,n))
            gxfm = F*current_xk+g
            cvxopt.blas.copy(gxfm, y)
            cvxopt.base.gemv(F, x, y, beta=1.0)

            # Need to store this to use with corrected objective hessian later
            if z is not None: 
                z_obj = z[0]
                z[0] = 0.0
                H = cvxopt.base.matrix(0.0, (n,n)) #this is the only line from cvxopt in this if statement

            for i, start, stop in ind:
                ymax = max(y[start:stop])
                y[start:stop] = cvxopt.base.exp(y[start:stop] - ymax)
                ysum = cvxopt.blas.asum(y, n=stop-start, offset=start)
                f[i] = ymax + math.log(ysum)
                cvxopt.blas.scal(1.0/ysum, y, n=stop-start, offset=start)
                cvxopt.base.gemv(F, y, Df, trans='T', m=stop-start, incy=mnl+1, offsetA=start, offsetx=start, offsety=i)
                if z is not None:
                    Fsc[:K[i], :] = F[start:stop, :] 
                    for k in range(start,stop):
                        cvxopt.blas.axpy(Df, Fsc, n=n, alpha=-1.0, incx=mnl+1, incy=Fsc.size[0], offsetx=i, offsety=k-start)
                        cvxopt.blas.scal(math.sqrt(y[k]), Fsc, inc=Fsc.size[0], offset=k-start)
                    cvxopt.blas.syrk(Fsc, H, trans='T', k=stop-start, alpha=z[i], beta=1.0)

            # ===================================================================
            # Need to overwrite the objective data
            # ===================================================================
            Kpen = formulation.solverOptions.penaltyConstant
                
            x_numpy = np.array(x[0:Nvars]).flatten()
            slacks_numpy = np.array(x[Nvars:]).flatten()
        
            val_obj = TS0 + np.dot(x_numpy,np.array(TS1)) + 0.5 * np.dot(np.dot(x_numpy.T,TS2),x_numpy) + np.dot(Kpen*slacks_numpy.T,slacks_numpy)
            
            grad_obj = []
            for i in range(0,len(x_numpy)):
                grad_obj.append(TS1[i] + np.dot(TS2[i],x_numpy))
            for i in range(0,len(slacks_numpy)):
                grad_obj.append(2*Kpen*slacks_numpy[i])
                
            f[0] = val_obj
            Df = np.array(Df)
            Df[0] = np.array(grad_obj)
            Df = cvxopt.matrix(Df)

            if z is None: 
                return f, Df
            else: 
                obj_H = []
                for rw in TS2:
                    obj_H.append(rw + [0.0]*len(slacks_numpy))
                for i in range(0,len(slacks_numpy)):
                    obj_H.append([0.0]*Nvars + [0.0]*i + [2.0*Kpen] + [0.0]*(len(slacks_numpy)-i-1) )

                H += z_obj*cvxopt.matrix(obj_H)
                return f, Df, H

        for i in range(0,len(ngpcConstraints)):
            conIndex += 1
            feval = constraintFunctionEvaluations[1][i] * units.dimensionless
            grad  = constraintGradients_ngp[i]
            log_feval = np.log(feval.magnitude)
            log_grad = []
            for ii in range(0,len(x_k)):
                if hasattr(grad[ii],'substitute'):
                    log_grad.append( (grad[ii].substitute(x_k_dict) * x_k[ii] / feval).magnitude )
                else:
                    log_grad.append( (grad[ii] * x_k[ii] / feval).magnitude )
            if ngpcConstraints[i].operator == '==':
                A.append(log_grad)
                b.append(-log_feval)
                moneInd.append(conIndex)
            else:
                G.append(log_grad)
                h.append(-log_feval)
                moniInd.append(conIndex)
            
        for i in range(0,len(runtimeConstraints)):
            con = runtimeConstraints[i]  
            for jjj in range(0,len(con.outputs)):
                conIndex += 1
                feval = constraintFunctionEvaluations[2][i]['con'][0][jjj] * units.dimensionless
                grad = constraintFunctionEvaluations[2][i]['con'][1][jjj]
                
                conVarnames = [cvv.name for cvv in con.inputs] + [con.outputs[jjj].name]
                grad_corrected = [0.0*feval.units/x_k[ii].units for ii in range(0,len(x_k))]
                for ii in range(0,len(grad)):
                    gv = grad[ii]
                    gn = conVarnames[ii]
                    grad_corrected[names.index(gn)] = gv
                log_feval = np.log(feval.magnitude) 
                log_grad = [ (grad_corrected[ii] * x_k[ii] / feval).magnitude for ii in range(0,len(x_k))]

                # cody
                
                if con.operators[jjj] == '==':
                    A.append(log_grad)
                    b.append(-log_feval)
                    moneInd.append(conIndex)
                else:
                    G.append(log_grad)
                    h.append(-log_feval)
                    moniInd.append(conIndex)
            
        numCons = len(b) + len(h) + len(posyCons)
        for i in range(0,len(b)):
            slackVec = [0.0]*numCons
            slackVec[i] = -1.0
            A[i].extend(slackVec)
        for i in range(len(b),len(b)+len(h)):
            slackVec = [0.0]*numCons
            slackVec[i] = -1.0
            G[i-len(b)].extend(slackVec)

        if len(b) == 0:
            A=None
            b=None
        else:
            A = cvxopt.matrix(A).T
            b = cvxopt.matrix(b)
    
        if len(h) == 0:
            G=None
            h=None
        else:
            G = cvxopt.matrix(G).T
            h = cvxopt.matrix(h)

        # G = cvxopt.matrix(G).T
        # h = cvxopt.matrix(h)
        # A = cvxopt.matrix(A).T
        # b = cvxopt.matrix(b)

        # writeData = {'Fpre':copy.deepcopy(Fnonlin),'Gpre':copy.deepcopy(G),'hpre':copy.deepcopy(h),'Apre':copy.deepcopy(A),'bpre':copy.deepcopy(b)}

        sol = cvxopt.solvers.cp(F=Fnonlin, G=G, h=h, A=A, b=b)
        # writeData['sol'] = copy.deepcopy(sol)
        # writeData['F'] = copy.deepcopy(Fnonlin)
        # writeData['G'] = copy.deepcopy(G)
        # writeData['h'] = copy.deepcopy(h)
        # writeData['A'] = copy.deepcopy(A)
        # writeData['b'] = copy.deepcopy(b)

        # setPickle()
        # with open('cvxDataDump_' + str(debuggingCounter) + '_' + str(iterate) + '.dl', "wb") as dill_file:
        #     dill.dump(writeData, dill_file)
        # unsetPickle()


        zlnL = list(sol['znl'])
        
        if len(zlnL) == len(posyCons):
            pass # valid
        elif len(zlnL) == len(posyCons)+1:
            zlnL = zlnL[1:]
        else:
            raise ValueError('Length of nonlinear sensitivities does not match expectation')
        
        dualStack = list(sol['zl']) + list(sol['y']) + zlnL
        indexStack = moniInd + moneInd + posyInd
        
        LCP_dualVariablesRaw = [dualStack[indexStack.index(vl)] for vl in range(0,len(dualStack))]
        LCP_slacks = abs(sol['x'][Nvars:])
        LCP_stepVars = [np.exp(sol['x'][i]+np.log(x_k[i].magnitude))*x_k[i].units - x_k[i] for i in range(0,Nvars)]
        
        # =====================================================
        # =====================================================
        # =====================================================
        for watchdog_Iterate in [0,1]:
            newObj = sol['primal objective'] * formulation.objective.units + formulation.objective.substitute(x_k_dict)
            
            slacks = LCP_slacks
            fullStep_x = LCP_stepVars #[rs.variables[nm] for nm in stepNames]
            
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

            if stepNumberIndex in range(0,len(formulation.solverOptions.baseStepSchedule)):
                stepSize = formulation.solverOptions.baseStepSchedule[stepNumberIndex]
            else:
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
                                        LCP_dualVariablesRaw, 
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
                
                # MF_output = meritFunction(  vrs,
                #                             formulation.solverOptions.rho,
                #                             formulation.solverOptions.mu_margin,
                #                             formulation.objective,
                #                             [gpcConstraints, ngpcConstraintsLHS, runtimeConstraints], 
                #                             stepSize, 
                #                             fullStep_x, 
                #                             x_k_dict, 
                #                             subbedOG, 
                #                             mu_k, 
                #                             LCP_dualVariablesRaw, 
                #                             constraint_operators, 
                #                             None, 
                #                             mode = 2)

                # phi_kp1 = MF_output[0]
                # D_kp1   = MF_output[1]
                # constraintFunctionEvaluations = MF_output[2]
                # newLagrangeMultipliers = MF_output[3]
                
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
                                                LCP_dualVariablesRaw, 
                                                constraint_operators, 
                                                None, 
                                                mode = 2)

                    phi_kp1 = MF_output[0]
                    D_kp1   = MF_output[1]
                    constraintFunctionEvaluations = MF_output[2]
                    newLagrangeMultipliers = MF_output[3]
                except Exception as e:
                    # Step was too big, one of the constraints became imaginary
                    # print('============\nA Step Failed\n=============')
                    # print(fullStep_x)
                    doWolfeChecking = False
                    res.messages.append('Attempted step of size %.2e at iteration %d is too large and is causing a constraint to return as imaginary.  Decreasing according to the tau schedule.'%(
                                        stepSize, stepNumberIndex))
                    stepSize *= tau
                    # print(stepSize)
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
                            LCP_dualVariablesRaw_cache = copy.deepcopy(LCP_dualVariablesRaw)
                            LCP_slacks_cache = copy.deepcopy(LCP_slacks)
                            LCP_stepVars_cache = copy.deepcopy(LCP_stepVars)
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
                                LCP_dualVariablesRaw = LCP_dualVariablesRaw_cache
                                LCP_slacks = LCP_slacks_cache
                                LCP_stepVars = LCP_stepVars_cache
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
        if iterate == 0:
            res.designPoints.append(x_k)   
        res.designPoints.append(x_kp1)
        res.lagrangeMultipliers.append(lamSubs)
        x_kp1_dict = dict(zip(names,x_kp1))
        prevObj =formulation.objective.substitute(x_k_dict) 
        newObj = formulation.objective.substitute(x_kp1_dict)
        res.objectiveFunctionValues.append(newObj)
        
        lamDict = dict(zip(lamNames,copy.deepcopy(lamSubs)))
        s_k    = [stepSize * fullStep_x[i] for i in range(0,Nvars)]
        gl_k = []
        gl_kp1 = []
        glb_k = []
        glb_kp1 = []

        rtc_LGk   = [0.0*units.dimensionless]*len(vrs)
        rtc_LGkp1 = [0.0*units.dimensionless]*len(vrs)
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
                    current_lam = LCP_dualVariablesRaw[constraintIndex]
                    if nm in rcNames:
                        conval_k = constraintFunctionEvaluations_previousIteration[2][ii]['con'][0][jj]
                        dc_dxk   = constraintFunctionEvaluations_previousIteration[2][ii]['con'][1][jj][rcNames.index(nm)]
                        rtc_LGk[i] += current_lam * dc_dxk * x_k[names.index(nm)] / conval_k
            
                        conval_kp1 = constraintFunctionEvaluations[2][ii]['con'][0][jj]
                        dc_dxkp1   = constraintFunctionEvaluations[2][ii]['con'][1][jj][rcNames.index(nm)]
                        rtc_LGkp1[i] += current_lam * dc_dxkp1 * x_kp1[names.index(nm)] / conval_kp1
                        
        combinedDict_k = copy.deepcopy(x_k_dict)
        combinedDict_kp1 = copy.deepcopy(x_kp1_dict)

        # In logspace
        for ky, vl in lamDict.items():
            combinedDict_k[ky] = vl
            combinedDict_kp1[ky] = vl
        for i in range(0,len(vrs)):
            if hasattr(objectiveGradient[i],'substitute'):
                ogv_k = objectiveGradient[i].substitute(combinedDict_k)
                ogv_kp1 = objectiveGradient[i].substitute(combinedDict_kp1)
            else: # is a quantity
                ogv_k = objectiveGradient[i]
                ogv_kp1 = objectiveGradient[i]
            
            lg_k   = 1/prevObj * ogv_k   * x_k[i]
            lg_kp1 = 1/newObj  * ogv_kp1 * x_kp1[i]
            
            lgb_k   = 1/prevObj * ogv_k   * x_k[i]
            lgb_kp1 = 1/newObj  * ogv_kp1 * x_kp1[i]

            for ii in range(0,len(constraintGradients)):
                cg = constraintGradients[ii]
                if ii < len(gpcConstraints):
                    conval_k   = constraintFunctionEvaluations_previousIteration[0][ii]
                    conval_kp1 = constraintFunctionEvaluations[0][ii]
                else:
                    conval_k   = constraintFunctionEvaluations_previousIteration[1][ii-len(gpcConstraints)]
                    conval_kp1 = constraintFunctionEvaluations[1][ii-len(gpcConstraints)]
                    if hasattr(cg[i],'substitute'):
                        lgb_k   += LCP_dualVariablesRaw[ii] * cg[i].substitute(combinedDict_k)   * x_k[i]   / conval_k
                        lgb_kp1 += LCP_dualVariablesRaw[ii] * cg[i].substitute(combinedDict_kp1) * x_kp1[i] / conval_kp1
                    else: # is quantity
                        lgb_k   += LCP_dualVariablesRaw[ii] * cg[i] * x_k[i]   / conval_k
                        lgb_kp1 += LCP_dualVariablesRaw[ii] * cg[i] * x_kp1[i] / conval_kp1
                
                if hasattr(cg[i],'substitute'):
                    lg_k   += LCP_dualVariablesRaw[ii] * cg[i].substitute(combinedDict_k)   * x_k[i]   / conval_k
                    lg_kp1 += LCP_dualVariablesRaw[ii] * cg[i].substitute(combinedDict_kp1) * x_kp1[i] / conval_kp1
                else: # is quantity
                    lg_k   += LCP_dualVariablesRaw[ii] * cg[i] * x_k  [i] / conval_k
                    lg_kp1 += LCP_dualVariablesRaw[ii] * cg[i] * x_kp1[i] / conval_kp1
                    
            lg_k   += rtc_LGk[i]
            lg_kp1 += rtc_LGkp1[i]
            
            lgb_k   += rtc_LGk[i]
            lgb_kp1 += rtc_LGkp1[i]
            
            gl_k.append(lg_k)
            gl_kp1.append(lg_kp1)  
            
            glb_k.append(lgb_k)
            glb_kp1.append(lgb_kp1)  

        magStep = sum([svl.magnitude**2 for svl in s_k])**0.5
        magStep_scaled = sum([(s_k[i]/x_kp1[i])**2 for i in range(0,len(s_k))])**0.5
        gradLag = stripUnitsMat(gl_kp1)
        res.magnitudeOfStep.append(magStep)
        res.magnitudeOfStepScaled.append(magStep)
        res.gradientOfLagrangian.append(gradLag)
        tc1 = all(abs(gradLag) < formulation.solverOptions.lagrangianGradientTolerance)
        tc2 = res.magnitudeOfStepScaled[-1] <= formulation.solverOptions.stepMagnitudeTolerance
        if formulation.solverOptions.relativeTolerance is not None and iterate>1:
            relchange = abs(res.objectiveFunctionValues[-1] - res.objectiveFunctionValues[-2])/res.objectiveFunctionValues[-1]
            tc3 = relchange < formulation.solverOptions.relativeTolerance
        else:
            tc3 = False


        if tc1:
            res.terminationStatus = ('Converged, gradient of Lagrangian (%e) within tolerance (%e) '%( max(abs(gradLag)),
                                               formulation.solverOptions.lagrangianGradientTolerance) )
            res.terminationFlag = True
            break
        elif tc2 and iterate>0:
            res.terminationStatus = ('Converged, magnitude of step (%.2e)'%(res.magnitudeOfStep[-1]) +
                                     ' is witin specified tolerance (%.2e)'%(formulation.solverOptions.stepMagnitudeTolerance))
            res.terminationFlag = True
            break
        elif tc3:
            res.terminationStatus = ('Converged, relative change in the objective funciton (%.2e)'%(relchange) +
                                     ' is witin specified tolerance (%.2e)'%(formulation.solverOptions.relativeTolerance))
            res.terminationFlag = True
            break
        else:
            #Update Hessian
            g_k = glb_k
            g_kp1 = glb_kp1
            xlog_k   = [np.log(vl.magnitude) for vl in x_k  ]
            xlog_kp1 = [np.log(vl.magnitude) for vl in x_kp1]
            s_kLog = [xlog_kp1[i] - xlog_k[i] for i in range(0,len(xlog_k))]

            y_k = [g_kp1[i] - g_k[i] for i in range(0,len(g_k))]
            sTy = sum([s_kLog[i] * y_k[i] for i in range(0,len(s_kLog))])
            Bs = [sum([Bmat[i][j] * s_kLog[j] for j in range(0,len(s_kLog))]) for i in range(0,len(Bmat))]
            sTBs = sum([s_kLog[i] * Bs[i] for i in range(0,len(s_kLog))])

            thetaCon = sTy >= 0.2*sTBs
            if thetaCon:
                theta = 1.0 * units.dimensionless
            else:
                theta = 0.8*sTBs / (sTBs - sTy)
            r_k = [theta * y_k[i] + (1-theta)*Bs[i] for i in range(0,len(y_k))]
            sTr = sum([s_kLog[i] * r_k[i] for i in range(0,len(s_kLog))])
            rrT_sTr = [[r_k[i] * r_k[j]/sTr for i in range(0,len(r_k))] for j in range(0,len(r_k))]
            sTB = [sum([s_kLog[i] * Bmat[i][j] for i in range(0,len(s_kLog))]) for j in range(0,len(Bmat))]
            BssTB_sTBs = [[Bs[i] * sTB[j]/sTBs for i in range(0,len(Bs))] for j in range(0,len(sTB))]
            Bmat_new = [[Bmat[i][j] - BssTB_sTBs[i][j] + rrT_sTr[i][j] for i in range(0,len(Bmat)) ] for j in range(0,len(Bmat[0]))]

            Bmat = Bmat_new
            res.hessianApproximations.append(Bmat)

            # print(x_k)
            x_k = x_kp1
            # print(x_kp1)
            lam = lam_kp1

            if formulation.solverOptions.progressFilename is not None:
                setPickle()
                with open(formulation.solverOptions.progressFilename + "_iteration-%d.dl"%(iterate+1), "wb") as dill_file:
                    dill.dump(res, dill_file)
                with open(formulation.solverOptions.progressFilename, "wb") as dill_file:
                    dill.dump(res, dill_file)
                unsetPickle()

    vDict = dict(zip(names, x_kp1))

    # res.solvedFormulation = copy.deepcopy(formulation)
    # res.solvedFormulation = formulation
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

    for i in range(0,len(slacks)):
        if abs(slacks[i]) > formulation.solverOptions.constraintTolerance:
            res.messages.append('WARNING: One or more constraint tolerances not met at the final solution (tol=%.2e)'%(formulation.solverOptions.constraintTolerance))

            
    return res