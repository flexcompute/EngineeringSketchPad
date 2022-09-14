import copy as copy
import corsairlite
from corsairlite.core.dataTypes import OptimizationOutput, Container
from corsairlite.optimization.formulation import Formulation
# from corsairlite.optimization.solverWrappers import cvxopt as solverTree
from corsairlite.core.pickleFunctions import setPickle, unsetPickle
import numpy as np
import cvxopt
import dill

import time 

def SignomialProgramStandardSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-5
    cvxopt.solvers.options['reltol'] = 1e-5
    cvxopt.solvers.options['feastol'] = 1e-5

    outputMessages = []
    if hasattr(formulation.solverOptions,'autoSelect'):
        if formulation.solverOptions.autoSelect == True:
            outputMessages.append('Auto-selected solve type [sp]')

    options = copy.deepcopy(formulation.solverOptions)
    initialFormulation = copy.deepcopy(formulation)
    formulation = formulation.toSignomialProgram(substituted=False)
    names = [vbl.name for vbl in formulation.variables_only]

    if formulation.solverOptions.x0 is not None:
        designPoint = formulation.solverOptions.x0
    else:
        designPoint = {}
        for vr in formulation.variables_only:
            designPoint[vr.name] = vr.guess

    res = OptimizationOutput()
    res.initialFormulation = copy.deepcopy(initialFormulation)
    res.solveType = formulation.solverOptions.solveType
    res.solver = formulation.solverOptions.solver.lower()
    res.solvedFormulation = copy.deepcopy(formulation)
    res.designPoints = [designPoint]
    res.intermediateFormulations = []
    res.intermediateSolves = []
    res.error = []

    err = 10.0 * formulation.solverOptions.tol
    previousObjective = float("inf") * initialFormulation.objective.units
    terminationStatus = 'Did not converge, maximum number of iterations reached (%d)'%(formulation.solverOptions.maxIter)
    newConstraints = []
    for ii in range(0,formulation.solverOptions.maxIter):
        numIter = ii + 1
        newFormulation = Formulation()
        newConstraints = copy.deepcopy(newConstraints)
        if ii == 0:
            signomialIndicies = []
            signomialIndiciesNC = []
            constraintTraceback = []
            cixNC = 0
            for cix in range(0,len(formulation.constraints)):
                con = formulation.constraints[cix]
                try:
                    conSet = con.toPosynomial(GPCompatible = True)
                except:
                    conSet = con.approximateConvexifiedSignomial(designPoint)
                    signomialIndicies.append(cix)
                for sc in conSet:
                    csl = sc.toPosynomial(GPCompatible = True)
                    for cn in csl:
                        if cix in signomialIndicies:
                            signomialIndiciesNC.append(cix)
                        else:
                            signomialIndiciesNC.append(None)
                        cn.to_base_units()
                        cn_subbed = cn.substitute()
                        newConstraints.append(cn_subbed)
                        constraintTraceback.append(cix)
    
            Fmat = [[] for i in range(0,len(names))]
            gvec = []
            K = []
    
            lookupTable = [[] for iii in range(0,len(formulation.constraints)+1)]
            rowCounter = -1
    
            Nterms = 0
            sbd_objective = formulation.objective.substitute(preserveType=True)
            for tm in sbd_objective.terms:
                rowCounter += 1
                lookupTable[0].append(rowCounter)
                Nterms += 1
                gvec.append(tm.constant.magnitude)
                exps = [0.0 for i in range(0,len(names))]
                for i in range(0,len(tm.variables)):
                    vr = tm.variables[i]
                    exps[names.index(vr.name)] = tm.exponents[i]
                for i in range(0,len(names)):
                    Fmat[i].append(exps[i]) 
            K.append(Nterms)

            for con_ix in range(0,len(newConstraints)):
                con = newConstraints[con_ix]
                responsibleConstraint = constraintTraceback[con_ix]
                Nterms = 0
                for tm in con.LHS.terms:
                    rowCounter += 1
                    lookupTable[responsibleConstraint+1].append(rowCounter)
                    Nterms += 1
                    gvec.append(tm.constant.magnitude)
                    exps = [0.0 for i in range(0,len(names))]
                    for i in range(0,len(tm.variables)):
                        vr = tm.variables[i]
                        exps[names.index(vr.name)] = tm.exponents[i]
                    for i in range(0,len(names)):
                        Fmat[i].append(exps[i]) 
                K.append(Nterms)
        
        else:
            for cix in range(0,len(formulation.constraints)):
                if cix in signomialIndicies:
                    dependentRows = lookupTable[cix+1]
                    con = formulation.constraints[cix]
                    newConstraintsIndex = np.where(np.array(signomialIndiciesNC)==cix)[0]
                    conSet = con.approximateConvexifiedSignomial(designPoint)
                    if len(conSet) > 1:
                        raise ValueError('Monomial Approximation of a signomial is returning via an unexpected method')
                    csl = conSet[0].toPosynomial(GPCompatible = True)
                    if len(csl) == 1: #signomial inequality
                        cn = csl[0]
                        cn.to_base_units()
                        cn_subbed = cn.substitute()
                        for iii in range(0,len(cn_subbed.LHS.terms)):
                            tm = cn_subbed.LHS.terms[iii]
                            gvec[dependentRows[iii]] = tm.constant.magnitude
                            for i in range(0,len(tm.variables)):
                                vr = tm.variables[i]
                                Fmat[names.index(vr.name)][dependentRows[iii]] = tm.exponents[i]
                    else: # signomial equality, there will be exactly two dependent rows
                        for iii in range(0,len(csl)):
                            cn = csl[iii]
                            cn.to_base_units()
                            cn_subbed = cn.substitute()
                            for tm in cn_subbed.LHS.terms:
                                gvec[dependentRows[iii]] = tm.constant.magnitude
                                for i in range(0,len(tm.variables)):
                                    vr = tm.variables[i]
                                    Fmat[names.index(vr.name)][dependentRows[iii]] = tm.exponents[i]
        
        F = cvxopt.matrix(Fmat)
        g = cvxopt.log(cvxopt.matrix(gvec))
        sol = cvxopt.solvers.gp(K, F, g) 
        
        res.intermediateSolves.append(sol)
        solvedObjective = cvxopt.exp( sol['primal objective']) * formulation.objective.units

        vrs = list(cvxopt.exp( sol['x'] ))
        fm_vrs = formulation.variables_only
        for i in range(0,len(vrs)):
            unt = fm_vrs[i].units
            vrs[i] = vrs[i] * unt

        solvedObjective = cvxopt.exp( sol['primal objective']) * formulation.objective.units
        designPoint = dict(zip(names, vrs))
        
        err = (solvedObjective - previousObjective)/solvedObjective
        res.error.append(err)
        previousObjective = solvedObjective
        res.designPoints.append(designPoint)

        if abs(err) <= formulation.solverOptions.tol:
            terminationStatus = 'Converged, relative error (%.2e) is less than specified tolerance (%.2e)'%(abs(err), formulation.solverOptions.tol)
            break
        else:
            if formulation.solverOptions.progressFilename is not None:
                setPickle()
                with open(formulation.solverOptions.progressFilename + "_iteration-%d.dl"%(ii+1), "wb") as dill_file:
                    dill.dump(res, dill_file)
                with open(formulation.solverOptions.progressFilename, "wb") as dill_file:
                    dill.dump(res, dill_file)
                unsetPickle()

        N = formulation.solverOptions.convergenceCheckIterations
        if ii > N:
            er0 = copy.deepcopy(res.error)
            er1 = er0[-N:]
            er = np.array([er1[i].magnitude for i in range(0,len(er1))])
            poly = np.polyfit(np.linspace(1,N,N),er,1)
            slope = poly[0]
            if slope > formulation.solverOptions.convergenceCheckSlope:
                terminationStatus = 'Did not converge to specified tolerance (%.2e), but achieved a tolerance of %.2e'%(formulation.solverOptions.tol, res.error[-1])
                break
    # -----------------------------------------
    res.sensitivities = None
    res.sensitivityInfo = Container()
    res.sensitivityInfo.signomialIndicies = signomialIndicies
    res.sensitivityInfo.signomialIndiciesNC = signomialIndiciesNC
    res.sensitivityInfo.designPoint = copy.deepcopy(designPoint)
    res.sensitivityInfo.lookupTable = lookupTable
    res.sensitivityInfo.sol = sol
    res.variables = designPoint
    res.objective = solvedObjective
    res.numberOfIterations = numIter
    res.terminationStatus = terminationStatus
    res.messages = outputMessages
    return res