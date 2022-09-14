import copy as copy
import cvxopt
import numpy as np
from corsairlite import units, Q_
from corsairlite.core.dataTypes import OptimizationOutput

# Required Options:
# options.solver
# options.solveType
# options.returnRaw

def GeometricProgramSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7
    cvxopt.solvers.options['abstol'] = 1e-4
    cvxopt.solvers.options['reltol'] = 1e-4
    cvxopt.solvers.options['feastol'] = 1e-4

    outputMessages = []
    if hasattr(formulation.solverOptions,'autoSelect'):
        if formulation.solverOptions.autoSelect == True:
            outputMessages.append('Auto-selected solve type [gp]')

    options = copy.deepcopy(formulation.solverOptions)
    initialFormulation = copy.deepcopy(formulation)
    formulation = initialFormulation.toGeometricProgram()
    names = [vbl.name for vbl in formulation.variables_only]
    Fmat = [[] for i in range(0,len(names))]
    gvec = []
    K = []
    Nterms = 0

    for tm in formulation.objective.terms:
        Nterms += 1
        gvec.append(tm.constant.magnitude)
        exps = [0.0 for i in range(0,len(names))]
        for i in range(0,len(tm.variables)):
            vr = tm.variables[i]
            exps[names.index(vr.name)] = tm.exponents[i]
        for i in range(0,len(names)):
            Fmat[i].append(exps[i]) 
    K.append(Nterms)

    for con in formulation.constraints:
        Nterms = 0
        for tm in con.LHS.terms:
            Nterms += 1
            gvec.append(tm.constant.magnitude)
            exps = [0.0 for i in range(0,len(names))]
            for i in range(0,len(tm.variables)):
                vr = tm.variables[i]
                exps[names.index(vr.name)] = tm.exponents[i]
            for i in range(0,len(names)):
                Fmat[i].append(exps[i]) 
        K.append(Nterms)

    F = cvxopt.matrix(Fmat)
    g = cvxopt.log(cvxopt.matrix(gvec))

    sol = cvxopt.solvers.gp(K, F, g)

    if options.returnRaw:
        vrs = list(cvxopt.exp( sol['x'] ))
        for i in range(0,len(vrs)):
            unt = formulation.variables_only[i].units
            vrs[i] = vrs[i] * unt
        obj = cvxopt.exp( sol['primal objective']) * formulation.objective.units
        vDict = dict(zip(names, vrs))
        return [sol,vDict,obj]

    if sol['status'] != 'optimal':
        if sol['status'] == 'unknown':
            unkCons = []
            unkCons.append(sol['relative gap'] <= formulation.solverOptions.unknownTolerance)
            unkCons.append(sol['primal infeasibility'] <= formulation.solverOptions.unknownTolerance)
            unkCons.append(sol['dual infeasibility'] <= formulation.solverOptions.unknownTolerance)
            unkCons.append(sol['primal slack'] <= formulation.solverOptions.unknownTolerance)
            unkCons.append(sol['dual slack'] <= formulation.solverOptions.unknownTolerance)
            if all(unkCons):
                outputMessages.append('Solver returned as unknown, but user specified tolerance of %.4e was met'%(formulation.solverOptions.unknownTolerance))
            else:
                raise RuntimeError("Solver did not return optimal solution.  Check formulation or pass in 'options.returnRaw=True' to bypass output.")

    res = OptimizationOutput()
    res.initialFormulation = copy.deepcopy(initialFormulation)
    res.solveType = options.solveType
    res.solver = options.solver.lower()
    res.solvedFormulation = copy.deepcopy(formulation)
    res.rawResult = copy.deepcopy(sol)

    vrs = list(cvxopt.exp( sol['x'] ))
    fm_vrs = formulation.variables_only
    for i in range(0,len(vrs)):
        unt = fm_vrs[i].units
        vrs[i] = vrs[i] * unt

    obj = cvxopt.exp( sol['primal objective']) * formulation.objective.units
    vDict = dict(zip(names, vrs))

    res.variables = vDict
    res.objective = obj
    for tm in initialFormulation.objective.terms:
        if isinstance(tm.substitute(), (int,float,Q_)):
            res.objective += tm.substitute()

    res.sensitivities = None
    res.messages = outputMessages
    
    return res