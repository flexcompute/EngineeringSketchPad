import copy as copy
import cvxopt
import numpy as np
from corsairlite import units, Q_
from corsairlite.core.dataTypes import OptimizationOutput

# Required Options:
# options.solver
# options.solveType
# options.returnRaw

def LinearProgramSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7

    outputMessages = []
    if hasattr(formulation.solverOptions,'autoSelect'):
        if formulation.solverOptions.autoSelect == True:
            outputMessages.append('Auto-selected solve type [lp]')

    options = copy.deepcopy(formulation.solverOptions)
    initialFormulation = copy.deepcopy(formulation)
    formulation = formulation.toLinearProgram()
    names = [vbl.name for vbl in formulation.variables_only]

    c_raw = [0.0 for i in range(0,len(names))]

    for tm in formulation.objective.terms:
        ix = names.index(tm.variables[0].name)
        if tm.exponents[0] == 1.0:
            c_raw[ix] = tm.constant.magnitude

    c = cvxopt.matrix(c_raw)

    G_raw = []
    h_raw = []

    A_raw = []
    b_raw = []

    for i in range(0,len(formulation.constraints)):
        con = formulation.constraints[i]
        if con.operator == '==':
            Arow = [0.0 for i in range(0,len(names))]
            for tm in con.LHS.terms:
                ix = names.index(tm.variables[0].name)
                Arow[ix] = tm.constant.magnitude
            A_raw.append(Arow)
            b_raw.append(con.RHS.terms[0].constant.magnitude)
        else:
            Grow = [0.0 for i in range(0,len(names))]
            for tm in con.LHS.terms:
                ix = names.index(tm.variables[0].name)
                Grow[ix] = tm.constant.magnitude
            G_raw.append(Grow)
            h_raw.append(con.RHS.terms[0].constant.magnitude)

    Gnew = [[0.0 for i in range(0,len(G_raw))] for ii in range(0,len(G_raw[0]))]
    for i in range(0,len(G_raw)):
        for j in range(0,len(G_raw[0])):
            Gnew[j][i] = G_raw[i][j]
            
    G = cvxopt.matrix(Gnew)
    h = cvxopt.matrix(h_raw)

    if len(A_raw) > 0:
        Anew = [[0.0 for i in range(0,len(A_raw))] for ii in range(0,len(A_raw[0]))]
        for i in range(0,len(A_raw)):
            for j in range(0,len(A_raw[0])):
                Anew[j][i] = A_raw[i][j]
            
        A = cvxopt.matrix(Anew)
        b = cvxopt.matrix(b_raw)

        sol = cvxopt.solvers.lp( c, G, h, A, b)
    else:
        sol = cvxopt.solvers.lp( c, G, h)

    if options.returnRaw:
        return sol

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

    vrs = list(sol['x'])
    fm_vrs = formulation.variables_only
    for i in range(0,len(vrs)):
        unt =fm_vrs[i].units
        vrs[i] = vrs[i] * unt

    obj = sol['primal objective'] * formulation.objective.units
    vDict = dict(zip(names, vrs))

    res.variables = vDict
    res.objective = obj
    for tm in initialFormulation.objective.terms:
        if isinstance(tm.substitute(), (int,float,Q_)):
            res.objective += tm.substitute()
    # -----------------------------------------
    res.sensitivities = None
    res.messages = outputMessages
    return res