import copy as copy
import cvxopt
import numpy as np
from corsairlite import units, Q_
from corsairlite.core.dataTypes import OptimizationOutput

# Required Options:
# options.solver
# options.solveType
# options.returnRaw

def QuadraticProgramSolver(formulation):
    cvxopt.solvers.options['show_progress'] = False
    cvxopt.solvers.options['maxiters'] = 1000
    cvxopt.solvers.options['abstol'] = 1e-7
    cvxopt.solvers.options['reltol'] = 1e-6
    cvxopt.solvers.options['feastol'] = 1e-7
    
    outputMessages = []
    if hasattr(formulation.solverOptions,'autoSelect'):
        if formulation.solverOptions.autoSelect == True:
            outputMessages.append('Auto-selected solve type [qp]')

    options = copy.deepcopy(formulation.solverOptions)
    initialFormulation = copy.deepcopy(formulation)
    formulation = formulation.toQuadraticProgram()
    names = [vbl.name for vbl in formulation.variables_only]

    # range(0,len(formulation.objective.terms))
    pRow = [0.0 for i in range(0,len(names))]
    Praw = [copy.deepcopy(pRow) for i in range(0,len(names))]
    qraw = [0.0 for i in range(0, len(names))]
    for tm in formulation.objective.terms:
        nvar = len(tm.variables)
        if nvar == 1:
            ix = names.index(tm.variables[0].name)
            if tm.exponents[0] == 1.0:
                qraw[ix] = tm.constant.magnitude
            if tm.exponents[0] == 2.0:
                Praw[ix][ix] = tm.constant.magnitude
        if nvar == 2:
            ix1 = names.index(tm.variables[0].name)
            ix2 = names.index(tm.variables[1].name)
            Praw[ix1][ix2] = tm.constant.magnitude/2.0
            Praw[ix2][ix1] = tm.constant.magnitude/2.0
    P = cvxopt.matrix(Praw)
    q = cvxopt.matrix(qraw)
    Graw = None
    Araw = None
    hraw = None
    braw = None
    Iix = 0
    Eix = 0
    for conNum in range(0,len(formulation.constraints)):
        con = formulation.constraints[conNum]
        if con.operator == '<=':
            if Graw is None:
                Graw = [[] for i in range(0,len(names))]
                hraw = []
            for i in range(0,len(Graw)):
                Graw[i].append(0.0)
            for tm in con.LHS.terms:
                ix = names.index(tm.variables[0].name)
                Graw[ix][Iix] = tm.constant.magnitude
            hraw.append(con.RHS.terms[0].constant.magnitude)
            Iix += 1

        if con.operator == '==':
            if Araw is None:
                Araw = [[] for i in range(0,len(names))]
                braw = []
            for i in range(0,len(Araw)):
                Araw[i].append(0.0)
            for tm in con.LHS.terms:
                ix = names.index(tm.variables[0].name)
                Araw[ix][Eix] = tm.constant.magnitude 
            braw.append(con.RHS.terms[0].constant.magnitude)
            Eix += 1      

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
        sol = cvxopt.solvers.qp( 2.0*P, q, G, h, A, b)
    elif A is not None:
        sol = cvxopt.solvers.qp( 2.0*P, q, A = A, b = b)
    elif G is not None:
        sol = cvxopt.solvers.qp( 2.0*P, q, G, h)
    else:
        sol = cvxopt.solvers.qp( 2.0*P, q)

    if options.returnRaw:
        return sol

    if sol['status'] != 'optimal':
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
        unt = fm_vrs[i].units
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