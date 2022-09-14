import copy as copy
import corsairlite
from corsairlite import units, Q_, _Unit
from corsairlite.optimization.term import Term
from corsairlite.optimization.constant import Constant
from corsairlite.optimization.formulation import Formulation
from corsairlite.core.dataTypes import OptimizationOutput, Container

# =====================================================================================================================
# The solve function
# =====================================================================================================================
def solve(formulation):
    formulation.generateSolverOptions()
    # options = copy.deepcopy(options)
# =====================================================================================================================
# Initial setup
# =====================================================================================================================
    if formulation.solverOptions.solveType is None: 
        formulation.solverOptions.autoSelect = True
        try:
            f = formulation.toLinearProgram()
            formulation.solverOptions.solveType = 'lp'
        except:
            try:
                f = formulation.toQuadraticProgram()
                formulation.solverOptions.solveType = 'qp'
            except:
                try:
                    f = formulation.toGeometricProgram()  
                    formulation.solverOptions.solveType = 'gp'
                except:
                    try:
                        f = formulation.toSignomialProgram()
                        formulation.solverOptions.solveType = 'pccp' # will automatically try sp
                    except:
                        formulation.solverOptions.solveType = 'sqp'
        # solveType = 'sqp'
        # solveType = 'sgp'

    formulation.generateSolverOptions()

    if formulation.solverOptions.solver is None:
        if formulation.solverOptions.solveType in ['lp','qp','gp','sp','pccp','sqp','lsqp','slcp']:
            formulation.solverOptions.solver = 'cvxopt'
# =====================================================================================================================
# Solve Linear Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'lp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.LinearProgramSolver(formulation)
            return res
# =====================================================================================================================
# Solve Quadratic Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'qp':
        if formulation.solverOptions.solver.lower() == 'cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.QuadraticProgramSolver(formulation)
            return res
# =====================================================================================================================
# Solve Geometric Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'gp':
        import numpy as np
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.GeometricProgramSolver(formulation)
            return res
# =====================================================================================================================
# Solve Signomial Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'sp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.SignomialProgramStandardSolver(formulation)
            return res
# =====================================================================================================================
# Solve Signomial Program using the Penalty Convex-Concave Programming Algorithm
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'pccp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.SignomialProgramPCCPSolver(formulation)
            return res
# =====================================================================================================================
# Solve Sequential Quadratic Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'sqp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.SequentialQuadraticProgramSolver(formulation)
            return res
# =====================================================================================================================
# Solve Logspace Sequential Quadratic Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'lsqp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.LogspaceSequentialQuadraticProgramSolver(formulation)
            return res
# =====================================================================================================================
# Solve Sequential Log Convex Program
# =====================================================================================================================
    if formulation.solverOptions.solveType.lower() == 'slcp':
        if formulation.solverOptions.solver.lower() =='cvxopt':
            import corsairlite.optimization.solverWrappers.cvxopt as cvxopt
            res = cvxopt.SequentialLogConvexProgramSolver(formulation)
            return res



            