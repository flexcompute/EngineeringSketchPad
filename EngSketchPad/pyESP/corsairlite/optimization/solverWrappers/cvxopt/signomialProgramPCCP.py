import copy as copy
from corsairlite.optimization import Variable, Constraint
from corsairlite.core.dataTypes import OptimizationOutput
from corsairlite.optimization.formulation import Formulation
from corsairlite.optimization.solverWrappers import cvxopt as solverTree

def SignomialProgramPCCPSolver(formulation):
    initialFormulation = copy.deepcopy(formulation)
    formulation = formulation.toSignomialProgram(substituted=False)
    # I'm on the fence about defaulting to this or not...  For now leave commented.
    # formulation.solverOptions.solveType = 'sp'
    # # formulation.generateSolverOptions()
    # try: # If this solves as an SP, it will be much faster, so do it.
    #     res = solverTree.SignomialProgramStandardSolver(formulation)
    #     res.messages.append('The [sp] method was used in place of [pccp] method')
    #     return res
    # except:
    #     formulation.solverOptions.solveType = 'pccp'

    slackVariableList = []
    slackCounter = 0
    newConstraints = []

    cons = formulation.constraints
    for con in cons:
        try:
            conSet = con.toPosynomial(GPCompatible = True)
            newConstraints += conSet
        except:
            con = con.toAllPositive()
            slackCounter += 1
            s = Variable(name = "__pccpSlack_%d"%(slackCounter), 
                         guess = 1.0, 
                         units = "-",   
                         description = "Slack number %d for PCCP"%(slackCounter))
            slackVariableList.append(s)
            if con.operator == '>=':
                newCon = Constraint(s * con.LHS, '>=', con.RHS)
                newConstraints += [newCon, s>=1]
            elif con.operator == '<=':
                newCon = Constraint(con.LHS, '<=', s*con.RHS)
                newConstraints += [newCon, s>=1]
            elif con.operator == '==':
                newCon = Constraint(s * con.LHS, '>=', con.RHS)
                newConstraints += [newCon]
                slackCounter += 1
                s2 = Variable(name = "__pccpSlack_%d"%(slackCounter), 
                             guess = 1.0, 
                             units = "-",   
                             description = "Slack number %d for PCCP"%(slackCounter))
                slackVariableList.append(s2)
                newCon = Constraint(con.LHS, '<=', s2*con.RHS)
                newConstraints += [newCon, s >=1, s2 >=1]
                
            else:
                raise ValueError('Invalid constraint operator')
            
    
    newObjective = copy.deepcopy(formulation.objective)

    for vr in slackVariableList:
        newObjective *= vr**formulation.solverOptions.penaltyExponent
    newFormulation = Formulation(newObjective, newConstraints)

    if formulation.solverOptions.x0 is not None:
        for vr in slackVariableList:
            formulation.solverOptions.x0[vr.name] = vr.guess

    # spOptions = copy.deepcopy(options)
    # spOptions.solveType = 'sp'
    # spOptions.returnRaw = False
    newFormulation.solverOptions = copy.deepcopy(formulation.solverOptions)
    newFormulation.solverOptions.solveType = 'sp'
    newFormulation.solverOptions.returnRaw = False
    # newFormulation.generateSolverOptions()
    try:
        res = solverTree.SignomialProgramStandardSolver(newFormulation)
    except Exception as e:
        if formulation.solverOptions.returnRaw:
            newFormulation.solverOptions.returnRaw = True
            res = solverTree.SignomialProgramStandardSolver(newFormulation)
            return res
        else:
            raise e

    for i in range(0,len(res.messages)):
        if res.messages[i] == 'Auto-selected solve type [sp]':
            res.messages[i] = 'Auto-selected solve type [pccp]'

    slackNames = [slackVariableList[i].name for i in range(0,len(slackVariableList))]
    for i in range(0,len(slackVariableList)):
        if abs(res.variables[slackNames[i]]-1) > formulation.solverOptions.slackTol:
            res.messages.append('WARNING: One or more slack tolerances not met (tol=%.2e)'%(formulation.solverOptions.slackTol))

    for i in range(0,len(slackVariableList)):
        res.objective /= res.variables[slackNames[i]]**formulation.solverOptions.penaltyExponent
        del res.variables[slackNames[i]]

    res.solveType = 'pccp'
    res.solver = formulation.solverOptions.solver.lower()
    # res.sensitivities = res.intermediateSolves[-1].sensitivities
            
    return res