import copy as copy
import numpy as np

def computeSensitivities(res):
# =====================================================================================================================
# CVXOPT
# =====================================================================================================================
    if res.solver == 'cvxopt':
    # -----------------------------------------------------------------------------------------------------------------
    # LP
    # -----------------------------------------------------------------------------------------------------------------
        if res.solveType == 'lp':
            sol = res.rawResult 
            vDict = res.variables
            initialFormulation = res.initialFormulation
            unsubbedLP = initialFormulation.toLinearProgram(substituted = False)

            cons   = copy.deepcopy(unsubbedLP.inequalityConstraints) + copy.deepcopy(unsubbedLP.equalityConstraints)
            objULP = copy.deepcopy(unsubbedLP.objective)
            cons = [objULP == 0*objULP.units] + cons
            expressions = [c.LHS for c in cons]
            la = np.array( list(sol['z']) + list(sol['y']) )
            if len(la) == len(expressions) - 1:
                # solver did not include objective sensitivity
                la = np.append(1.0, la)

            term_sensitivity = []
            for i in range(0,len(expressions)):
                expr = expressions[i]
                rw = []
                for tm in unsubbedLP.terms:
                    checkList = []
                    for ii in range(0,len(expr.terms)):
                        checkList.append(tm.isEqual(expr.terms[ii]))
                    if any(checkList):
                        p = la[i]
                        rw.append(p)
                    else:
                        rw.append(0.0)
                term_sensitivity.append(rw)

            term_sensitivity = np.array(term_sensitivity)

            sensitivityDict = {}
            for cst in unsubbedLP.constants:
                cstRow = []
                for tm in unsubbedLP.terms:
                    cstRow.append(tm.sensitivity(cst,vDict).magnitude)

                dterm_dcst = [cstRow for i in range(0,len(expressions))]
                dterm_dcst = np.array(dterm_dcst)
                variableSensitivity = term_sensitivity * dterm_dcst
                sens = np.sum(np.sum(variableSensitivity, axis=1))
                sensitivityDict[cst.name] = sens * cst.value.magnitude / res.objective.magnitude
            res.sensitivities = sensitivityDict
    # -----------------------------------------------------------------------------------------------------------------
    # QP
    # -----------------------------------------------------------------------------------------------------------------
        if res.solveType == 'qp':
            sol = res.rawResult 
            vDict = res.variables
            initialFormulation = res.initialFormulation
            unsubbedQP = initialFormulation.toQuadraticProgram(substituted = False)

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
            if len(la) == len(expressions) - 1:
                # solver did not include objective sensitivity
                la = np.append(1.0, la)

            term_sensitivity = []
            for i in range(0,len(expressions)):
                expr = expressions[i]
                rw = []
                for tm in unsubbedQP.terms:
                    checkList = []
                    for ii in range(0,len(expr.terms)):
                        checkList.append(tm.isEqual(expr.terms[ii]))
                    if any(checkList):
                        p = la[i]
                        rw.append(p)
                    else:
                        rw.append(0.0)
                term_sensitivity.append(rw)

            term_sensitivity = np.array(term_sensitivity)

            sensitivityDict = {}
            for cst in unsubbedQP.constants:
                cstRow = []
                for tm in unsubbedQP.terms:
                    cstRow.append(tm.sensitivity(cst,vDict).magnitude)

                dterm_dcst = [cstRow for i in range(0,len(expressions))]
                dterm_dcst = np.array(dterm_dcst)
                variableSensitivity = term_sensitivity * dterm_dcst
                sens = np.sum(np.sum(variableSensitivity, axis=1))
                sensitivityDict[cst.name] = sens * cst.value.magnitude / res.objective.magnitude
            res.sensitivities = sensitivityDict
    # -----------------------------------------------------------------------------------------------------------------
    # GP
    # -----------------------------------------------------------------------------------------------------------------
        if res.solveType == 'gp':
            sol = res.rawResult 
            vDict = res.variables
            initialFormulation = res.initialFormulation
            unsubbedGP = initialFormulation.toGeometricProgram(substituted = False)

            cons   = copy.deepcopy(unsubbedGP.constraints)
            objUGP = copy.deepcopy(unsubbedGP.objective)
            cons = [objUGP == 0*objUGP.units] + cons
            expressions = [c.LHS for c in cons]
            la = sol['znl']
            if len(la) == len(expressions) - 1:
                # solver did not include objective sensitivity
                la = np.append(1.0, la)
            
            usgpConstants = unsubbedGP.constants
            usgpTerms = unsubbedGP.terms

            term_sensitivity = []
            for i in range(0,len(expressions)):
                expr = expressions[i]
                exprTerms = expr.terms
                rw = [0.0] * len(usgpTerms)
                for tm_idx in range(0,len(usgpTerms)):
                    tm = usgpTerms[tm_idx]
                    checkList = [False] * len(exprTerms)
                    for ii in range(0,len(exprTerms)):
                        checkList[ii] = tm.isEqual(exprTerms[ii])
                    if any(checkList):
                        p = la[i]
                        e = expr.substitute(vDict)
                        av = p/e
                        rw[tm_idx] = av.magnitude
                term_sensitivity.append(rw)
            term_sensitivity = np.array(term_sensitivity)

            sensitivityDict = {}
            for cst in usgpConstants:
                cstRow = []
                for tm in usgpTerms:
                    cstRow.append(tm.sensitivity(cst,vDict).magnitude)

                dterm_dcst = [cstRow for i in range(0,len(expressions))]
                dterm_dcst = np.array(dterm_dcst)
                variableSensitivity = term_sensitivity * dterm_dcst
                sens = np.sum(np.sum(variableSensitivity, axis=1)) 
                sensitivityDict[cst.name] = sens * cst.value.magnitude
            res.sensitivities = sensitivityDict
        # -----------------------------------------------------------------------------------------------------------------
        # SP and PCCP
        # -----------------------------------------------------------------------------------------------------------------
        if res.solveType in ['sp','pccp']:
            unsubbedSP = res.solvedFormulation
            signomialIndicies = res.sensitivityInfo.signomialIndicies
            signomialIndiciesNC = res.sensitivityInfo.signomialIndiciesNC
            lookupTable = res.sensitivityInfo.lookupTable
            designPoint = res.sensitivityInfo.designPoint
            sol = res.sensitivityInfo.sol
            cons = []
            for cix in range(0,len(unsubbedSP.constraints)):
                if cix in signomialIndicies:
                    dependentRows = lookupTable[cix+1]
                    con = unsubbedSP.constraints[cix]
                    newConstraintsIndex = np.where(np.array(signomialIndiciesNC)==cix)[0]
                    conSet = con.approximateConvexifiedSignomial(designPoint)
                    if len(conSet) > 1:
                        raise ValueError('Monomial Approximation of a signomial is returning via an unexpected method')
                    csl = conSet[0].toPosynomial(GPCompatible = True)
                    for cn in csl:
                        cons.append(cn)
                else:
                    cons.append(unsubbedSP.constraints[cix])
            objUSP = copy.deepcopy(unsubbedSP.objective)
            cons = [objUSP == 0*objUSP.units] + cons
            expressions = [c.LHS for c in cons]
            la = sol['znl']
            if len(la) == len(expressions) - 1:
                # solver did not include objective sensitivity
                la = np.append(1.0, la)

            term_sensitivity = []
            for i in range(0,len(expressions)):
                expr = expressions[i]
                rw = []
                for tm in unsubbedSP.terms:
                    checkList = []
                    for ii in range(0,len(expr.terms)):
                        checkList.append(tm.isEqual(expr.terms[ii]))
                    if any(checkList):
                        p = la[i]
                        e = expr.substitute(designPoint)
                        av = p/e
                        rw.append(av.magnitude)
                    else:
                        rw.append(0.0)
                term_sensitivity.append(rw)

            term_sensitivity = np.array(term_sensitivity)

            sensitivityDict = {}
            for cst in unsubbedSP.constants:
                cstRow = []
                for tm in unsubbedSP.terms:
                    cstRow.append(tm.sensitivity(cst,designPoint).magnitude)

                dterm_dcst = [cstRow for i in range(0,len(expressions))]
                dterm_dcst = np.array(dterm_dcst)
                variableSensitivity = term_sensitivity * dterm_dcst
                sens = np.sum(np.sum(variableSensitivity, axis=1)) 
                sensitivityDict[cst.name] = sens * cst.value.magnitude
            res.sensitivities = sensitivityDict
# =====================================================================================================================
# MOSEK
# =====================================================================================================================
    if res.solver == 'mosek':
        pass