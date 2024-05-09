from corsairlite import units, Q_, _Unit
from collections import OrderedDict
from collections.abc import Mapping
import operator
# from pandas import DataFrame
import numpy as np
import inspect
import copy
# ==========================================================================================================
# ==========================================================================================================
# class AdditiveBreakdown(object):
#     def __init__(self, units):
#         inputUnits = units
#         units = corsairlite.units
#         self.components = {}
#         self.units = inputUnits

#     def __repr__(self):
#         pstr = ''
#         lengths = [len(key) for key in self.components.keys()]
#         maxlen = max(lengths)
#         pstr = pstr + 'Total : %.4f %s\n'%(self.total.to(self.units).magnitude, str(self.units))

#         sorted_x = sorted(self.components.items(), key=operator.itemgetter(1))
#         sorted_dict = OrderedDict(sorted_x)
#         kys = list(reversed(sorted_dict.keys()))
#         vls = list(reversed(sorted_dict.values()))
#         # Does alphabetical instead
#         # kys = sorted(self.components.keys())
#         # vls = [self.components[kys[i]] for i in range(0,len(kys))]
#         pitems = OrderedDict(zip(kys, vls))

#         for key,val in pitems.iteritems():
#             for i in range(0,maxlen - len(key) + 4):
#                 pstr = pstr + ' '
#             pstr = pstr + '%s : %.4f %s\n'%(key,val.to(self.units).magnitude,str(self.units))

#         return pstr

#     @property
#     def total(self):
#         total = 0.0
#         for key,val in self.components.iteritems():
#             total = total + val
#         return total
# ==========================================================================================================
# ==========================================================================================================
# class MultiplicativeBreakdown(object):
#     def __init__(self, units):
#         inputUnits = units
#         units = corsairlite.units
#         self.components = {}
#         self.units = inputUnits

#     def __repr__(self):
#         pstr = ''
#         lengths = [len(key) for key in self.components.keys()]
#         maxlen = max(lengths)
#         pstr = pstr + 'Total : %.4f %s\n'%(self.total.to(self.units).magnitude, str(self.units))

#         sorted_x = sorted(self.components.items(), key=operator.itemgetter(1))
#         sorted_dict = OrderedDict(sorted_x)
#         kys = sorted_dict.keys()
#         vls = sorted_dict.values()
#         # Does alphabetical instead
#         # kys = sorted(self.components.keys())
#         # vls = [self.components[kys[i]] for i in range(0,len(kys))]
#         pitems = OrderedDict(zip(kys, vls))

#         for key,val in pitems.iteritems():
#             for i in range(0,maxlen - len(key) + 4):
#                 pstr = pstr + ' '
#             pstr = pstr + '%s : %.4f %s\n'%(key,val.to(self.units).magnitude,str(self.units))

#         return pstr

#     @property
#     def total(self):
#         total = 1.0
#         for key,val in self.components.iteritems():
#             total = total * val
#         return total
# ==========================================================================================================
# ==========================================================================================================
def determineCategory(val):
    # if hasattr(val, '__dict__'):
    if isinstance(val, Container):
        return 'Container'
    if isinstance(val, str):
        return 'String'
    if isinstance(val, bool):
        return 'Boolean'
    if isinstance(val, float) or isinstance(val, int):
        return 'Scalar'
    if isinstance(val, units.Quantity):
        if hasattr(val.magnitude, "__len__"):
            return 'Vector'
        else:
            return 'Scalar'
    else:
        if hasattr(val, "__len__"):
            return 'Vector'

    # raise ValueError('Could not categorize')

class Container(object):
    def __init__(self):
        pass

    def __eq__(self, other):
        if not isinstance(other, Container):
            raise ValueError('Comparison between Container object and non-Container object is not defined')

        keys1 = self.__dict__.keys()
        keys2 = other.__dict__.keys()

        if len(keys1) != len(keys2):
            return False

        for ky1 in keys1:
            if ky1 not in keys2:
                return False
        for ky2 in keys2:
            if ky2 not in keys1:
                return False

        # At this point, we know both containers have the same keys in them
        for ky in keys1:
            val1 = self.__dict__[ky]
            val2 = other.__dict__[ky]

            try:
                val1 == val2
                if val1 != val2:
                    return False
                couldNotDetermine = False
            except:
                couldNotDetermine = True

            if couldNotDetermine:
                cat1 = determineCategory(val1)
                cat2 = determineCategory(val2)
                if cat1 != cat2:
                    return False

                # Both values are known to be of the same type
                cat = cat1
                if cat == 'String':
                    if val1 != val2:
                        return False
                if cat == 'Boolean':
                    if val1 != val2:
                        return False
                if cat == 'Vector':
                    try:
                        chk = val1 != val2
                        if chk.any():
                            return False
                    except:
                        chk = val1 != val2
                        if chk:
                            return False
                if cat == 'Scalar':
                    if val1 != val2:
                        return False
                if cat == 'Container':
                    if val1 != val2:
                        return False
        return True

    def __ne__(self, other):
        return not self.__eq__(other)
# ==========================================================================================================
# ==========================================================================================================
# class SweepTable(DataFrame):
#     def __init__(self, data=None, index=None, columns=None, dtype=None, copy=False):
#         super(SweepTable, self).__init__(data, index, columns, dtype, copy)

#     def cutSlice(self, independent, dependent, cutValues):
#         dfTemp = self
#         for key, val in cutValues.iteritems():
#             dfTemp = dfTemp.loc[dfTemp[key] == val]

#         dfTemp = dfTemp.sort_values(by=independent)

#         res = {}
#         if isinstance(independent,(list,tuple)):
#             for idp in independent:
#                 res[idp] = dfTemp[idp].values
#         else:
#             res[independent] = dfTemp[independent].values

#         if isinstance(dependent,(list,tuple)):
#             for dep in dependent:
#                 res[dep] = dfTemp[dep].values
#         else:
#             res[dependent] = dfTemp[dependent].values

#         return res
# ==========================================================================================================
# ==========================================================================================================
# class TypeCheckedDict(dict):
#     def __init__(self, checkItem, itemDict = None):
#         super(TypeCheckedDict, self).__init__()
#         self.checkItem = checkItem

#         if itemDict is not None:
#             if isinstance(itemDict, dict):
#                 for ky, vl in itemDict.iteritems():
#                     print ky
#                     print vl
#                     self.__setitem__(ky,vl)
#             else:
#                 raise ValueError('Input to itemDict is not iterable')

#     def __setitem__(self, key, val):
#         if isinstance(val, self.checkItem):
#             super(TypeCheckedDict, self).__setitem__(key, val)
#         else:
#             raise ValueError('Input must be an instance of the defined type')

#     def __getitem__(self, key):
#         return super(TypeCheckedDict, self).__getitem__(key)

#     def update(self, other=None, **kwargs):
#         if other is not None:
#             for k, v in other.items() if isinstance(other, Mapping) else other:
#                 self.__setitem__(k,v)
#         for k, v in kwargs.items():
#             self.__setitem__(k,v)
# ==========================================================================================================
# ==========================================================================================================
class TypeCheckedList(list):
    def __init__(self, checkItem, itemList = None):
        super(TypeCheckedList, self).__init__()
        self.checkItem = checkItem

        if itemList is not None:
            if isinstance(itemList, list) or isinstance(itemList, tuple):
                for itm in itemList:
                    self.append(itm)
            else:
                raise ValueError('Input to itemList is not iterable')

    def __setitem__(self, key, val):
        if isinstance(val, self.checkItem):
            super(TypeCheckedList, self).__setitem__(key, val)
        else:
            raise ValueError('Input must be an instance of the defined type')

    def __setslice__(self, i, j, sequence):
        performSet = True
        for val in sequence:
            if not isinstance(val, self.checkItem):
                performSet = False
                break

        if performSet:
            super(TypeCheckedList, self).__setslice__(i,j,sequence)
        else:
            raise ValueError('All values in the input must be an instance of the defined type')

    def append(self, val):
        if isinstance(val, self.checkItem):
            super(TypeCheckedList, self).append(val)
        else:
            raise ValueError('Input must be an instance of the defined type')
# ==========================================================================================================
# ==========================================================================================================
class ValueCheckedList(list):
    def __init__(self, expectedLength, expectedUnits, itemList = None):
        super(ValueCheckedList, self).__init__()
        if isinstance(expectedUnits, str):
            ut = getattr(units, expectedUnits)
        elif isinstance(expectedUnits, _Unit):
            ut = expectedUnits
        else:
            raise ValueError('Units not valid')

        self.expectedLength = expectedLength
        self.expectedUnits = ut

        if itemList is not None:
            if isinstance(itemList,np.ndarray):
                itemList = itemList * units.dimensionless

            if isinstance(itemList, Q_):
                if hasattr(itemList.magnitude, '__len__'):
                    tmp = []
                    for itm in itemList:
                        tmp.append(itm)
                    itemList = tmp

            if isinstance(itemList, list) or isinstance(itemList, tuple):
                for itm in itemList:
                    self.append(itm)
            else:
                raise ValueError('Input to itemList is not a valid iterable')

    def __eq__(self, other):
        if not isinstance(other, ValueCheckedList):
            raise ValueError('Comparison between ValueCheckedList and non-ValueCheckedList is not defined')

        if len(other) != len(self):
            return False

        for i in range(0,len(self)):
            chk = other[i] != self[i]
            try:
                if chk.any():
                    return False
            except:
                if chk:
                    return False

        return True

    def __ne__(self, other):
        return not self.__eq__(other)

    def __setitem__(self, key, val):
        val = self.correctUnits(val)
        if self.checkValid(val):
            super(ValueCheckedList, self).__setitem__(key, val)
        else:
            raise ValueError('The input does not match either the expected units or the expected length')

    def __setslice__(self, i, j, sequence):
        performSet = True
        correctedSequence = []
        for val in sequence:
            val = self.correctUnits(val)
            if not self.checkValid(val):
                performSet = False
                break
            correctedSequence.append(val)

        if performSet:
            super(ValueCheckedList, self).__setslice__(i,j,correctedSequence)
        else:
            raise ValueError('The input does not match either the expected units or the expected length')

    def append(self, val):
        val = self.correctUnits(val)
        if self.checkValid(val):
            super(ValueCheckedList, self).append(val)
        else:
            raise ValueError('The input does not match either the expected units or the expected length')

    def correctUnits(self, val):
        if not isinstance(val, Q_):
            val = val * units.dimensionless

        unt = 1.0 * self.expectedUnits
        unt = unt.to_base_units()
        valtmp = val.to_base_units()
        if valtmp.units == unt.units:
            val = val.to(self.expectedUnits)
        else:
            raise ValueError('Input has invalid units')

        return val

    def checkValid(self, val):
        correctUnits = False
        correctLength = False

        if self.expectedLength == None:
            correctLength = True

        if isinstance(val, Q_):
            if val.units == self.expectedUnits:
                correctUnits = True
        else:
            if units.dimensionless == self.expectedUnits:
                correctUnits = True

        if not hasattr(val, '__len__'):
            if self.expectedLength == 0:
                correctLength = True

        if isinstance(val, Q_):
            if not hasattr(val.magnitude, '__len__'):
                if self.expectedLength == 0:
                    correctLength = True
            else:
                vec = val.magnitude
        elif hasattr(val, '__len__'):
            vec = val
        else:
            pass

        if len(vec) == self.expectedLength:
            correctLength = True

        return correctUnits and correctLength
# ==========================================================================================================
# ==========================================================================================================
# class VariableContainer(Container):
#     def __init__(self):
#         super(VariableContainer, self).__init__()

#     @property
#     def NumberOfVariables(self):
#         dct = self.VariableDictionary

#         Nvars = 0
#         for key, val in dct.iteritems():
#             key = key.replace('_true','')
#             isVector  = False
#             isScalar  = False
#             isString  = False
#             isBoolean = False

#             if isinstance(val, str):
#                 isString = True
#             if isinstance(val, bool):
#                 isBoolean = True
#             if isinstance(val, float) or isinstance(val, int):
#                 isScalar = True
#             if isinstance(val, units.Quantity):
#                 if hasattr(val.magnitude, "__len__"):
#                     isVector = True
#                 else:
#                     isScalar = True
#             else:
#                 if hasattr(val, "__len__"):
#                     isVector = True

#             if isString:
#                 Nvars += 1
#             if isBoolean:
#                 Nvars += 1
#             if isVector and not isString:
#                 vec = np.asarray(val)
#                 Nvars += vec.size
#             if isScalar and not isBoolean:
#                 Nvars += 1

#         return Nvars

#     @property
#     def VariableDictionary(self):
#         dct = self.__dict__
#         rtnDct = {}

#         for key, val in dct.iteritems():
#             key = key.replace('_true','')
#             isVector  = False
#             isScalar  = False
#             isString  = False
#             isBoolean = False
#             isClass   = False

#             if hasattr(val, '__dict__'):
#                 for k, v in val.__dict__.iteritems():
#                     if isinstance(v, VariableContainer):
#                         isClass = True
#             if isinstance(val, str):
#                 isString = True
#             if isinstance(val, bool):
#                 isBoolean = True
#             if isinstance(val, float) or isinstance(val, int):
#                 isScalar = True
#             if isinstance(val, units.Quantity):
#                 if hasattr(val.magnitude, "__len__"):
#                     isVector = True
#                 else:
#                     isScalar = True
#             else:
#                 if hasattr(val, "__len__"):
#                     isVector = True

#             if isString:
#                 rtnDct[key] = val
#             if isBoolean:
#                 rtnDct[key] = val
#             if isVector and not isString:
#                 rtnDct[key] = val
#             if isScalar and not isBoolean:
#                 rtnDct[key] = val
#             if isClass:
#                 for k, v in val.__dict__.iteritems():
#                     if isinstance(v, VariableContainer):
#                         classDct = v.VariableDictionary
#                         break
#                 interDct = {}
#                 for ik, iv in classDct.iteritems():
#                     topLevelKey = key + '.' + ik
#                     interDct[topLevelKey] = iv
#                 rtnDct.update(interDct)

#         return OrderedDict(sorted(rtnDct.items()))
# ==========================================================================================================
# ==========================================================================================================
class OptimizationOutput(object):
    def __init__(self):
        self.initialFormulation = None
        self.solveType = None
        self.solver = None
        self.solvedFormulation = None
        self.rawResult = None
        self.objective = None
        self.variables = None

    def computeSensitivities(self):
        from corsairlite.optimization.computeSensitivities import computeSensitivities
        computeSensitivities(self)


    def result(self, Ndecimal=2):
        solveDict = {
                        'lp':'Linear Program',
                        'qp':'Quadratic Program',
                        'gp':'Geometric Program',
                        'sp':'Signomial Program--Standard',
                        'pccp':'Signomial Program--Penalty Convex Concave Programming',
                        'ga':'Genetic Algorithm',
                        'nlp': 'Non-Linear Program',
                        'pso': 'Particle Swarm Optimization',
                        'sa': 'Simulated Annealing',
                        'nms': 'Nelder-Mead Simplex',
                        'sqp': 'Sequential Quadratic Program',
                        'lsqp': 'Logspace Sequential Quadratic Program',
                        'slcp': 'Sequential Log Convex Program'
        }

        fstr = '%.' + str(Ndecimal) + 'f'
        pstr = ''
        pstr += '\n'
        # pstr += '\n'
        pstr += 'Objective\n'
        pstr += '---------\n'
        pstr += '   ' + fstr%(self.objective.magnitude) + ' ' + str(self.objective.units) + '\n'
        pstr += '\n'
        pstr += 'Variables\n'
        pstr += '---------\n'
        # names = [self.solvedFormulation.variables[i].name for i in range(0,len(self.solvedFormulation.variables))]
        namesUnsorted = list(self.variables.keys())
        names = sorted(namesUnsorted)
        vls = [self.variables[names[i]].magnitude for i in range(0,len(names))]
        for i in range(0,len(vls)):
            v = vls[i]
            if float(v).is_integer():
                vls[i] = '%d'%(v) + ' '*(1+Ndecimal)
            else:
                vls[i] = fstr%(v)
        uts = [str(self.variables[names[i]].units) for i in range(0,len(names))]
        for i in range(0,len(uts)):
            if uts[i] == 'dimensionless':
                uts[i] = '[-]'

        # allVars = copy.deepcopy(self.solvedFormulation.variables)
        allVars = copy.deepcopy(self.initialFormulation.variables_only)
        avnms = [allVars[i].name for i in range(0,len(allVars))]
        ifVars = []
        for nm in names:
            tfl = [avnms[i] == nm for i in range(0,len(avnms))]
            ix = tfl.index(True)
            ifVars.append(allVars[ix])
        varDict = dict(zip(names, ifVars))
        descrs = [varDict[names[i]].description for i in range(0,len(names))]

        widths = []
        widths.append(max([len(s) for s in names]))
        widths.append(max([len(s) for s in vls]))
        widths.append(max([len(s) for s in uts]))
        widths.append(max([len(s) for s in descrs]))

        for i in range(0,len(names)):
            pstr += '   '
            pstr += names[i].ljust(widths[0])
            pstr += '  :  '
            pstr += vls[i].rjust(widths[1])
            pstr += '   '
            pstr += uts[i].center(widths[2])
            pstr += '   '
            pstr += descrs[i].ljust(widths[3])
            pstr += '\n'
        pstr += '\n'

        from corsairlite.optimization.constant import Constant
        vrs = [self.initialFormulation.variables[i] for i in range(0,len(self.initialFormulation.variables))]
        csts = []
        for v in vrs:
            if isinstance(v, Constant):
                csts.append(v)
        if len(csts) != 0:
            pstr += 'Constants\n'
            pstr += '---------\n'
            cstNames = [csts[i].name for i in range(0,len(csts))]
            vls = [csts[i].value.magnitude for i in range(0,len(csts))]
            for i in range(0,len(vls)):
                v = vls[i]
                if float(v).is_integer():
                    vls[i] = '%d'%(v) + ' '*(1+Ndecimal)
                else:
                    vls[i] = fstr%(v)
            uts = [str(csts[i].units) for i in range(0,len(csts))]
            for i in range(0,len(uts)):
                if uts[i] == 'dimensionless':
                    uts[i] = '[-]'
            descrs = [csts[i].description for i in range(0,len(csts))]

            widths = []
            widths.append(max([len(s) for s in cstNames]))
            widths.append(max([len(s) for s in vls]))
            widths.append(max([len(s) for s in uts]))
            widths.append(max([len(s) for s in descrs]))

            for i in range(0,len(csts)):
                pstr += '   '
                pstr += cstNames[i].ljust(widths[0])
                pstr += '  :  '
                pstr += vls[i].rjust(widths[1])
                pstr += '   '
                pstr += uts[i].center(widths[2])
                pstr += '   '
                pstr += descrs[i].ljust(widths[3])
                pstr += '\n'
            pstr += '\n'

        pstr += 'Sensitivities\n'
        pstr += '-------------\n'

        if len(csts) == 0 and self.solveType in ['gp','lp','qp','sp','pccp']:
            pstr += '   No Constants in the Formulation\n'
        elif self.solveType in ['gp','lp','qp','sp','pccp'] and self.sensitivities is None:
            pstr += '   Sensitivities not computed\n'
        elif len(csts) != 0 and self.solveType in ['gp','lp','qp','sp','pccp']:
            # pstr += 'Sensitivities\n'
            # pstr += '-------------\n'
            # if list(self.sensitivities.keys()) == ['__error_message']:
            #     if isinstance(self.sensitivities['__error_message'],str):
            #         pstr += self.sensitivities['__error_message']
            # else:
            cstNames = [csts[i].name for i in range(0,len(csts))]
            sensVals = [self.sensitivities[cN] for cN in cstNames]
            odr = np.flipud(np.argsort(abs(np.array(sensVals)))).tolist()
            cstNames = [cstNames[i] for i in odr]
            sensVals = [sensVals[i] for i in odr]
            descrs   = [csts[i].description for i in odr]

            for i in range(0,len(sensVals)):
                v = sensVals[i]
                if v > 10**(-Ndecimal):
                    sgn = '+ '
                elif v < -10**(-Ndecimal):
                    sgn = '- '
                else:
                    sgn = '  '
                if float(v).is_integer():
                    sensVals[i] = sgn+'%d'%(abs(v)) + ' '*(1+Ndecimal)
                else:
                    sensVals[i] = sgn+fstr%(abs(v))

            widths = []
            widths.append(max([len(s) for s in cstNames]))
            widths.append(max([len(s) for s in sensVals]))
            widths.append(max([len(s) for s in descrs]))

            for i in range(0,len(csts)):
                pstr += '   '
                pstr += cstNames[i].ljust(widths[0])
                pstr += '  :  '
                pstr += sensVals[i].rjust(widths[1])
                pstr += '   '
                pstr += descrs[i].ljust(widths[2])
                pstr += '\n'
        elif self.solveType in ['sqp','lsqp']:
            pstr += '   Sensitivities not currently implemented for optimization type %s \n'%(solveDict[self.solveType])
        else:
            pstr += '   Sensitivities not available for optimizations of type %s \n'%(solveDict[self.solveType])
        pstr += '\n'

        pstr += 'Solve Report\n'
        pstr += '------------\n'
        if self.solveType in solveDict:
            pstr += '   Solve Method           :  %s \n'%(solveDict[self.solveType])
        else:
            pstr += '   Solve Method           :  %s \n'%('Not Specified')

        if self.solveType in ['gp','lp','qp']:
            cs = 'Convex'
        elif self.solveType in ['sp','pccp','sqp','lsqp','slcp']:
            cs = 'Iterative'
        elif self.solveType in ['ga','pso','sa','nms']:
            cs = 'Heuristic'
        else:
            cs = 'Unknown'
        pstr += '   Classification         :  %s \n'%(cs)

        pstr += '   Solver                 :  %s \n'%(self.solver)

        if self.solveType in ['sp','pccp','nms','sqp','lsqp','slcp']:
            pstr += '   Number of Iterations   :  %s \n'%(self.numberOfIterations)

        if self.solveType == 'ga':
            pstr += '   Number of Generations  :  %s \n'%(self.numberOfGenerations)
            pstr += '   Population Size        :  %s \n'%(self.populationSize)

        if self.solveType == 'pso':
            pstr += '   Number of Time Steps   :  %s \n'%(self.numberOfTimeSteps)
            pstr += '   Swarm Size             :  %s \n'%(self.swarmSize)

        if self.solveType == 'sa':
            pstr += '   Number of Trials       :  %s \n'%(self.numberOfTrials)

        if self.solveType in ['sp','pccp','ga','pso','nms','sqp','lsqp','slcp']:
            pstr += '   Termination Status     :  %s \n'%(self.terminationStatus)

        if self.solveType in ['sa','ga','pso','nms']:
            pstr += '   Feasibility            :  %s \n'%(self.feasibility)

        if hasattr(self,'messages'):
            if len(self.messages) > 0:
                pstr += '\n'
                pstr += 'Messages\n'
                pstr += '--------\n'
                nsp = len(str(len(self.messages)))
                for i in range(0,len(self.messages)):
                    msg = self.messages[i]
                    for ky,vl in solveDict.items():
                        msg = msg.replace('['+ky+']',vl)
                    pstr += '   %s. '%(str(i+1).rjust(nsp)) + msg
                    pstr += '\n'





        # if self.solveType in ['ga','pso','nms']:
        #     pstr += '   Termination Status  :  %s \n'%(self.terminationStatus)
        #     pstr += '   Feasibility         :  %s \n'%(self.feasibility)
        #     pstr += '\n'

        # if self.solveType == 'sa':
        #     pstr += '   Feasibility  :  %s \n'%(self.feasibility)
        #     pstr += '\n'

        # pstr += 'Solve Type\n'
        # pstr += '----------\n'
        # if self.solveType in solveDict:
        #     pstr += '   Solve Method  :  %s \n'%(solveDict[self.solveType])
        # else:
        #     pstr += '   Solve Method  :  %s \n'%('Not Specified')
        # pstr += '   Solver        :  %s'%(self.solver)
        # pstr += '\n'
        # pstr += '\n'

        return pstr
# ==========================================================================================================
# ==========================================================================================================
