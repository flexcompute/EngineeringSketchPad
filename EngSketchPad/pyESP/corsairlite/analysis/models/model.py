import copy
from corsairlite import units, Q_
from corsairlite.core.dataTypes import  TypeCheckedList, Variable

errorString = 'This function is calling to the base Model class and has not been defined.'

class Model(object):

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def __init__(self):
        
        # List of the inputs and outputs
        self.inputs = TypeCheckedList(Variable, [])
        self.outputs = TypeCheckedList(Variable, [])
        
        # # Store all the calls to the black box in case you need them, stores as an AnalysisData object,
        # # but must be included in the manual writing of the self.BlackBox function
        # # Disabled in corsairlite to remove the pandas dependency
        # self.data = AnalysisData()
        
        # A simple description of the model
        self.description = None

        # Defines the order of derivative available in the black box
        self.availableDerivative = 0

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    # These models must be defined in each individual model, just placeholders here
    def BlackBox(*args, **kwargs):
        raise AttributeError(errorString)

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def Gradient_FD(self, x_star, delta_x = 1e-5):
        if self.data.availableDerivative >= 1:
            rs = self.BlackBox(*x_star)
            return rs[1]
        else:
            # Finite difference scheme that is the default if no derivative is available
            # Note: this only works with scalar values in x_star
            inputNames  = [self.inputs[i].name for i in range(0,len(self.inputs))]
            outputNames = [self.outputs[i].name for i in range(0,len(self.outputs))]

            x_star = self.sanitizeInputs(*x_star)
            opt_xstar = self.BlackBox(*x_star)
            if len(self.outputs) == 1:
                opt_xstar = [opt_xstar]
            opt_xstar = self.checkOutputs(**dict(zip(outputNames, opt_xstar)))
            if len(self.outputs) == 1:
                opt_xstar = [opt_xstar]

            dataDict = dict(zip(inputNames, x_star))
            dataDict.update(dict(zip(outputNames,opt_xstar)))
            self.data.addData([dataDict])

            derivativeVectors = []
            for ii in range(0,len(self.outputs)):
                dVec = []
                for i in range(0,len(self.inputs)):
                    x_eval = copy.deepcopy(x_star)
                    x_eval[i] *= (1+delta_x)
                    opt_xeval = self.BlackBox(*x_eval)
                    if len(self.outputs) == 1:
                        opt_xeval = [opt_xeval]
                    opt_xeval = self.checkOutputs(**dict(zip(outputNames, opt_xeval)))
                    if len(self.outputs) == 1:
                        opt_xeval = [opt_xeval]
                    dataDict = dict(zip(inputNames, x_eval))
                    dataDict.update(dict(zip(outputNames,opt_xeval)))
                    self.data.addData([dataDict])
                    
                    d_dx = (opt_xeval[ii] - opt_xstar[ii])/(delta_x * x_star[i])
                    dVec.append(d_dx)
                derivativeVectors.append(dVec)

            return derivativeVectors


# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def parseInputs(*args, **kwargs):
        args = list(args)       # convert tuple to list
        self = args.pop(0)      # pop off the self argument

        inputNames = [self.inputs[i].name for i in range(0,len(self.inputs))]
        # outputNames = [self.outputs[i].name for i in range(0,len(self.outputs))]
        # derivativeNames = [self.derivativePrefix + '%s'%(i) for i in range(1,self.availableDerivative+1)]
        # nameList = inputNames + outputNames + derivativeNames

        # ------------------------------
        # ------------------------------
        if len(args) + len(kwargs.values()) == 1:
            if len(args) == 1:
                inputData = args[0]
            if len(kwargs.values()) == 1:
                inputData = list(kwargs.values())[0]

            if len(inputNames) == 1:
                try:
                    rs = self.sanitizeInputs(inputData)
                    return [dict(zip(inputNames,[rs]))], -self.availableDerivative-1, {} # one input being passed in
                except:
                    pass #otherwise, proceed
        
            if isinstance(inputData, (list,tuple)):
                dataRuns = []
                for idc in inputData:
                    if isinstance(idc, dict):
                        sips = self.sanitizeInputs(**idc)
                        if len(inputNames) == 1:
                            sips = [sips]
                        runDictS = dict(zip(inputNames,sips))
                        dataRuns.append(runDictS)                # the BlackBox([{'x1':x1, 'x2':x2},{'x1':x1, 'x2':x2},...]) case
                    elif isinstance(idc,(list,tuple)):
                        if len(idc) == len(inputNames):
                            sips = self.sanitizeInputs(*idc)
                            if len(inputNames) == 1:
                                sips = [sips]
                            runDictS = dict(zip(inputNames,sips))
                            dataRuns.append(runDictS)                # the BlackBox([ [x1, x2], [x1, x2],...]) case
                        else:
                            raise ValueError('Entry in input data list has improper length')
                    else:
                        raise ValueError("Invalid data type in the input list.  Note that BlackBox([x1,x2,y]) must be passed in as BlackBox([[x1,x2,y]]) or "+
                                                "BlackBox(*[x1,x2,y]) or BlackBox({'x1':x1,'x2':x2,'y':y}) or simply BlackBox(x1, x2, y) to avoid processing singularities.  "+ 
                                                "Best practice is BlackBox({'x1':x1,'x2':x2,'y':y})")
                return dataRuns, self.availableDerivative, {}

            elif isinstance(inputData, dict):
                if set(list(inputData.keys())) == set(inputNames):
                    try:
                        inputLengths = [len(inputData[kw]) for kw in inputData.keys()]
                    except:
                        sips = self.sanitizeInputs(**inputData)
                        if len(inputNames) == 1:
                            sips = [sips]
                        return [dict(zip(inputNames,sips))], -self.availableDerivative-1, {} # the BlackBox(*{'x1':x1, 'x2':x2}) case, somewhat likely

                    if not all([inputLengths[i] == inputLengths[0] for i in range(0,len(inputLengths))]):
                        sips = self.sanitizeInputs(**inputData)
                        return [dict(zip(inputNames,sips))], -self.availableDerivative-1, {} # the BlackBox(*{'x1':x1, 'x2':x2}) case where vectors for x1... are passed in (likely to fail on previous line)
                    else:
                        try:
                            sips = self.sanitizeInputs(**inputData)
                            return [dict(zip(inputNames,sips))], -self.availableDerivative-1, {} # the BlackBox(*{'x1':x1, 'x2':x2}) case where vectors all inputs have same length intentionally (likely to fail on previous line)
                        except:
                            dataRuns = []
                            for i in range(0,inputLengths[0]):
                                runDict = {}
                                for ky,vl in inputData.items():
                                    runDict[ky] = vl[i]
                                sips = self.sanitizeInputs(**runDict)
                                if len(inputNames) == 1:
                                    sips = [sips]
                                runDictS = dict(zip(inputNames,sips))
                                dataRuns.append(runDictS)
                            return dataRuns, self.availableDerivative, {} # the BlackBox({'x1':x1_vec, 'x2':x2_vec}) case, most likely
                else:
                    raise ValueError('Keywords did not match the exptected list')
            else:
                raise ValueError('Got unexpected data type %s'%(str(type(inputData))))
        # ------------------------------
        # ------------------------------
        else:
            if any([ list(kwargs.keys())[i] in inputNames for i in range(0,len(list(kwargs.keys()))) ]):
                # some of the inputs are defined in the kwargs
                if len(args) >= len(inputNames):
                    raise ValueError('A keyword input is defining an input, but there are too many unkeyed arguments for this to occour.  Check the inputs.')
                else:
                    if len(args) != 0:
                        availableKeywords = inputNames[-len(args):]
                    else: 
                        availableKeywords = inputNames

                    valList = args + [None]*(len(inputNames)-len(args))
                    for ky in availableKeywords:
                        ix = inputNames.index(ky)
                        valList[ix] = kwargs[ky]

                    if any([valList[i]==None for i in range(0,len(valList))]):
                        raise ValueError('Kewords did not properly fill in the remaining arguments. Check the inputs.')

                    sips = self.sanitizeInputs(*valList)
                    if len(inputNames) == 1:
                        sips = [sips]

                    remainingKwargs = copy.deepcopy(kwargs)
                    for nm in inputNames:
                        del remainingKwargs[nm]

                    return [dict(zip(inputNames,sips))], -self.availableDerivative-1, remainingKwargs # Mix of args and kwargs define inputs
            else:
                # all of the inputs are in the args
                try:
                    sips = self.sanitizeInputs(*args[0:len(inputNames)])
                    if len(inputNames) == 1:
                        sips = [sips]

                    remainingKwargs = copy.deepcopy(kwargs)
                    remainingKwargs['remainingArgs'] = args[len(inputNames):]
                    return [dict(zip(inputNames,sips))], -self.availableDerivative-1, remainingKwargs # all inputs are in args
                except:
                    runCases, returnMode, extra_singleInput = self.parseInputs(args[0])
                    remainingKwargs = copy.deepcopy(kwargs)
                    remainingKwargs['remainingArgs'] = args[len(inputNames):]
                    return runCases, returnMode, remainingKwargs

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def sanitizeInputs(self, *args, **kwargs):
        nameList = [self.inputs[i].name for i in range(0,len(self.inputs))]
        if len(args) + len(kwargs.values()) > len(nameList):
            raise ValueError('Too many inputs')
        if len(args) + len(kwargs.values()) < len(nameList):
            raise ValueError('Not enough inputs')
        inputDict = {}
        for i in range(0,len(args)):
            rg = args[i]
            inputDict[nameList[i]] = rg

        for ky, vl in kwargs.items():
            if ky in nameList:
                inputDict[ky] = vl
            else:
                raise ValueError('Unexpected input keyword argument %s in the inputs'%(ky))

        opts = []

        for i in range(0,len(nameList)):
            name = nameList[i]
            nameCheck = self.inputs[i].name
            unts = self.inputs[i].units
            size = self.inputs[i].size
            
            if name != nameCheck:
                raise RuntimeError('Something went wrong and values are not consistent.  Check your defined inputs.')

            ipval = inputDict[name]

            if isinstance(ipval, Variable):
                ipval = ipval.value

            if unts is not None:
                if not isinstance(ipval, Q_):
                    ipval = ipval * units.dimensionless
                try:
                    ipval_correctUnits = ipval.to(unts)
                except:
                    raise ValueError('Could not convert %s of %s to %s'%(name, str(ipval),str(unts)))
            else:
                ipval_correctUnits = ipval

            if size is not None:
                try:
                    szVal = ipval_correctUnits.magnitude
                except:
                    szVal = ipval_correctUnits
                if size == 0:
                    if not isinstance(szVal, (int, float)):
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))
                elif isinstance(size, int):
                    ck1 = len(szVal) != size
                    ck2 = len(np.array(szVal).shape) != 1
                    ck3 = np.array(szVal).shape[0] != size
                    if ck1 or ck2 or ck3:
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))
                else:
                    ck1 = np.array(szVal).shape != tuple(size)
                    if ck1:
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))

            opts.append(ipval_correctUnits)
        if len(opts) == 1:
            opts = opts[0]

        return opts

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def checkOutputs(self, *args, **kwargs):
        nameList = [self.outputs[i].name for i in range(0,len(self.outputs))]
        if len(args) + len(kwargs.values()) > len(nameList):
            raise ValueError('Too many inputs')
        if len(args) + len(kwargs.values()) < len(nameList):
            raise ValueError('Not enough inputs')
        inputDict = {}
        for i in range(0,len(args)):
            rg = args[i]
            inputDict[nameList[i]] = rg

        for ky, vl in kwargs.items():
            if ky in nameList:
                inputDict[ky] = vl
            else:
                raise ValueError('Unexpected input keyword argument %s in the inputs'%(ky))

        opts = []

        for i in range(0,len(nameList)):
            name = nameList[i]
            nameCheck = self.outputs[i].name
            unts = self.outputs[i].units
            size = self.outputs[i].size
            
            if name != nameCheck:
                raise RuntimeError('Something went wrong and values are not consistent.  Check your defined inputs.')

            ipval = inputDict[name]

            if isinstance(ipval, Variable):
                ipval = ipval.value

            if unts is not None:
                if isinstance(ipval, (int,float)):
                    ipval = ipval * units.dimensionless
                try:
                    ipval_correctUnits = ipval.to(unts)
                except:
                    raise ValueError('Could not convert %s of %s to %s'%(name, str(ipval),str(unts)))
            else:
                ipval_correctUnits = ipval

            if size is not None:
                try:
                    szVal = ipval_correctUnits.magnitude
                except:
                    szVal = ipval_correctUnits
                if size == 0:
                    if not isinstance(szVal, (int, float)):
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))
                elif isinstance(size, int):
                    ck1 = len(szVal) != size
                    ck2 = len(np.array(szVal).shape) != 1
                    ck3 = np.array(szVal).shape[0] != size
                    if ck1 or ck2 or ck3:
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))
                else:
                    ck1 = np.array(szVal).shape != tuple(size)
                    if ck1:
                        raise ValueError('Size of %s did not match the expected size %s'%(name, str(size)))

            opts.append(ipval_correctUnits)
        if len(opts) == 1:
            opts = opts[0]

        return opts

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def getSummary(self, whitespace = 6):
        pstr = ''
        pstr += 'Model Description\n'
        pstr += '=================\n'
        descr_str = self.description.__repr__()
        pstr += descr_str[1:-1] + '\n\n'
        
        longestName = 0
        longestUnits = 0
        longestSize = 0
        for ipt in self.inputs:
            nml = len(ipt.name)
            if nml > longestName:
                longestName = nml
            if ipt.units is None:
                unts = 'None'
            else:
                unts = ipt.units._repr_html_()
                unts = unts.replace('<sup>','^')
                unts = unts.replace('</sup>','')
                unts = unts.replace('\[', '[')
                unts = unts.replace('\]', ']')
            unl = len(unts)
            if unl > longestUnits:
                longestUnits = unl
            if ipt.size is None:
                lsz = 4
            else:
                if type(ipt.size) == list:
                    lsz = len(ipt.size.__repr__())
                else:
                    lsz = len(str(ipt.size))
            if lsz > longestSize:
                longestSize = lsz
        namespace = max([4, longestName]) + whitespace
        unitspace = max([5, longestUnits]) + whitespace
        sizespace = max([4, longestSize]) + whitespace
        fulllength = namespace+unitspace+sizespace+11
        pstr += 'Inputs\n'
        pstr += '='*fulllength
        pstr += '\n'
        pstr += 'Name'.ljust(namespace)
        pstr += 'Units'.ljust(unitspace)
        pstr += 'Size'.ljust(sizespace)
        pstr += 'Description'
        pstr += '\n'
        pstr += '-'*(namespace - whitespace)
        pstr += ' '*whitespace
        pstr += '-'*(unitspace - whitespace)
        pstr += ' '*whitespace
        pstr += '-'*(sizespace - whitespace)
        pstr += ' '*whitespace
        pstr += '-----------'
        pstr += '\n'
        for ipt in self.inputs:
            pstr += ipt.name.ljust(namespace)
            if ipt.units is None:
                unts = 'None'
            else:
                unts = ipt.units._repr_html_()
                unts = unts.replace('<sup>','^')
                unts = unts.replace('</sup>','')
                unts = unts.replace('\[', '[')
                unts = unts.replace('\]', ']')
            pstr += unts.ljust(unitspace)
            if ipt.size is None:
                lnstr = 'None'
            else:
                lnstr = '%s'%(ipt.size.__repr__())
            pstr += lnstr.ljust(sizespace)
            pstr += ipt.description
            pstr += '\n'
        pstr += '\n'
        
        longestName = 0
        longestUnits = 0
        longestSize = 0
        for opt in self.outputs:
            nml = len(opt.name)
            if nml > longestName:
                longestName = nml
            if opt.units is None:
                unts = 'None'
            else:
                unts = opt.units._repr_html_()
                unts = unts.replace('<sup>','^')
                unts = unts.replace('</sup>','')
                unts = unts.replace('\[', '[')
                unts = unts.replace('\]', ']')
            unl = len(unts)
            if unl > longestUnits:
                longestUnits = unl
            if opt.size is None:
                lsz = 4
            else:
                if type(opt.size) == list:
                    lsz = len(opt.size.__repr__())
                else:
                    lsz = len(str(opt.size))
            if lsz > longestSize:
                longestSize = lsz
        namespace = max([4, longestName]) + whitespace
        unitspace = max([5, longestUnits]) + whitespace
        sizespace = max([4, longestSize]) + whitespace
        fulllength = namespace+unitspace+sizespace+11
        pstr += 'Outputs\n'
        pstr += '='*fulllength
        pstr += '\n'
        pstr += 'Name'.ljust(namespace)
        pstr += 'Units'.ljust(unitspace)
        pstr += 'Size'.ljust(sizespace)
        pstr += 'Description'
        pstr += '\n'
        pstr += '-'*(namespace - whitespace)
        pstr += ' '*whitespace
        pstr += '-'*(unitspace - whitespace)
        pstr += ' '*whitespace
        pstr += '-'*(sizespace - whitespace)
        pstr += ' '*whitespace
        pstr += '-----------'
        pstr += '\n'
        for opt in self.outputs:
            pstr += opt.name.ljust(namespace)
            if opt.units is None:
                unts = 'None'
            else:
                unts = opt.units._repr_html_()
                unts = unts.replace('<sup>','^')
                unts = unts.replace('</sup>','')
                unts = unts.replace('\[', '[')
                unts = unts.replace('\]', ']')
            pstr += unts.ljust(unitspace)
            if opt.size is None:
                lnstr = 'None'
            else:
                lnstr = '%s'%(opt.size.__repr__())
            pstr += lnstr.ljust(sizespace)
            pstr += opt.description
            pstr += '\n'
        pstr += '\n'
        
        return pstr

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    
    @property
    def summary(self):
        return self.getSummary()

# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def __repr__(self):
        pstr = 'AnalysisModel( ['
        for i in range(0,len(self.outputs)):
            pstr += self.outputs[i].name
            pstr += ','
        pstr = pstr[0:-1]
        pstr += ']'
        pstr += ' == '
        pstr += 'f('
        for ipt in self.inputs:
            pstr += ipt.name
            pstr += ', '
        pstr = pstr[0:-2]
        pstr += ')'
        pstr += ' , '
        pstr = pstr[0:-2]
        pstr += '])'
        return pstr
# ---------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------
    def setupDataContainer(self):
        self.data.inputs = copy.deepcopy(self.inputs)
        self.data.outputs = copy.deepcopy(self.outputs)
        self.data.availableDerivative = self.availableDerivative

# import os
# import copy
# import numpy as np
# from corsairlite import units
# from corsairlite.core.dataTypes.variable import Variable
# from corsairlite.analysis.models.model import Model
    
# class SignomialTest(Model):
#     def __init__(self):
#         # Set up all the attributes by calling Model.__init__
#         super(SignomialTest, self).__init__()
        
#         #Setup Inputs
#         self.inputs.append(Variable('x', None, units.dimensionless,  'Independent Variable'))
        
#         #Setup Outputs
#         self.outputs.append(Variable('y', None, units.dimensionless,  'Dependent Variable'))

#         #Set the data storage inputs and outputs
#         self.data.inputs = copy.deepcopy(self.inputs)
#         self.data.outputs = copy.deepcopy(self.outputs)
    
#         #Simple model description
#         self.description = 'This model evaluates the function: max([-6*x-6, x**4-3*x**2])'
        
#         self.availableDerivative = 1
#         self.data.availableDerivative = self.availableDerivative
    
#     #standard function call is y(, dydx, ...) = self.BlackBox({'x1':x1, 'x2':x2})
#     def BlackBox(*args, **kwargs):
#         args = list(args)       # convert tuple to list
#         self = args.pop(0)      # pop off the self argument
#         runCases, returnMode, extra = self.parseInputs(*args, **kwargs)
#         print(extra)
#         x = [ runCases[i]['x'] for i in range(0,len(runCases)) ]
#         x = np.array([ self.sanitizeInputs(xval).to('').magnitude for xval in x ])
#         y = np.maximum(-6*x-6, x**4-3*x**2)
        
#         dydx = 4*x**3 - 6*x
#         ddy_ddx = 12*x**2 - 6
#         gradientSwitch = -6*x-6 > x**4-3*x**2
#         dydx[gradientSwitch] = -6
#         ddy_ddx[gradientSwitch] = 0

#         self.data.addData({'x': x.tolist(), 'y': y.tolist(), '__D1':[[[dydx[i]]] for i in range(0,len(dydx))]})

# #         for i in range(0,len(runCases)):
# #             self.data.addData({'x': x[i], 'y': y[i], '__D1':[[dydx[i]]]})

# #         for i in range(0,len(runCases)):
# #             self.data.addData(*[x[i], y[i], [[dydx[i]]]])

# #         addDataMatrix = []
# #         for i in range(0,len(runCases)):
# #             addDataMatrix.append({'x': x[i], 'y': y[i], '__D1':[[dydx[i]]]})
# #         self.data.addData(addDataMatrix)
        
# #         addDataMatrix = []
# #         for i in range(0,len(runCases)):
# #             addDataMatrix.append([ x[i],  y[i],[[dydx[i]]]])
# #         self.data.addData(addDataMatrix)
    
#         y = [ self.checkOutputs(yval) for yval in y ]
#         dydx = [dydx[i] * units.dimensionless for i in range(0,len(dydx))]
        
#         # etc...
#         # Mode -1: y with no case wrapper, single point run
#         # Mode 0: y only
#         # Mode 1: y dydx
#         # etc...
        
#         if returnMode < 0:
#             returnMode = -1*(returnMode + 1)
#             if returnMode == 0:
#                 return y[0]
#             if returnMode == 1:
#                 return y[0], dydx
#         else:
#             if returnMode == 0:
#                 opt = []
#                 for i in range(0,len(y)):
#                     opt.append([ y[i] ])
#                 return opt
#             if returnMode == 1:
#                 opt = []
#                 for i in range(0,len(y)):
#                     opt.append([ [y[i]], [ [[dydx[i]]] ] ])
#                 return opt

# s = SignomialTest()
# ivals = [[x] for x in np.linspace(-2,2,11)]

# # bbo = s.BlackBox(**{'x':0.5})
# # bbo = s.BlackBox({'x':0.5})
# # bbo = s.BlackBox(**{'x':0.5, 'optn':True})
# # bbo = s.BlackBox(*[0.5], **{'optn1': True, 'optn2': False})
# # bbo = s.BlackBox(*[0.5,True], **{'optn': False})
# # bbo = s.BlackBox({'x':[x for x in np.linspace(-2,2,11)]})
# # bbo = s.BlackBox([{'x':x} for x in np.linspace(-2,2,11)])
# # bbo = s.BlackBox([ [x] for x in np.linspace(-2,2,11)])
# # bbo = s.BlackBox([ [x] for x in np.linspace(-2,2,11)], True, optn=False)
# # bbo = s.BlackBox([ [x] for x in np.linspace(-2,2,11)], optn1=True, optn2=False)
# print(bbo)
# # y = [bbo[i][0][0] for i in range(0,len(bbo))] #runcase, 0th order output, first output value
# # dydx = [bbo[i][1][0][0][0] for i in range(0,len(bbo))] #runcase, 1st order output, first output value, dx1, dx1
# # print(y)
# # print(dydx)
# # print(s)
# # print(s.data.dataAll)
# # s.Gradient_FD({'x':0.5})