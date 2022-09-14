import copy
from corsairlite import units
from corsairlite.optimization import Variable
from corsairlite.analysis.models import Model
from corsairlite.core.dataTypes import TypeCheckedList
from corsairlite.optimization import Constraint

class RuntimeConstraint(object):
    def __init__(self, outputs = None, operators = None, inputs = None, analysisModel = None):        
        if inputs is None:
            self.inputs = TypeCheckedList(Variable, [])
        else:
            self.inputs = inputs

        if outputs is None:
            self.outputs = TypeCheckedList(Variable, [])
        else:
            self.outputs = outputs

        self.operators = operators
        # self.analysisModel = copy.deepcopy(analysisModel)
        self.analysisModel = analysisModel


        self.scaledConstraintVector = None
        self.scaledVariablesVector  = None
    
    @property
    def analysisModel(self):
        return self._analysisModel
    @analysisModel.setter
    def analysisModel(self,val):
        if val is None:
            self._analysisModel = None
        else:
            if isinstance(val, Model):
                self._analysisModel = val
            else:
                raise ValueError('Invalid analysisModel.  Must be an instance of corsairlite.analysis.model.Model')

    @property
    def inputs(self):
        return self._inputs
    @inputs.setter
    def inputs(self, val):
        if isinstance(val, (list, tuple)):
            if all(isinstance(x, Variable) for x in val):
                self._inputs = val
            else:
                raise ValueError('Input list must be a list or tuple.  TypeCheckedList and list types are preferred.')
        else:
            raise ValueError('Input list must be a list or tuple.  TypeCheckedList and list types are preferred.')

    @property
    def operators(self):
        return self._operators
    @operators.setter
    def operators(self, val):
        if isinstance(val, (list, tuple)):
            if all(x in ['==','<=','>='] for x in val):
                if len(val) == len(self.outputs):
                    self._operators = val
                else:
                    raise ValueError('Operator list is not of correct length')
            else:
                raise ValueError('Operator list has an invalid value')
        else:
            raise ValueError('Operator list must be a list or tuple.')
            
    @property
    def outputs(self):
        return self._outputs
    @outputs.setter
    def outputs(self, val):
        if isinstance(val, (list, tuple)):
            if all(isinstance(x, Variable) for x in val):
                self._outputs = val
            else:
                raise ValueError('Output list must be a list or tuple.  TypeCheckedList and list types are preferred.')
        else:
            raise ValueError('Output list must be a list or tuple.  TypeCheckedList and list types are preferred.')

    def returnConstraintSet(self, xStar, fitMode, fitType=None):
        inputNames = [self.inputs[i].name for i in range(0,len(self.inputs))]
        outputNames = [self.outputs[i].name for i in range(0,len(self.outputs))]
        xvars = self.inputs
        yvars = self.outputs
        
        if len(xvars) != len(self.inputs):
            raise ValueError('Input length does not match')            
        if len(yvars) != len(self.outputs):
            raise ValueError('Output length does not match')
        # -------------------------------------------------
        if fitMode == 'raw':
            if self.analysisModel.availableDerivative >= 1:
                rs = self.analysisModel.BlackBox(*xStar)
                yStar = rs[0]
                dyStar_dxStar = rs[1]
                if len(self.outputs) == 1:
                    yStar = [yStar]
                    dyStar_dxStar = [dyStar_dxStar]
            else:
                yStar = self.analysisModel.BlackBox(*xStar)
                dyStar_dxStar = self.analysisModel.Gradient_FD(xStar)
                if len(self.outputs) == 1:
                    yStar = [yStar]
            return yStar, dyStar_dxStar
        # -------------------------------------------------
        elif fitMode.lower() == 'sqp':
            # xStar will be stacked [xStar, yStar]
            inputValues = copy.deepcopy(xStar)
            if self.scaledVariablesVector is not None:
                inputValues = [ inputValues[i] * self.scaledVariablesVector[i] for i in range(0,len(inputValues)) ]

            xStar = inputValues[0:len(inputNames)]
            yValues = inputValues[len(inputNames):]
            if self.analysisModel.availableDerivative >= 1:
                rs = self.analysisModel.BlackBox(*xStar)
                yStar = rs[0]
                dyStar_dxStar = rs[1]
                if len(self.outputs) == 1:
                    yStar = [yStar]
                    dyStar_dxStar = [dyStar_dxStar]
            else:
                yStar = self.analysisModel.BlackBox(*xStar)
                dyStar_dxStar = self.analysisModel.Gradient_FD(xStar)
                if len(self.outputs) == 1:
                    yStar = [yStar]

            sqp_feval = []
            sqp_grad = []
            for i in range(0,len(yValues)):
                sqp_grad_i = []
                if self.operators[i] == '>=':
                    sqp_feval.append(yStar[i] - yValues[i])
                    sqp_grad_i = dyStar_dxStar[i] + [-1*units.dimensionless]
                else:
                    sqp_feval.append(yValues[i] - yStar[i])
                    sqp_grad_i = [-1*dyStar_dxStar[i][ii] for ii in range(0,len(dyStar_dxStar[i]))] + [1*units.dimensionless]
                sqp_grad.append(sqp_grad_i)

            if self.scaledVariablesVector is not None:
                # sqp_feval is untouched
                for ii in range(0,len(sqp_grad)):
                    sqp_grad[ii] = [ sqp_grad[ii][i]*self.scaledVariablesVector[i] for i in range(0,len(sqp_grad[ii]))]                 

            if self.scaledConstraintVector is not None:
                sqp_feval = [ sqp_feval[i]/self.scaledConstraintVector[i] for i in range(0,len(self.outputs))] 
                for ii in range(0,len(sqp_grad)):
                    sqp_grad[ii] = [ sqp_grad[ii][i]/self.scaledConstraintVector[ii] for i in range(0,len(sqp_grad[ii]))] 

            return sqp_feval, sqp_grad
        # -------------------------------------------------
        elif fitMode.lower() == 'lsqp':
            sgpDict = self.returnConstraintSet(xStar,'sgp')
            return sgpDict['con']
    
        elif fitMode.lower() == 'sgp':
            # xStar will be stacked [xStar, yStar]
            inputValues = copy.deepcopy(xStar)
            if self.scaledVariablesVector is not None:
                inputValues = [ inputValues[i] * self.scaledVariablesVector[i] for i in range(0,len(inputValues)) ]

            xStar = inputValues[0:len(inputNames)]
            yValues = inputValues[len(inputNames):]
            if self.analysisModel.availableDerivative >= 1:
                rs = self.analysisModel.BlackBox(*xStar)
                yStar = rs[0]
                dyStar_dxStar = rs[1]
                if len(self.outputs) == 1:
                    yStar = [yStar]
                    dyStar_dxStar = [dyStar_dxStar]
            else:
                yStar = self.analysisModel.BlackBox(*xStar)
                dyStar_dxStar = self.analysisModel.Gradient_FD(xStar)
                if len(self.outputs) == 1:
                    yStar = [yStar]


            sgp_feval = []
            sgp_grad = []
            for i in range(0,len(yValues)):
                sgp_grad_i = []
                if self.operators[i] == '>=':
                    sgp_feval.append(yStar[i]/yValues[i])
                    sgp_grad_i = [ dyStar_dxStar[i][j]/yValues[i] for j in range(0,len(self.inputs)) ] + [-1 * yStar[i]/yValues[i]**2]
                else:
                    sgp_feval.append(yValues[i]/yStar[i])
                    sgp_grad_i = [ -yValues[i]*dyStar_dxStar[i][j]/yStar[i]**2 for j in range(0,len(self.inputs)) ] + [ 1 / yStar[i]]
                sgp_grad.append(sgp_grad_i)

            if self.scaledVariablesVector is not None:
                # sgp_feval is untouched
                for ii in range(0,len(sgp_grad)):
                    sgp_grad[ii] = [ sgp_grad[ii][i]*self.scaledVariablesVector[i] for i in range(0,len(sgp_grad[ii]))]  

            retDict = {}
            retDict['raw'] = (yStar, dyStar_dxStar)
            retDict['con'] = (sgp_feval, sgp_grad)

            return retDict

        # -------------------------------------------------
        elif fitMode == 'local':
        # -------------------------------------------------
            if fitType == 'affine':
                if self.analysisModel.availableDerivative >= 1:
                    rs = self.analysisModel.BlackBox(*xStar)
                    yStar = rs[0]
                    dyStar_dxStar = rs[1]
                    if len(self.outputs) == 1:
                        yStar = [yStar]
                        dyStar_dxStar = [dyStar_dxStar]
                else:
                    yStar = self.analysisModel.BlackBox(*xStar)
                    dyStar_dxStar = self.analysisModel.Gradient_FD(xStar)
                    if len(self.outputs) == 1:
                        yStar = [yStar]

                constraintSet = []
                for i in range(0,len(self.outputs)):
                    conLHS = yvars[i]
                    conOperator = self.operators[i]
                    conRHS = yStar[i]
                    for ii in range(0,len(xvars)):
                        conRHS += dyStar_dxStar[i][ii] * (xvars[i] - xStar[i])
                    con = Constraint(conLHS, conOperator, conRHS)
                    constraintSet.append(con)
                    
                return constraintSet
        # -------------------------------------------------
            elif fitType == 'monomial':
                yStar, dyStar_dxStar = self.returnConstraintSet(xStar,'lsqp')
                constraintSet = []
                for i in range(0,len(self.outputs)):
                    vrs = xvars + [yvars[i]]
                    mon = yStar[i]
                    for ii in range(0,len(vrs)):
                        a_n = xStar[ii]/yStar[i] * dyStar_dxStar[i][ii]
                        mon *= (vrs[ii]/xStar[ii])**a_n.to('').magnitude
                    con = Constraint(mon, '<=', 1*units.dimensionless)
                    constraintSet.append(con)

                return constraintSet
        # -------------------------------------------------
            else:
                raise ValueError('Invalid entry for "fitType", must be one of ["affine","monomial"]')
        # -------------------------------------------------
        elif fitMode == 'fit':
            fitdata = None # need to generate this somehow
            # from corsairlite.fitting.convex.maxAffine import MaxAffineFitObject as MAFO
            # from corsairlite.fitting.convex.softMaxAffine import SoftMaxAffineFitObject as SMAFO
            # from corsairlite.fitting.convex.implicitSoftMaxAffine import ImplicitSoftMaxAffineFitObject as ISMAFO
            # from corsairlite.fitting.nonconvex.differenceOfMaxAffine import DifferenceOfMaxAffineFitObject as DMAFO
            # from corsairlite.fitting.nonconvex.softDifferenceOfMaxAffine import SoftDifferenceOfMaxAffineFitObject as SDMAFO
            if fitType == 'affine':
                from corsairlite.fitting.convex.affine import AffineFitObject as AFO
                fitObject
                print('affine') 
            elif fitType == 'maxaffine':
                print('maxaffine') 
            elif fitType == 'monomial':
                print('monomial') 
            elif fitType == 'monomial_maxaffine':
                print('monomial_maxaffine') 
            elif fitType == 'posynomial_softmaxaffine':
                print('posynomial_softmaxaffine') 
            elif fitType == 'posynomial_implicitsoftmaxaffine':
                print('posynomial_implicitsoftmaxaffine') 
            elif fitType == 'signomial_softdifferenceofmaxaffine':
                print('signomial_softdifferenceofmaxaffine')
            else:
                raise ValueError('Invalid entry for "fitType", must be one of [...]')
                
        else:
            raise ValueError('Invalid mode, must be either "raw", "local", or "fit".')
            
    def __repr__(self):
        pstr = 'RuntimeConstraintSet([ '
        for i in range(0,len(self.outputs)):
            pstr += self.outputs[i].name
            pstr += ' %s '%(self.operators[i])
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

    def substitute(self, substitutions = {}):
        inputNames = [self.inputs[i].name for i in range(0,len(self.inputs))]
        outputNames = [self.outputs[i].name for i in range(0,len(self.outputs))]
        opt = copy.deepcopy(self)

        for i in range(0,len(inputNames)):
            nm = inputNames[i]
            if nm in list(substitutions.keys()):
                newVal = substitutions[nm]
                if isinstance(newVal, Variable):
                    opt.inputs[i] = newVal
                else:
                    raise ValueError('Substitution for %s is not a Variable'%(nm))

        for i in range(0,len(outputNames)):
            nm = outputNames[i]
            if nm in list(substitutions.keys()):
                newVal = substitutions[nm]
                if isinstance(newVal, Variable):
                    opt.outputs[i] = newVal
                else:
                    raise ValueError('Substitution for %s is not a Variable'%(nm))

        return opt




