from corsairlite.optimization import Variable
# from corsairlite.analysis.model import Model
# from corsairlite.analysis.fits.runtimeFit import RuntimeFit
# from corsairlite.analysis.modelVariableDescription import ModelVariableDescription 
from corsairlite.core.dataTypes import TypeCheckedList

class RuntimeExpression(object):
    def __init__(self, analysisModel = None, inputs = None, output = None, fittingOptions= None):
        # self.analysisModel = analysisModel
        if inputs is None:
            self.inputs = TypeCheckedList(Variable, [])
        else:
            self.inputs = inputs
        self.output = output
        # self.fittingOptions = fittingOptions
    
    # @property
    # def analysisModel(self):
    #     return self._analysisModel
    # @analysisModel.setter
    # def analysisModel(self,val):
    #     if val is None:
    #         self._analysisModel = None
    #     else:
    #         if isinstance(val, Model):
    #             self._analysisModel = val
    #         else:
    #             raise ValueError('Invalid analysisModel.  Must be an instance of corsairlite.analysis.model.Model')

    @property
    def inputs(self):
        return self._inputs
    @inputs.setter
    def inputs(self, val):
        if isinstance((list, tuple), val):
            if all(isinstance(x, Variable) for x in val):
                self._inputs = val
            else:
                raise ValueError('Input list must be a list or tuple.  TypeCheckedList and list types are preferred.')
        else:
            raise ValueError('Input list must be a list or tuple.  TypeCheckedList and list types are preferred.')

    @property
    def output(self):
        return self._output
    @output.setter
    def output(self,val):
        if val is None:
            self._output = None
        else:
            if isinstance(val, Variable):
                self._output = val
            else:
                raise ValueError('Invalid output.  Must be a Variable.')

    # @property
    # def fittingFunction(self):
    #     return self._fittingFunction
    # @fittingFunction.setter
    # def fittingFunction(self, val):
    #     if val is None:
    #         self._fittingFunction = None
    #     else:
    #         if isinstance(val, RuntimeFit):
    #             self._fittingFunction = val
    #         else:
    #             raise ValueError('Invalid fittingFunction.  Must be a fitting function from the corsairlite.analysis.fits.runtime library')


    def __repr__(self):
        pstr = ''
        pstr += self.output.name
        pstr += '('
        for ipt in self.inputs:
            pstr += ipt.name
            pstr += ', '
        pstr = pstr[0:-2]
        pstr += ')'
        return pstr