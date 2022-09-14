import copy
import random
import numpy as np
import multiprocessing
import math
import dill
import corsairlite
from corsairlite import units, Q_
from corsairlite.core.dataTypes import OptimizationOutput, TypeCheckedList
from corsairlite.optimization import Variable, Constant, Term, RuntimeConstraint

# Required Options:
# options.solver
# options.solveType
# options.returnRaw

# =====================================================================================================================
# Pickles variable
# =====================================================================================================================
def variable__getstate__(self):
    dct = {}
    dct['name'] = self.name
    dct['guess'] = self.guess.magnitude
    dct['guess_units'] = str(self.guess.units)
    dct['units'] = str(self.units)
    dct['description'] = self.description
    dct['type_'] = self.type_
    if self.bounds is not None:
        dct['bounds'] = [self.bounds[0].magnitude, self.bounds[1].magnitude]
        dct['bounds_units_lower'] = str(self.bounds[0].units)
        dct['bounds_units_upper'] = str(self.bounds[1].units)
    else:
        dct['bounds'] = None
    dct['size'] = self.size
    dct['transformation'] = self.transformation 
    return dct
# =====================================================================================================================
# Unpickles variable
# =====================================================================================================================
def variable__setstate__(self,dct):
    from corsairlite.core import pint_custom
    units = pint_custom.get_application_registry()
    self.units = units(dct['units']).units
    self.name = dct['name']
    self.size = dct['size']
    self.guess = dct['guess'] * corsairlite.units(dct['guess_units'])
    self.description = dct['description']
    self.type_ = dct['type_']
    if dct['bounds'] is not None:
        self.bounds = [dct['bounds'][0] * corsairlite.units(dct['bounds_units_lower']), dct['bounds'][1] * corsairlite.units(dct['bounds_units_upper'])]
    else:
        self.bounds = None
    self.transformation = dct['transformation']
# =====================================================================================================================
# Pickles constant
# =====================================================================================================================
def constant__getstate__(self):
    dct = {}
    dct['name'] = self.name
    dct['value'] = self.value.magnitude
    dct['value_units'] = str(self.value.units)
    dct['units'] = str(self.units)
    dct['description'] = self.description
    dct['size'] = self.size
    return dct
# =====================================================================================================================
# Unpickles constant
# =====================================================================================================================
def constant__setstate__(self,dct):
    from corsairlite.core import pint_custom
    units = pint_custom.get_application_registry()
    self.units = corsairlite.units(dct['units']).units
    self.name = dct['name']
    self.size = dct['size']
    self.value = dct['value'] * corsairlite.units(dct['value_units'])
    self.description = dct['description']
# =====================================================================================================================
# Pickles term
# =====================================================================================================================
def term__getstate__(self):
    dct = {}
    dct['variables'] = self.variables
    dct['constant'] = self.constant.magnitude
    dct['constant_units'] = str(self.constant.units)
    dct['exponents'] = self.exponents
    hcs = []
    hcus = []
    for hc in self.hiddenConstants:
        if isinstance(hc,(int,float)):
            hcs.append(hc)
            hcus.append(None)
        elif isinstance(hc,Q_):
            hcs.append(hc.magnitude)
            hcus.append(str(hc.units))
        else:
            raise ValueError('Invalid hiddenConstant type')
    dct['hiddenConstants'] = hcs
    dct['hiddenConstants_units'] = hcus
    dct['hiddenConstantExponents'] = self.hiddenConstantExponents
    dct['runtimeExpressions'] = self.runtimeExpressions
    dct['runtimeExpressionExponents'] = self.runtimeExpressionExponents
    return dct
# =====================================================================================================================
# Unpickles term
# =====================================================================================================================
def term__setstate__(self,dct):
    from corsairlite.core import pint_custom
    units = pint_custom.get_application_registry()
    self.variables = dct['variables']
    self.constant = dct['constant'] * corsairlite.units(dct['constant_units'])
    self.exponents = dct['exponents']
    hcl = TypeCheckedList((int,float, Q_), [])
    for i in range(0,len(dct['hiddenConstants'])):
        hc = dct['hiddenConstants'][i]
        un = dct['hiddenConstants_units'][i]
        if un is None:
            hcl.append(hc)
        else:
            hcl.append(hc * corsairlite.units(un))
    self.hiddenConstants = []
    self.hiddenConstantExponents = dct['hiddenConstantExponents']
    self.runtimeExpressions = dct['runtimeExpressions']
    self.runtimeExpressionExponents = dct['runtimeExpressionExponents']

# =====================================================================================================================
# Pickles RuntimeConstraint
# =====================================================================================================================
def rtc__getstate__(self):
    dct = {}
    dct['inputs'] = self.inputs
    dct['outputs'] = self.outputs
    dct['operators'] = self.operators
    dct['analysisModel'] = self.analysisModel
    dct['scaledConstraintVector'] = self.scaledConstraintVector
    dct['scaledVariablesVector'] = self.scaledVariablesVector
    return dct
# =====================================================================================================================
# Unpickles RuntimeConstraint
# =====================================================================================================================
def rtc__setstate__(self,dct):
     self.inputs = dct['inputs'] 
     self.outputs = dct['outputs'] 
     self.operators = dct['operators'] 
     self.analysisModel = dct['analysisModel'] 
     self.scaledConstraintVector = dct['scaledConstraintVector'] 
     self.scaledVariablesVector = dct['scaledVariablesVector'] 

# =====================================================================================================================
# The pickle functions
# =====================================================================================================================
def setPickle():
    Variable.__getstate__ = variable__getstate__
    Variable.__setstate__ = variable__setstate__
    Constant.__getstate__ = constant__getstate__
    Constant.__setstate__ = constant__setstate__
    Term.__getstate__ = term__getstate__
    Term.__setstate__ = term__setstate__
    RuntimeConstraint.__getstate__ = rtc__getstate__
    RuntimeConstraint.__setstate__ = rtc__setstate__

def unsetPickle():
    delattr(Variable,'__getstate__')
    delattr(Variable,'__setstate__')
    delattr(Constant,'__getstate__')
    delattr(Constant,'__setstate__')
    delattr(Term,'__getstate__')
    delattr(Term,'__setstate__')
    delattr(RuntimeConstraint,'__getstate__')
    delattr(RuntimeConstraint,'__setstate__')