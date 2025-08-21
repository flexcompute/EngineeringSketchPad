
from .caps import *

from .legacy import capsConvert
# from .legacy import capsProblem, capsAnalysis, capsGeometry, capsBound
# from .legacy import CapsProblem, CapsAnalysis, CapsGeometry, CapsBound

from .problem import Problem
from .problem import capsProblem, capsAnalysis, capsGeometry, capsBound
from .problem import CapsProblem, CapsAnalysis, CapsGeometry, CapsBound

def setLegacy(LEGACY):
    from . import problem
    problem.LEGACY = LEGACY