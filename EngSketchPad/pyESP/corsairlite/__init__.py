##  Order is important in this file!!!!

from .core.pint_custom import UnitRegistry
from .core import pint_custom
units = UnitRegistry()
pint_custom.set_application_registry(units)
Q_ = units.Quantity
from .core.pint_custom.unit import _Unit

import functools, re

# Add parsing of UDUNITS-style power
units.preprocessors.append(
    functools.partial(
        re.sub,
        r"(?<=[A-Za-z])(?![A-Za-z])(?<![0-9\-][eE])(?<![0-9\-])(?=[0-9\-])",
        "**",
    )
)

import os
pathHere = os.path.realpath(__file__)  
path_to_corsairlite = pathHere[0:-21] # removes the trailing '/corsairlite/__init__.py' giving the package location
if path_to_corsairlite[-1] != 'r':
    path_to_corsairlite = pathHere[0:-20] #need to account for running off of .py and not .pyc the first time

# import shutil
# import warnings
# if shutil.which('avl') is not None:
#     path_to_AVL = shutil.which('avl')
# else:
#     warnings.warn('Could not find an AVL executable')
# if shutil.which('xfoil') is not None:
#     path_to_XFOIL = shutil.which('xfoil')
# else:
#     warnings.warn('Could not find an XFOIL executable')
# if shutil.which('serveCSM') is not None:
#     path_to_ESP_serveCSM = shutil.which('serveCSM')
# else:
#     warnings.warn('Could not find the serveCSM executable')


# import core
# import configuration
# import analysis
# import optimization

# from IPython.display import display
# from corsairlite.core.locationChecks import inNotebook
# from corsairlite.core.locationChecksimport inIPython
# def help(helpItem):
# 	if inIPython():
# 		print helpItem.__doc__
# 	else:
# 		return display(helpItem.__doc__)
