###################################################################
#                                                                 #
# pyESP --- Python interface to serveESP/timPython                #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# Copyright (C) 2021  John F. Dannenhoffer, III (Syracuse University)
#
# This library is free software; you can redistribute it and/or
#    modify it under the terms of the GNU Lesser General Public
#    License as published by the Free Software Foundation; either
#    version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
#    License along with this library; if not, write to the Free Software
#    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#     MA  02110-1301  USA

import ctypes
import os
import sys
import atexit

from   pyEGADS import egads

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

# load the shared library
if sys.platform.startswith('darwin'):
    _esp = ctypes.CDLL(_ESP_ROOT + "/lib/python.so")
elif sys.platform.startswith('linux'):
    _esp = ctypes.CDLL(_ESP_ROOT + "/lib/python.so")
elif sys.platform.startswith('win32'):
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        _esp = ctypes.CDLL(_ESP_ROOT + "\\lib\\python.dll")
    else:
        _esp = ctypes.CDLL(_ESP_ROOT + "\\lib\\python.dll", winmode=0)
else:
    raise IOError("Unknown platform: " + sys.platform)

# ======================================================================

def GetModl():
    """
    GetModl - get serveESP's active MODL

    inputs:
        (None)
    outputs:
        (None)
    """
    _esp.timGetModl.argtypes = [ctypes.POINTER(ctypes.c_void_p)]
    _esp.timGetModl.restype  =  ctypes.c_int

    modl   = ctypes.c_void_p()
    status =_esp.timGetModl(modl)
    _processStatus(status, "GetModl")
    
    return modl

# ======================================================================

def SetModl(modl):
    """
    SetModl - set serveESP's active MODL

    inputs:
        modl     OpenCSM MODL
    outputs:
        (None}
    """
    _esp.timSetModl.argtypes = [ctypes.c_void_p]
    _esp.timSetModl.restype  =  ctypes.c_int

    status = _esp.timSetModl(modl._modl)
    _processStatus(status, "SetModl")

    # make sure modl is not cleaned up when python exits
    if modl._finalize:
        modl._finalize.detach()
    
    modl._pyowned  = False
    modl._finalize = None

    # since we are passing the MODL back to serveESP, take ownership
    #     of the EGADS context
    if modl._context:
        modl._context.py_to_c(takeOwnership=True)

    return None

# ======================================================================

def ViewModl(modl):
    """
    ViewModl - make modl active and view it in serveESP

    inputs:
        modl     OpenCSM MODL
    outputs:
        (None}
    """
    _esp.timViewModl.argtypes = [ctypes.c_void_p]
    _esp.timViewModl.restype  =  ctypes.c_int

    status = _esp.timViewModl(modl._modl)
    _processStatus(status, "ViewModl")

    return None

# ======================================================================

# ======================================================================
# Helper routines
# ======================================================================

def _processStatus(status, routine):
    """
    _processStatus - raise error when status < 0

    inputs:
        status      return status
        routine     routine where error occurred
    outputs:
        (None)
    """

    if   (status < 0):
        raise EspError("esp error", routine)

# ======================================================================
# EspError class
# ======================================================================

class EspError(Exception):

    # Constructor or Initializer
    def __init__(self, value, routine):
        self.value   = value
        self.routine = routine

    # __str__ is to print() the value
    def __str__(self):
        return("EspError: "+repr(self.value)+" detected during call to "+self.routine)

# ======================================================================

