###################################################################
#                                                                 #
# pyESP --- Python interface to serveESP/timPyscript              #
#                                                                 #
#              Written by John Dannenhoffer @ Syracuse University #
#                                                                 #
###################################################################

# Copyright (C) 2022  John F. Dannenhoffer, III (Syracuse University)
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
from   pyCAPS  import caps

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

# load the shared library
if sys.platform.startswith('darwin'):
    _esp  = ctypes.CDLL(_ESP_ROOT + "/lib/pyscript.so")
    _ocsm = ctypes.CDLL(_ESP_ROOT + "/lib/libocsm.dylib")
elif sys.platform.startswith('linux'):
    _esp  = ctypes.CDLL(_ESP_ROOT + "/lib/pyscript.so")
    _ocsm = ctypes.CDLL(_ESP_ROOT + "/lib/libocsm.so")
elif sys.platform.startswith('win32'):
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        _esp  = ctypes.CDLL(_ESP_ROOT + "\\lib\\pyscript.dll")
        _ocsm = ctypes.CDLL(_ESP_ROOT + "\\lib\\ocsm.dll")
    else:
        _esp  = ctypes.CDLL(_ESP_ROOT + "\\lib\\pyscript.dll", winmode=0)
        _ocsm = ctypes.CDLL(_ESP_ROOT + "\\lib\\ocsm.dll")
else:
    raise IOError("Unknown platform: " + sys.platform)

# ======================================================================

def GetCaps(esp):
    """
    GetCaps - get serveESP's active CAPS

    inputs:
        esp         ESP structure        
    outputs:
        problem     CAPS problem
    """
    _esp.timGetCaps.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                                ctypes.c_void_p]
    _esp.timGetCaps.restype  =  ctypes.c_int

    problem = ctypes.c_void_p()
    status  =_esp.timGetCaps(problem, esp)
    _processStatus(status, "GetCaps")
    
    return caps.capsObj(ctypes.cast(problem,caps.c_capsObj), None, deleteFunction=None)

# ======================================================================

def SetCaps(problem, esp):
    """
    SetCaps - set serveESP's active CAPS

    inputs:
        problem     CAPS problem
        esp         ESP structure
    outputs:
        (None}
    """
    _esp.timSetCaps.argtypes = [ctypes.c_void_p]
    _esp.timSetCaps.restype  =  ctypes.c_int

    status = _esp.timSetCaps(ctypes.cast(problem._obj,caps.c_void_p))
    _processStatus(status, "SetCaps")

    # make sure problem is not cleaned up when python exits
    if problem._finalize:
        problem._finalize.detach()
    
    problem._finalize = None

    return None

# ======================================================================

def GetModl(esp):
    """
    GetModl - get serveESP's active MODL

    inputs:
        esp         ESP structure        
    outputs:
        modl        OpenCSm MODL
    """
    _esp.timGetModl.argtypes = [ctypes.POINTER(ctypes.c_void_p),
                                ctypes.c_void_p]
    _esp.timGetModl.restype  =  ctypes.c_int

    modl   = ctypes.c_void_p()
    status =_esp.timGetModl(modl, esp)
    _processStatus(status, "GetModl")
    
    return modl

# ======================================================================

def SetModl(modl, esp):
    """
    SetModl - set serveESP's active MODL

    inputs:
        modl        OpenCSM MODL
        esp         ESP structure
    outputs:
        (None}
    """
    _esp.timSetModl.argtypes = [ctypes.c_void_p,
                                ctypes.c_void_p]
    _esp.timSetModl.restype  =  ctypes.c_int

    status = _esp.timSetModl(modl._modl, esp)
    _processStatus(status, "SetModl")

    # make sure modl is not cleaned up when python exits
    if modl._pyowned:
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
        modl        OpenCSM MODL
    outputs:
        (None}
    """
    _esp.timViewModl.argtypes = [ctypes.c_void_p]
    _esp.timViewModl.restype  =  ctypes.c_int

    status = _esp.timViewModl(modl._modl)
    _processStatus(status, "ViewModl")

    return None

# ======================================================================

def GetEsp(timName):
    """
    ocsm.GetEsp - get the ESP structure

    inputs:
        (None)
    outputs:
        ESP         ESP structure
    """

    _ocsm.tim_getEsp.argtypes = [ctypes.c_char_p]
    _ocsm.tim_getEsp.restype  =  ctypes.c_void_p

    if (isinstance(timName, str)):
        timName = timName.encode()
            
    ESP = _ocsm.tim_getEsp(timName)

    return ESP
    
# ======================================================================

def TimLoad(timName, esp, data):
    """
    ocsm.TimLoad - load a TIM

    inputs:
        timName     name of TIM
        esp         ESP structure
        data        data used to initialize TIM
    outputs:
        (None)
    """

    _ocsm.tim_load.argtypes = [ctypes.c_char_p,
                               ctypes.c_void_p,
                               ctypes.c_char_p]
    _ocsm.tim_load.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    if (isinstance(data, str)):
        data = data.encode()

    status = _ocsm.tim_load(timName, esp, data)
    _processStatus(status, "TimLoad")

    return

# ======================================================================

def TimQuit(timName):
    """
    ocsm.TimQuit - quits a TIM

    inputs:
        timName     name of TIM
    outputs:
        (None)
    """

    _ocsm.tim_quit.argtypes = [ctypes.c_char_p]
    _ocsm.tim_quit.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    status = _ocsm.tim_quit(timName)
    _processStatus(status, "TimQuit")

    return

# ======================================================================

def TimSave(timName):
    """
    ocsm.TimSave - saves data and quits a TIM

    inputs:
        timName     name of TIM
    outputs:
        (None)
    """

    _ocsm.tim_save.argtypes = [ctypes.c_char_p]
    _ocsm.tim_save.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    status = _ocsm.tim_save(timName)
    _processStatus(status, "TimSave")

    return

# ======================================================================

def TimMesg(timName, mesg):
    """
    ocsm.TimMesg - sends a message to a TIM

    inputs:
        timName     name of TIM
        mesg        message to send
    outputs:
        (None)
    """

    _ocsm.tim_mesg.argtypes = [ctypes.c_char_p,
                               ctypes.c_char_p]
    _ocsm.tim_mesg.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    if (isinstance(mesg, str)):
        mesg = mesg.encode()

    status = _ocsm.tim_mesg(timName, mesg)
    _processStatus(status, "TimMesg")

    return

# ======================================================================

def TimHold(timName, overlay):
    """
    ocsm.TimHold - hold timName until overlay finished

    inputs:
        timName     name of TIM to be held
        overlay     name of overlay
    outputs:
        (None)
    """

    _ocsm.tim_hold.argtypes = [ctypes.c_char_p,
                               ctypes.c_char_p]
    _ocsm.tim_hold.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    if (isinstance(overlay, str)):
        overlay = overlay.encode()

    status = _ocsm.tim_hold(timName, overlay)
    _processStatus(status, "TimHold")

    return

# ======================================================================

def TimBcst(timName, text):
    """
    ocsm.TimBcst - broadcasts a message to all browsers

    inputs:
        timName     name of (any) TIM
        text        text to broadcast
    outputs:
        (None)
    """

    _ocsm.tim_bcst.argtypes = [ctypes.c_char_p,
                               ctypes.c_char_p]
    _ocsm.tim_bcst.restype  =  ctypes.c_int

    if (isinstance(timName, str)):
        timName = timName.encode()

    if (isinstance(text, str)):
        text = text.encode()

    status = _ocsm.tim_bcst(timName, text)
    _processStatus(status, "TimBcst")

    return

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

