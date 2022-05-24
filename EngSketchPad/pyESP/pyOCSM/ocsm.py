###################################################################
#                                                                 #
# pyOCSM --- Python version of OpenCSM API                        #
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

# get the value of _ESP_ROOT
try:
    _ESP_ROOT = os.environ["ESP_ROOT"]
except:
    raise RuntimeError("ESP_ROOT must be set -- Please fix the environment...")

# load the shared library
if sys.platform.startswith('darwin'):
    _ocsm = ctypes.CDLL(_ESP_ROOT + "/lib/libocsm.dylib")
elif sys.platform.startswith('linux'):
    _ocsm = ctypes.CDLL(_ESP_ROOT + "/lib/libocsm.so")
elif sys.platform.startswith('win32'):
    if sys.version_info.major == 3 and sys.version_info.minor < 8:
        _ocsm = ctypes.CDLL(_ESP_ROOT + "\\lib\\ocsm.dll")
    else:
        _ocsm = ctypes.CDLL(_ESP_ROOT + "\\lib\\ocsm.dll", winmode=0)
else:
    raise IOError("Unknown platform: " + sys.platform)

# ======================================================================

def Version():
    """
    Version - return current version

    inputs:
        (None)
    outputs:
        imajor      major version number
        iminor      minor version number
    """
    _ocsm.ocsmVersion.argtypes = [ctypes.POINTER(ctypes.c_int),
                                  ctypes.POINTER(ctypes.c_int)]
    _ocsm.ocsmVersion.restype  =  ctypes.c_int

    imajor = ctypes.c_int()
    iminor = ctypes.c_int()

    status = _ocsm.ocsmVersion(imajor, iminor)
    _processStatus(status, "Version")

    return (imajor.value, iminor.value)

# ======================================================================

def SetOutLevel(ilevel):
    """
    ocsm.SetOutLevel - set output level

    inputs:
        ilevel      0  errors only
                    1  nominal (default)
                    2  debug
    outputs:
        (None)
    """
    _ocsm.ocsmSetOutLevel.argtypes = [ctypes.c_int]
    _ocsm.ocsmSetOutLevel.restype  =  ctypes.c_int

    status = _ocsm.ocsmSetOutLevel(ilevel)
    _processStatus(status, "SetOutLevel")

    return

# ======================================================================

def PrintEgo(theEgo):
    """
    ocsm.PrintEgo - print the contents of an EGADS ego

    inputs:
        theEgo      the EGADS ego
    outputs:
        (None)
    """
    _ocsm.ocsmPrintEgo.argtypes = [egads.c_ego]

    _ocsm.ocsmPrintEgo(theEgo.py_to_c())

    return

# ======================================================================

def _ocsmFree(modl=None):
    """
    ocsm._ocsmFree - free a MODL (private method)

    inputs:
        (None)
    output:
        (None)
    """
    _ocsm.ocsmFree.argtypes = [ctypes.c_void_p]
    _ocsm.ocsmFree.restype  =  ctypes.c_int

    status = _ocsm.ocsmFree(modl)
    _processStatus(status, "Free")

    return

# ======================================================================
# ocsm constants
# ======================================================================

# Branch activities
ACTIVE       = 300
SUPPRESSED   = 301
INACTIVE     = 302
DEFERRED     = 303

# Body types
SOLID_BODY   = 400
SHEET_BODY   = 401
WIRE_BODY    = 402
NODE_BODY    = 403
NULL_BODY    = 404

# Parameter types
DESPMTR      = 500
CFGPMTR      = 501
CONPMTR      = 502
LOCALVAR     = 503
OUTPMTR      = 504

# Selector types
NODE         = 600
EDGE         = 601
FACE         = 602
BODY         = 603

# ======================================================================
# Ocsm class
# ======================================================================

class Ocsm(object):
    """
    official Python bindings for the OpenCSM API

    written by John Dannenhoffer (author of OpenCSM)
    """

    def __init__(self, init):
        """
        ocsm.Load - create a MODL by reading a .csm file

        inputs:
            init        file to be read (with .csm) or pointer to MODL
        output:
            (None)
        """
        _ocsm.ocsmLoad.argtypes = [ctypes.c_char_p,
                                   ctypes.POINTER(ctypes.c_void_p)]
        _ocsm.ocsmLoad.restype  =  ctypes.c_int

        _ocsm.ocsmLoadFromModel.argtypes = [ctypes.c_void_p,
                                            ctypes.POINTER(ctypes.c_void_p)]
        _ocsm.ocsmLoadFromModel.restype  =  ctypes.c_int

        if (isinstance(init, str)):
            filename = init.encode()

            self._modl    = ctypes.c_void_p()
            self._mesgCB  = None
            self._sizeCB  = None
            self._context = None
            status        = _ocsm.ocsmLoad(filename, self._modl)
            _processStatus(status, "Load")

            # create finalizer to ensure _ocsmFree is called during garbage collection
            self._pyowned  = True
            self._finalize = egads.finalize(self, _ocsmFree, self._modl)

        elif (isinstance(init, egads.ego)):
            emodel = init.py_to_c()

            self._modl    = ctypes.cast(emodel, ctypes.c_void_p)
            self._mesgCB  = None
            self._sizeCB  = None
            self._context = None
            status        = _ocsm.ocsmLoadFromModel(emodel, self._modl)
            _processStatus(status, "LoadFromModel")

            # create finalizer to ensure _ocsmFree is called during garbage collection
            self._pyowned  = True
            self._finalize = egads.finalize(self, _ocsmFree, self._modl)

        else:
            self._modl    = ctypes.c_void_p()
            self._mesgCB  = None
            self._sizeCB  = None
            self._context = None

            self._modl = init

            # we do not want to call _ocsmFree since we do not own the MODL
            self._pyowned  = False
            self._finalize = None

        return

# ======================================================================

    def __del__(self):
        """
        ocsm.__del__ - free up the MODL

        Note: good programming practice would be to call ocsm.Free()
              yourself rather than waiting for Python to do the
              garbage collection for you!!!
        """

        # explicitly call the finalizer to invoke _ocsmFree
        if (self._finalize is not None):
            self._finalize()

        return

# ======================================================================

    def __copy__(self):
        """
        ocsm.__copy__ - do not use this.  Use ocsm.Copy to make a partial copy
        """
        raise OcsmError("Do not use this.  Use ocsm.Copy to make a partial copy to inputs only", "__copy__")

# ======================================================================

    def __deepcopy__(self):
        """
        ocsm.__deepcopy__ - do not use this.  Use ocsm.Copy to make a partial copy
        """
        raise OcsmError("Do not use this.  Use ocsm.Copy to make a partial copy to inputs only", "__deepcopy__")

# ======================================================================

    def Load(self, filename):
        """
        ocsm.Load - create a MODL by using __init__ method
        """
        raise OcsmError("Use __init__ method instead", "Load")

# ======================================================================

    def LoadDict(self, dictname):
        """
        ocsm.LoadDist - load dictionary from dictname

        inputs:
            dictname    file that contains dictionary
        outputs:
            (None)
        """
        _ocsm.ocsmLoadDict.argtypes = [ctypes.c_void_p,
                                       ctypes.c_char_p]
        _ocsm.ocsmLoadDict.restype  =  ctypes.c_int

        if (isinstance(dictname, str)):
            dictname = dictname.encode()

        status    = _ocsm.ocsmLoadDict(self._modl, dictname)
        _processStatus(status, "LoadDict")

        return

# ======================================================================

    def UpdateDespmtrs(self, filename):
        """
        ocsm.UpdateDespmtrs - update CFGPMTRs and DESPMTRs from filename

        inputs:
            filename    file that contains DESPMTRs
        outputs:
            (None)
        """
        _ocsm.ocsmUpdateDespmtrs.argtypes = [ctypes.c_void_p,
                                             ctypes.c_char_p]
        _ocsm.ocsmUpdateDespmtrs.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmUpdateDespmtrs(self._modl, filename)
        _processStatus(status, "UpdateDespmtrs")

        return

# ======================================================================

    def GetFilelist(self):
        """
        ocsm.GetFilelist - get a list of all .csm, .cpc, and .udc files

        inputs:
            filelist    bar-separated list of files
        outputs:
            (None)
        """
        _ocsm.ocsmGetFilelist.argtypes = [ctypes.c_void_p,
                                          ctypes.POINTER(ctypes.c_char_p)]
        _ocsm.ocsmGetFilelist.restype  =  ctypes.c_int

        pfilelist = ctypes.c_char_p()

        status = _ocsm.ocsmGetFilelist(self._modl, pfilelist)
        _processStatus(status, "GetFilelist")

        filelist = pfilelist.value.decode("utf-8")
        egads.free(pfilelist)

        return filelist

# ======================================================================

    def Save(self, filename):
        """
        ocsm.Save - save a MODL to a file

        inputs:
            filename    file to be written (with extension)
        output:
            (None)
        """
        _ocsm.ocsmSave.argtypes = [ctypes.c_void_p,
                                   ctypes.c_char_p]
        _ocsm.ocsmSave.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmSave(self._modl, filename)
        _processStatus(status, "Save")

        return

# ======================================================================

    def SaveDespmtrs(self, filename):
        """
        ocsm.SaveDespmtrs - save CFGPMTRs and DESPMTRs to file

        inputs:
            filename    file to be written
        outputs:
            (None)
        """
        _ocsm.ocsmSaveDespmtrs.argtypes = [ctypes.c_void_p,
                                           ctypes.c_char_p]
        _ocsm.ocsmSaveDespmtrs.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmSaveDespmtrs(self._modl, filename)
        _processStatus(status, "SaveDespmtrs")

        return

# ======================================================================

    def Copy(self):
        """
        ocsm.Copy - copy a MODL (but not the Bodys)

        inputs:
            (None)
        outputs:
            (None)
        side effect:
            a new Ocsm MODL is created
        """
        _ocsm.ocsmCopy.argtypes = [ctypes.c_void_p,
                                   ctypes.POINTER(ctypes.c_void_p)]
        _ocsm.ocsmCopy.restype  =  ctypes.c_int

        new = ctypes.c_void_p()
        status = _ocsm.ocsmCopy(self._modl, new)
        _processStatus(status, "Copy")

        cls              = self.__class__
        newModl          = cls.__new__(cls)
        newModl._modl    = new
        newModl._mesgCB  = self._mesgCB
        newModl._sizeCB  = self._sizeCB
        newModl._context = self._context
        newModl._pyowned = True

        # create a new finalizer for the new Ocsm instance
        if (self._finalize is not None):
            newModl._finalize = egads.finalize(newModl, _ocsmFree, newModl._modl)

        return newModl

# ======================================================================

    def Free(self):
        """
        ocsm.Free - free up all strorage associated with a MODL

        inputs:
            (None)
        outputs:
            (None)
        """

        # explicitly call the finalizer to invoke _ocsmFree
        self._finalize()

        # set member variables to None to possibly invole their garbage collection
        self._finalize = None
        self._modl     = None
        self._context  = None
        self._pyowned  = False

        return

# ======================================================================

    def Info(self):
        """
        ocsm.Info - get info about a MODL

        inputs:
            (None)
        outputs:
            nbrch       number of Branches
            npmtr       number of Parameters
            nbody       number of Bodys
        """
        _ocsm.ocsmInfo.argtypes = [ctypes.c_void_p,
                                   ctypes.POINTER(ctypes.c_int),
                                   ctypes.POINTER(ctypes.c_int),
                                   ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmInfo.restype  =  ctypes.c_int

        nbrch  = ctypes.c_int()
        npmtr  = ctypes.c_int()
        nbody  = ctypes.c_int()
        status = _ocsm.ocsmInfo(self._modl, nbrch, npmtr, nbody)
        _processStatus(status, "Info")

        return (nbrch.value, npmtr.value, nbody.value)

# ======================================================================

    def Check(self):
        """
        ocsm.Check - check that Branches are properly ordered

        inputs:
            (None)
        outputs:
            (None)
        """
        _ocsm.ocsmCheck.argtypes = [ctypes.c_void_p]
        _ocsm.ocsmCheck.restype  =  ctypes.c_int

        status = _ocsm.ocsmCheck(self._modl)
        _processStatus(status, "Check")

        return

# ======================================================================

    def RegMesgCB(self, callback):
        """
        ocsm.RegMesgCB - register a callback function for exporting messages

        inputs:
            callback    callback function (modl, text)
        outputs:
            (None)
        """
        _ocsm.ocsmRegMesgCB.argtypes = [ctypes.c_void_p,
                                        ctypes.CFUNCTYPE(None, ctypes.c_char_p)]
        _ocsm.ocsmRegMesgCB.restype  =  ctypes.c_int

        MESGFUNC = ctypes.CFUNCTYPE(None, ctypes.c_char_p)
        self._mesgCB = MESGFUNC(callback)

        status = _ocsm.ocsmRegMesgCB(self._modl, self._mesgCB)
        _processStatus(status, "RegMesgCB")

        return

# ======================================================================

    def RegSizeCB(self, callback):
        """
        ocsm.RegSizeCB - register a callback function for DESPMTR size changes

        inputs:
            callback    callback function (modl, ipmtr, nrow, ncol)
        outputs:
            (None)
        """
        _ocsm.ocsmRegSizeCB.argtypes = [ctypes.c_void_p,
                                        ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int)]
        _ocsm.ocsmRegSizeCB.restype  =  ctypes.c_int

        SIZEFUNC = ctypes.CFUNCTYPE(None, ctypes.c_void_p, ctypes.c_int, ctypes.c_int, ctypes.c_int)
        self._sizeCB = SIZEFUNC(callback)

        status = _ocsm.ocsmRegSizeCB(self._modl, self._sizeCB)
        _processStatus(status, "RegSizeCB")

        return

# ======================================================================

    def Build(self, buildTo, nbody):
        """
        ocsm.Build - build Bodys by executing the MODL up to a given Branch

        inputs:
            buildTo     last Branch to execute (or 0 for all, or -1 for no recycling)
            nbody       number of entres allowed in body
        outputs:
            builtTo     last Branch executed successfully
            nbody       number of Bodys on the stack
            body        array of Bodys on the stack
        """
        _ocsm.ocsmBuild.argtypes = [ctypes.c_void_p,
                                    ctypes.c_int,
                                    ctypes.POINTER(ctypes.c_int),
                                    ctypes.POINTER(ctypes.c_int),
                                    ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmBuild.restype  =  ctypes.c_int

        builtTo = ctypes.c_int()
        nbody_  = ctypes.c_int(nbody)

        if (nbody == 0):
            body = None
        else:
            body = (ctypes.c_int * nbody)()

        status = _ocsm.ocsmBuild(self._modl, buildTo, builtTo, nbody_, body)
        _processStatus(status, "Build")

        if (self._context is None and self._pyowned == True):

            # create an egads.Context that owns the ego and which will be properly
            #    closed when the self._context is garbage collected
            self._context = egads.c_to_py(self.GetEgo(0, BODY, 2).py_to_c(), deleteObject=True)

            # create a new finalizer after the context is created
            if (self._finalize is not None):
                self._finalize.detach()

                self._finalize = egads.finalize(self, _ocsmFree, self._modl)

        # if Bodys are returned, convert into a list
        if (body is not None):
            body = list(body[0:nbody_.value])

        return (builtTo.value, nbody_.value, body)

# ======================================================================

    def Perturb(self, npmtrs, ipmtrs, irows, icols, values):
        """
        ocsm.Perturb - create a perturbed MODL

        inputs:
            npmtrs   number of perturbed Parameters (or 0 to remove)
            ipmtrs   list of Parameter indices (1:npmtr)
            irows    list of row       indices (1:nrow)
            icols    list of column    indices (1:ncol)
            values   list of perturbed values
        outputs:
            (None)
        """
        _ocsm.ocsmPerturb.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmPerturb.restype  =  ctypes.c_int

        if (npmtrs > 0):
            ipmtrs_ = (ctypes.c_int    * (npmtrs))(*ipmtrs)
            irows_  = (ctypes.c_int    * (npmtrs))(*irows )
            icols_  = (ctypes.c_int    * (npmtrs))(*icols )
            values_ = (ctypes.c_double * (npmtrs))(*values)

            status = _ocsm.ocsmPerturb(self._modl, npmtrs, ipmtrs_, irows_, icols_, values_)
            _processStatus(status, "Perturb")
        else:
            status = _ocsm.ocsmPerturb(self._modl, 0, None, None, None, None)
            _processStatus(status, "Perturb")

        return

# ======================================================================

    def UpdateTess(self, ibody, filename):
        """
        ocsm.UpdateTess - update a tessellation from a file

        inputs:
            filename    name of .tess file
        outputs:
            (None)
        """
        _ocsm.ocsmUpdateTess.argtypes = [ctypes.c_void_p,
                                         ctypes.c_int,
                                         ctypes.c_char_p]
        _ocsm.ocsmUpdateTess.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmUpdateTess(self._modl, ibody, filename)
        _processStatus(status, "UpdateTess")

        return

# ======================================================================

    def NewBrch(self, iafter, type, filename, linenum, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9):
        """
        ocsm.NewBrch - create a new Branch

        inputs:
            iafter      Branch index (0-nbrch) after which to add
            type        Branch type
            filename    filename when Branch is defined
            linenum     linenumber where Branch is defined (bias-1)
            arg1        Argument 1
            arg2        Argument 2
            arg3        Argument 3
            arg4        Argument 4
            arg5        Argument 5
            arg6        Argument 6
            arg7        Argument 7
            arg8        Argument 8
            arg9        Argument 9
        outputs:
            (None)
        """
        _ocsm.ocsmNewBrch.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmNewBrch.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()
        if (isinstance(arg1, str)):
            arg1 = arg1.encode()
        if (isinstance(arg2, str)):
            arg2 = arg2.encode()
        if (isinstance(arg3, str)):
            arg3 = arg3.encode()
        if (isinstance(arg4, str)):
            arg4 = arg4.encode()
        if (isinstance(arg5, str)):
            arg5 = arg5.encode()
        if (isinstance(arg6, str)):
            arg6 = arg6.encode()
        if (isinstance(arg7, str)):
            arg7 = arg7.encode()
        if (isinstance(arg8, str)):
            arg8 = arg8.encode()
        if (isinstance(arg9, str)):
            arg9 = arg9.encode()

        status = _ocsm.ocsmNewBrch(self._modl, iafter, type, filename, linenum, arg1, arg2, arg3, arg4, arg5, arg6, arg7, arg8, arg9)
        _processStatus(status, "NewBrch")

        return

# ======================================================================

    def GetBrch(self, ibrch):
        """
        ocsm.GetBrch - get info about a Branch

        inputs:
            ibrch       Branch index (1:nbrch)
        outputs:
            type        Branch type
            bclass      Branch class
            actv        Branch activity
            ichld       ibrch of child (or 0 if root)
            ileft       ibrch of left parent (or 0)
            irite       ibrch of rite parent (or 0)
            narg        number of Arguments
            nattr       number of Attributes
        """
        _ocsm.ocsmGetBrch.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmGetBrch.restype  =  ctypes.c_int

        type   = ctypes.c_int()
        bclass = ctypes.c_int()
        actv   = ctypes.c_int()
        ichld  = ctypes.c_int()
        ileft  = ctypes.c_int()
        irite  = ctypes.c_int()
        narg   = ctypes.c_int()
        nattr  = ctypes.c_int()

        status = _ocsm.ocsmGetBrch(self._modl, ibrch, type, bclass, actv, ichld, ileft, irite, narg, nattr)
        _processStatus(status, "GetBrch")

        return (type.value, bclass.value, actv.value, ichld.value, ileft.value, irite.value, narg.value, nattr.value)

# ======================================================================

    def SetBrch(self, ibrch, actv):
        """
        ocsm.SetBrch - set activity for a Branch

        inputs:
            ibrch       Branch index (1:nbrch)
            actv        Branch activity
        outputs:
            (none)
        """
        _ocsm.ocsmSetBrch.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int]
        _ocsm.ocsmSetBrch.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetBrch(self._modl, ibrch, actv)
        _processStatus(status, "SetBrch")

        return

# ======================================================================

    def DelBrch(self,ibrch):
        """
        ocsm.DelBrch - delete a Branch (or whole Sketch if SKBEG)

        inputs:
            ibrch       Branch index (1:nbrch)
        outputs:
            (None)
        """
        _ocsm.ocsmDelBrch.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int]
        _ocsm.ocsmDelBrch.restype  =  ctypes.c_int

        status = _ocsm.ocsmDelBrch(self._modl, ibrch)
        _processStatus(status, "DelBrch")

        return

# ======================================================================

    def PrintBrchs(self, filename):
        """
        ocsm.PrintBrchs - print Branches to a file

        inputs:
            filename    file to which output is appended (or "" for stdout)
        outputs:
            (None)
        """
        _ocsm.ocsmPrintBrchs.argtypes = [ctypes.c_void_p,
                                         ctypes.c_char_p]
        _ocsm.ocsmPrintBrchs.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmPrintBrchs(self._modl, filename)
        _processStatus(status, "PrintBrchs")

        return

# ======================================================================

    def GetArg(self, ibrch, iarg):
        """
        ocsm.GetArg - get an Argument for a Branch

        inputs:
            ibrch       Branch   index (1:nbrch)
            iarg        Argument index (1:narg)
        outputs:
            defn        Argument definition
            value       Argument value
            dot         Arfument velocity
        """
        _ocsm.ocsmGetArg.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_char_p,
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetArg.restype  =  ctypes.c_int

        defn  = ctypes.create_string_buffer(257)
        value = ctypes.c_double()
        dot   = ctypes.c_double()

        status = _ocsm.ocsmGetArg(self._modl, ibrch, iarg, defn, value, dot)
        _processStatus(status, "GetArg")

        return (defn.value.decode("utf-8"), value.value, dot.value)

# ======================================================================

    def SetArg(self, ibrch, iarg, defn):
        """
        ocsm.SetArg - set an Argument for a Branch

        inputs:
            ibrch       Branch   index (1:nbrch)
            iarg        Argument index (1:narg)
            defn        Argument definition
        outputs:
            (None)
        """
        _ocsm.ocsmSetArg.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_char_p]
        _ocsm.ocsmSetArg.restype  =  ctypes.c_int

        if (isinstance(defn, str)):
            defn = defn.encode()

        status = _ocsm.ocsmSetArg(self._modl, ibrch, iarg, defn)
        _processStatus(status, "SetArg")

        return

# ======================================================================

    def RetAttr(self, ibrch, iattr):
        """
        ocsm.RetAttr - return an Attribute for a Branch by index

        inputs:
            ibrch       Branch    index (1:nbrch)
            iattr       Attribute index (1:nbrch)
        outputs:
            aname       Attribute name
            avalue      Attribute value
        """
        _ocsm.ocsmRetAttr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmRetAttr.restype  =  ctypes.c_int

        aname  = ctypes.create_string_buffer(257)
        avalue = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmRetAttr(self._modl, ibrch, iattr, aname, avalue)
        _processStatus(status, "RetAttr")

        return (aname.value.decode("utf-8"), avalue.value.decode("utf-8"))

# ======================================================================

    def GetAttr(self, ibrch, aname):
        """
        ocsm.RetAttr - return an Attribute for a Branch by name

        inputs:
            ibrch       Branch    index (1:nbrch)
            aname       Attribute name
        outputs:
            avalue      Attribute value
        """
        _ocsm.ocsmGetAttr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmGetAttr.restype  =  ctypes.c_int

        if (isinstance(aname, str)):
            aname = aname.encode()

        avalue = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmGetAttr(self._modl, ibrch, aname, avalue)
        _processStatus(status, "GetAttr")

        return avalue.value.decode("utf-8")

# ======================================================================

    def SetAttr(self, ibrch, aname, avalue):
        """
        ocsm.SetAttr - set an Attribute for a Branch

        inputs:
            ibrch       Branch index (1:nbrch) or 0 for global
            aname       Attribute name
            avalue      Attribute value (or blank to delete)
        outputs:
            (None)
        """
        _ocsm.ocsmSetAttr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmSetAttr.restype  =  ctypes.c_int

        if (isinstance(aname, str)):
            aname = aname.encode()
        if (isinstance(avalue, str)):
            avalue = avalue.encode()

        status = _ocsm.ocsmSetAttr(self._modl, ibrch, aname, avalue)
        _processStatus(status, "SetAttr")

        return

# ======================================================================

    def RetCsys(self, ibrch, icsys):
        """
        ocsm.RetCsys - return a Csystem for a Branch by index

        inputs:
            ibrch       Branch  index (1:nbrch)
            icsys       Csystem index (1:nattr)
        outputs:
            cname       Csystem name
            cvalue      Csystem value
        """
        _ocsm.ocsmRetCsys.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmRetCsys.restype  =  ctypes.c_int

        cname  = ctypes.create_string_buffer(257)
        cvalue = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmRetCsys(self._modl, ibrch, icsys, cname, cvalue)
        _processStatus(status, "RetCsys")

        return (cname.value.decode("utf-8"), cvalue.value.decode("utf-8"))

# ======================================================================

    def GetCsys(self, ibrch, cname):
        """
        ocsm.GetCsys - return a Csystem for a Branch by name

        inputs:
            ibrch       Branch index  (1:nbrch)
            cname       Csystem name
        outputs:
            cvalue      Csystem value
        """
        _ocsm.ocsmGetCsys.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmGetCsys.restype  =  ctypes.c_int

        if (isinstance(cname, str)):
            cname = cname.encode()

        cvalue = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmGetCsys(self._modl, ibrch, cname, cvalue)
        _processStatus(status, "GetCsys")

        return cvalue.value.decode("utf-8")

# ======================================================================

    def SetCsys(self, ibrch, cname, cvalue):
        """
        ocsm.SetCsys - set a Csystem for a Branch

        inputs:
            ibrch       Branch index (1:nbrch)
            cname       Csystem name
            cvalue      Csystem value
        outputs:
            (None)
        """
        _ocsm.ocsmSetCsys.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p,
                                      ctypes.c_char_p]
        _ocsm.ocsmSetCsys.restype  =  ctypes.c_int

        if (isinstance(cname, str)):
            cname = cname.encode()
        if (isinstance(cvalue, str)):
            cvalue = cvalue.encode()

        status = _ocsm.ocsmSetCsys(self._modl, ibrch, cname, cvalue)
        _processStatus(status, "SetCsys")

        return

# ======================================================================

    def PrintAttrs(self, filename):
        """
        ocsm.PrintAttrs - print global Attributes to file

        inputs:
            filename    file to which output is appended (or "" for stdout)
        outputs:
            (None)
        """
        _ocsm.ocsmPrintAttrs.argtypes = [ctypes.c_void_p,
                                         ctypes.c_char_p]
        _ocsm.ocsmPrintAttrs.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmPrintAttrs(self._modl, filename)
        _processStatus(status, "PrintAttrs")

        return

# ======================================================================

    def GetName(self, ibrch):
        """
        ocsm.GetName - get the name of a Branch

        inputs:
            ibrch       Branch index (1:nbrch)
        outputs:
            name        Brchch name
        """
        _ocsm.ocsmGetName.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p]
        _ocsm.ocsmGetName.restype  =  ctypes.c_int

        name = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmGetName(self._modl, ibrch, name)
        _processStatus(status, "GetName")

        return name.value.decode("utf-8")

# ======================================================================

    def SetName(self, ibrch, name):
        """
        ocsm.SetName - set the name for a Branch

        inputs:
            ibrch       Branch index (1:nbrch)
            name        Branch name
        outputs:
            (None)
        """
        _ocsm.ocsmSetName.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_char_p]
        _ocsm.ocsmSetName.restype  =  ctypes.c_int

        if (isinstance(name, str)):
            name = name.encode()

        status = _ocsm.ocsmSetName(self._modl, ibrch, name)
        _processStatus(status, "SetName")

        return

# ======================================================================

    def GetSketch(self, ibrch, maxlen):
        """
        ocsm.GetSketch - get string data associated with a Sketch

        not implemented
        """
        _ocsm.ocsmGetSketch.argtypes = [ctypes.c_void_p,
                                        ctypes.c_int,
                                        ctypes.c_int,
                                        ctypes.c_char_p,
                                        ctypes.c_char_p,
                                        ctypes.c_char_p,
                                        ctypes.c_char_p]
        _ocsm.ocsmGetSketch.restype  =  ctypes.c_int

        _processStatus(-298, "GetSketch")

        return # (begs, vars, cons, segs)

# ======================================================================

    def SolveSketch(self, vars_in, cons, vars_out):
        """
        ocsm.SolveSketch - solve for new Sketch variables

        not implemented
        """
        _ocsm.ocsmSolveSketch.argtypes = [ctypes.c_void_p,
                                          ctypes.c_char_p,
                                          ctypes.c_char_p,
                                          ctypes.c_char_p]
        _ocsm.ocsmSolveSketch.restype  =  ctypes.c_int

        _processStatus(-298, "SolveSketch")

        return # vars_out

# ======================================================================

    def SaveSketch(self, ibrch, vars, cons, segs):
        """
        ocsm.SaveSketch - overwrite Branches associated with a Sketch

        not implemented
        """
        _ocsm.ocsmSaveSketch.argtypes = [ctypes.c_void_p,
                                         ctypes.c_int,
                                         ctypes.c_char_p,
                                         ctypes.c_char_p,
                                         ctypes.c_char_p]
        _ocsm.ocsmSaveSketch.restype  =  ctypes.c_int

        _processStatus(-298, "SaveSketch")

        return

# ======================================================================

    def NewPmtr(self, name, type, nrow, ncol):
        """
        ocsm.NewPmtr - create a new Parameter

        inputs:
            name        Parameter name
            type        ocsm.DESPMTR
                        ocsm.CFGPMTR
                        ocsm.CONPMTR
                        ocsm.LOCALVAR
                        ocsm.OUTPMTR
            nrow        number of rows
            ncol        number of columns
        outputs:
            (None)
        """
        _ocsm.ocsmNewPmtr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_char_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int]
        _ocsm.ocsmNewPmtr.restype = ctypes.c_int

        if (isinstance(name, str)):
            name = name.encode()

        ipmtr  = ctypes.c_int()
        status = _ocsm.ocsmNewPmtr(self._modl, name, type, nrow, ncol)
        _processStatus(status, "NewPmtr")

        return

# ======================================================================

    def DelPmtr(self, ipmtr):
        """
        ocsm.DelPmtr = delete a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
        outputs:
            (None)
        """
        _ocsm.ocsmDelPmtr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int]
        _ocsm.ocsmDelPmtr.restype  =  ctypes.c_int

        status = _ocsm.ocsmDelPmtr(self._modl, ipmtr)
        _processStatus(status, "DelPmtr")

        return

# ======================================================================

    def FindPmtr(self, name, type, nrow, ncol):
        """
        ocsm.FindPmtr - find (or create) a Parameter

        inputs:
            name        Parameter name
            type        ocsm.DESPMTR
                        ocsm.CFGPMTR
                        ocsm.CONPMTR
                        ocsm.LOCALVAR
                        ocsm.OUTPMTR
            nrow        number of rows
            ncol        number of columns
        outputs:
            (None)
        """
        _ocsm.ocsmFindPmtr.argtypes = [ctypes.c_void_p,
                                       ctypes.c_char_p,
                                       ctypes.c_int,
                                       ctypes.c_int,
                                       ctypes.c_int,
                                       ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmFindPmtr.restype = ctypes.c_int

        if (isinstance(name, str)):
            name = name.encode()

        ipmtr  = ctypes.c_int()
        status = _ocsm.ocsmFindPmtr(self._modl, name, type, nrow, ncol, ipmtr)
        _processStatus(status, "FindPmtr")

        return ipmtr.value

# ======================================================================

    def GetPmtr(self, ipmtr):
        """
        ocsm.GetPmtr - get info about a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
        outputs:
            type        ocsm.DESPMTR, ocsm.CFGPMTR, ocsm.CONPMTR,
                            ocsm.LOCALVAL, or ocsm.OUTPMTR
            nrow        number of rows
            ncol        number of columns
            name        Parameter name
        """
        _ocsm.ocsmGetPmtr.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.c_char_p]
        _ocsm.ocsmGetPmtr.restype  =  ctypes.c_int

        type  = ctypes.c_int()
        nrow  = ctypes.c_int()
        ncol  = ctypes.c_int()
        name  = ctypes.c_char_p(b"\0"*65)
        status = _ocsm.ocsmGetPmtr(self._modl, ipmtr, type, nrow, ncol, name)
        _processStatus(status, "GetPmtr")

        return (type.value, nrow.value, ncol.value, name.value.decode("utf-8"))

# ======================================================================

    def PrintPmtrs(self, filename):
        """
        ocsm.PrintPmtrs - print DESPMTRs and OUTPMTRs to file

        inputs:
            filename    file to which output is appended (or "" for stdout)
        outputs:
            (None)
        """
        _ocsm.ocsmPrintPmtrs.argtypes = [ctypes.c_void_p,
                                         ctypes.c_char_p]
        _ocsm.ocsmPrintPmtrs.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmPrintPmtrs(self._modl, filename)
        _processStatus(status, "PrintPmtrs")

        return

# ======================================================================

    def GetValu(self, ipmtr, irow, icol):
        """
        ocsm.GetValue - get the Value of a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row       index (1:nrow)
            icol        column    index (1:ncol)
        outputs:
            value       Parameter Value
            dot         Parameter velocity
        """
        _ocsm.ocsmGetValu.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetValu.restype  =  ctypes.c_int

        value  = ctypes.c_double()
        dot    = ctypes.c_double()
        status = _ocsm.ocsmGetValu(self._modl, ipmtr, irow, icol, value, dot)
        _processStatus(status, "GetValu")

        return (value.value, dot.value)

# ======================================================================

    def GetValuS(self, ipmtr):
        """
        ocsm.GetValue - get the Value of a string Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
        outputs:
            str         Parameter Value
        """
        _ocsm.ocsmGetValuS.argtypes = [ctypes.c_void_p,
                                       ctypes.c_int,
                                       ctypes.c_char_p]
        _ocsm.ocsmGetValuS.restype  =  ctypes.c_int

        str = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmGetValuS(self._modl, ipmtr, str)
        _processStatus(status, "GetValuS")

        return str.value.decode("utf-8")

# ======================================================================

    def SetValu(self, ipmtr, irow, icol, defn):
        """
        ocsm.SetValu - set a Value for a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row       index (1:nrow)
            icol        column    indes (1:ncol)
            defn        definieiotn of Value
        outputs:
            (None)
        """
        _ocsm.ocsmSetValu.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_char_p]
        _ocsm.ocsmSetValu.restype  =  ctypes.c_int

        if (isinstance(defn, str)):
            defn = defn.encode()

        status = _ocsm.ocsmSetValu(self._modl, ipmtr, irow, icol, defn)
        _processStatus(status, "SetValu")

        return

# ======================================================================

    def SetValuD(self, ipmtr, irow, icol, value):
        """
        ocsm.SetValuD - set the (double) value of a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row       index (1:nrow)
            icol        column    index (1:ncol)
            value       value to set
        outputs:
            (None)
        """
        _ocsm.ocsmSetValuD.argtypes = [ctypes.c_void_p,
                                       ctypes.c_int,
                                       ctypes.c_int,
                                       ctypes.c_int,
                                       ctypes.c_double]
        _ocsm.ocsmSetValuD.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetValuD(self._modl, ipmtr, irow, icol, value)
        _processStatus(status, "SetValuD")

        return

# ======================================================================

    def GetBnds(self, ipmtr, irow, icol):
        """
        ocsm.GetBnds - get the Bounds of a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row       index (1:nrow)
            icol        column    index (1:ncol)
        outputs:
            lbound      lower Bound
            ubound      upper Bound
        """
        _ocsm.ocsmGetBnds.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetBnds.restype  =  ctypes.c_int

        lbound = ctypes.c_double()
        ubound = ctypes.c_double()

        status = _ocsm.ocsmGetBnds(self._modl, ipmtr, irow, icol, lbound, ubound)
        _processStatus(status, "GetBnds")

        return (lbound.value, ubound.value)

# ======================================================================

    def SetBnds(self, ipmtr, irow, icol, lbound, ubound):
        """
        ocsm.SetBnds - set teh Bounds of a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row       index (1:nrow)
            icol        column    index (1:ncol)
            lbound      lower Bound
            ubound      upper Bound
        outputs:
            (None)
        """
        _ocsm.ocsmSetBnds.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_double,
                                      ctypes.c_double]
        _ocsm.ocsmSetBnds.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetBnds(self._modl, ipmtr, irow, icol, lbound, ubound)
        _processStatus(status, "SetBnds")

        return

# ======================================================================

    def SetDtime(self, dtime):
        """
        ocsm.SetDtime - set sensitivity FD time step (or select analytic)

        inputs:
            dtime       time step (or 0 to choose analytic)
        outputs:
            (None)
        """
        _ocsm.ocsmSetDtime.argtypes = [ctypes.c_void_p,
                                       ctypes.c_double]
        _ocsm.ocsmSetDtime.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetDtime(self._modl, dtime)
        _processStatus(status, "SetDtime")

        return

# ======================================================================

    def SetVel(self, ipmtr, irow, icol, defn):
        """
        ocsm.SetVel - set the velocity for a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row index       (1:nrow)
            icol        column index    (1:ncol)
            defn        definition of Velocity
        outputs:
            (None)
        """
        _ocsm.ocsmSetVel.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_char_p]
        _ocsm.ocsmSetVel.restype  =  ctypes.c_int

        if (isinstance(defn, str)):
            defn = defn.encode()

        status = _ocsm.ocsmSetVel(self._modl, ipmtr, irow, icol, defn)
        _processStatus(status, "SetVel")

        return

# ======================================================================

    def SetVelD(self, ipmtr, irow, icol, dot):
        """
        ocsm.SetVelD - set the (double) velocity for a Parameter

        inputs:
            ipmtr       Parameter index (1:npmtr)
            irow        row index       (1:nrow)
            icol        column index    (1:ncol)
            dot         velocity to set
        outputs:
            (None)
        """
        _ocsm.ocsmSetVelD.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_double]
        _ocsm.ocsmSetVelD.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetVelD(self._modl, ipmtr, irow, icol, dot)
        _processStatus(status, "SetVelD")

        return

# ======================================================================

    def GetUV(self, ibody, seltype, iselect, npnt, xyz):
        """
        ocsm.GetUV - get the parametric coordinates on an Edge or Face

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.EDGE or ocsm.FACE
            iselect     iedge or iface
            npnt        number of points
            xyz         x[0], y[0], z[0], x[1], ... z[npnt-1]
        outputs:
            uv          u[0], v[0], u[1], ... v[npnt-1]
        """
        _ocsm.ocsmGetUV.argtypes = [ctypes.c_void_p,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.c_int,
                                    ctypes.POINTER(ctypes.c_double),
                                    ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetUV.restype  =  ctypes.c_int

        xyz_ = (ctypes.c_double * (3*npnt))(*xyz)
        uv   = (ctypes.c_double * (2*npnt))()

        status = _ocsm.ocsmGetUV(self._modl, ibody, seltype, iselect, npnt, xyz_, uv)
        _processStatus(status, "GetUV")

        if (seltype == EDGE):
            return list(uv[0:npnt])
        else:
            return list(uv[0:2*npnt])

# ======================================================================

    def GetEgo(self, ibody, seltype, iselect):
        """
        ocsm.GetEgo - get the EGO associated with the Body or its parts

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.BODY, ocsm.NODE, ocsm.EDGE, or ocsm.FACE
            iselect     if ocsm.BODY, 0 for Body
                                      1 for Tessellation
                                      2 for context
                                      3 for EBody
                                      4 for Tessellation on EBody
                        if ocsm.NODE, Node index (1:nnode)
                        if ocsm.EDGE, Node index (1:nedge)
                        if ocsm.FACE, Node index (1:nface)
        outputs:
            theEgo      associated EGO
        """
        _ocsm.ocsmGetEgo.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.POINTER(egads.c_ego)]
        _ocsm.ocsmGetEgo.restype  =  ctypes.c_int

        theEgo = egads.c_ego()
        status = _ocsm.ocsmGetEgo(self._modl, ibody, seltype, iselect, ctypes.byref(theEgo))
        _processStatus(status, "GetEgo")

        # since this is a reference, pyEGADS does not need to manage its memory
        return egads.c_to_py(theEgo, deleteObject=False)

# ======================================================================

    def SetEgo(self, ibody, iselect, theEgo):
        """
        ocsm.SetEgo - set a tessellation or EBody for a Body

        inputs:
            ibody       Body index (1:nbody)
            iselect     1 for Tessellation
                        3 for EBody
                        4 for Tessellation on EBody
            theEgo      associated EGO
        outputs:
            (None)
        """
        _ocsm.ocsmSetEgo.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     egads.c_ego]
        _ocsm.ocsmSetEgo.restype  =  ctypes.c_int

        status = _ocsm.ocsmSetEgo(self._modl, ibody, iselect, theEgo.py_to_c(takeOwnership=True))
        _processStatus(status, "SetEgo")

        return

# ======================================================================

    def FindEnt(self, ibody, seltype, entID):
        """
        ocsm.FindEnt - find entity number given faceID or edgeID */

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.EDGE or ocsm.FACE
            entID       edgeID or faceID
        outputs:
            ient        entity number
        """
        _ocsm.ocsmFindEnt.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmFindEnt.restype  =  ctypes.c_int

        entID_ = (ctypes.c_int * 5)(*entID)
        ient   =  ctypes.c_int()

        status = _ocsm.ocsmFindEnt(self._modl, ibody, seltype, entID_, ient)
        _processStatus(status, "FindEnt")

        return ient.value

# ======================================================================

    def GetXYZ(self, ibody, seltype, iselect, npnt, uv):
        """
        ocsm.GetXYZ - get the coordinates of a Node, Edge, or Face

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.NODE, ocsm.EDGE, or ocsm.FACE
            iselect     inode, iedge, or iface
            npnt        number of points
            uv          u[0], v[0], u[1], ... v[npnt-1]
        outputs:
            xyz         x[0], y[0], z[0], x[1], ... z[npnt-1]
        """
        _ocsm.ocsmGetXYZ.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetXYZ.restype  =  ctypes.c_int

        if   (seltype == NODE):
            npnt = 1
            uv   =  ctypes.c_double
        elif (seltype == EDGE):
            uv_  = (ctypes.c_double * (  npnt))(*uv)
        else:
            uv_  = (ctypes.c_double * (2*npnt))(*uv)
        xyz = (ctypes.c_double * (3*npnt))()

        status = _ocsm.ocsmGetXYZ(self._modl, ibody, seltype, iselect, npnt, uv_, xyz)
        _processStatus(status, "GetXYZ")

        return list(xyz[0:3*npnt])

# ======================================================================

    def GetNorm(self, ibody, iface, npnt, uv):
        """
        ocsm.GetNorm - get the unit normals for a Face

        inputs:
            ibody       Body index (1:nbody)
            iface       Face index (1:nface)
            npnt        number of points
            uv          u[0], v[0], u[1], ... v[npnt-1]
        outputs:
            norm        normx[0], normy[0], normz[0], normx[1], ... normz[npnt-1]
        """
        _ocsm.ocsmGetNorm.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetNorm.restype  =  ctypes.c_int

        uv_  = (ctypes.c_double * (2*npnt))(*uv)
        norm = (ctypes.c_double * (3*npnt))()

        status = _ocsm.ocsmGetNorm(self._modl, ibody, iface, npnt, uv_, norm)
        _processStatus(status, "GetNorm")

        return list(norm[0:3*npnt])

# ======================================================================

    def GetVel(self, ibody, seltype, iselect, npnt, uv):
        """
        ocsm.GetVel - get the velocities on a Node, Edge, or Face

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.NODE, ocsm.EDGE, or ocsm.FACE
            iselect     inode, iedge, or iface
            npnt        number of points
            uv          u[0], v[0], u[1], ... v[npnt-1]
        outputs:
            vel         dx[0], dy[0], dz[0], dx[1], ... dz[npnt-1]
        """
        _ocsm.ocsmGetVel.argtypes = [ctypes.c_void_p,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.c_int,
                                     ctypes.POINTER(ctypes.c_double),
                                     ctypes.POINTER(ctypes.c_double)]
        _ocsm.ocsmGetVel.restype  =  ctypes.c_int

        if   (seltype == NODE):
            npnt = 1
            uv_  = None
        elif (seltype == EDGE):
            if (uv == None):
                uv_ = None
            else:
                uv_ = (ctypes.c_double * (  npnt))(*uv)
        else:
            if (uv == None):
                uv_  = None
            else:
                uv_  = (ctypes.c_double * (2*npnt))(*uv)
        vel = (ctypes.c_double * (3*npnt))()

        status = _ocsm.ocsmGetVel(self._modl, ibody, seltype, iselect, npnt, uv_, vel)
        _processStatus(status, "GetVel")

        return list(vel[0:3*npnt])

# ======================================================================

    def SetEgg(self, eggname):
        """
        ocsm.SetEgg - set up alternative tessellation by an external grid generator

        not implemented
        """

        return

# ======================================================================

    def GetTessNpnt(self, ibody, seltype, iselect):
        """
        ocsm.GetTessNpnt - get the number of tess points in an Edge or Face

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.EDGE or ocsm.FACE
            iselect     Edge or Face index (bias-1)
        outputs:
            npnt        number of points
        """
        _ocsm.ocsmGetTessNpnt.argtypes = [ctypes.c_void_p,
                                          ctypes.c_int,
                                          ctypes.c_int,
                                          ctypes.c_int,
                                          ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmGetTessNpnt.restype  =  ctypes.c_int

        npnt = ctypes.c_int()

        status = _ocsm.ocsmGetTessNpnt(self._modl, ibody, seltype, iselect, npnt)
        _processStatus(status, "GetTessNpnt")

        return npnt.value

# ======================================================================

    def GetTessVel(self, ibody, seltype, iselect):
        """
        ocsm.GetTessVel - get the tessellation velocities on a Node, Edge, or Face

        inputs:
            ibody       Body index (1:nbody)
            seltype     ocsm.NODE, ocsm.EDGE, or ocsm.FACE
            iselect     Node, Edge, or Face index (bias-1)
        outputs:
            dxyz        velocities
        """
        _ocsm.ocsmGetTessVel.argtypes = [ctypes.c_void_p,
                                         ctypes.c_int,
                                         ctypes.c_int,
                                         ctypes.c_int,
                                         ctypes.POINTER(ctypes.POINTER(ctypes.c_double))]
        _ocsm.ocsmGetTessVel.restype  =  ctypes.c_int

        if   (seltype == NODE):
            npnt = 1
        elif (seltype == EDGE):
            npnt = self.GetTessNpnt(ibody, seltype, iselect)
        elif (seltype == FACE):
            npnt = self.GetTessNpnt(ibody, seltype, iselect)
        else:
            npnt = 1

        dxyz = ctypes.POINTER(ctypes.c_double)()

        status = _ocsm.ocsmGetTessVel(self._modl, ibody, seltype, iselect, ctypes.byref(dxyz))
        _processStatus(status, "GetTessVel")

        return list(dxyz[0:3*npnt])

# ======================================================================

    def GetBody(self, ibody):
        """
        ocsm.GetBody - get info about a Body

        inputs:
            ibody       Body index (1:nbody)
        outputs:
            type        Branch type
            ichild      ibody of child (or 0 if root)
            ileft       ibody of left parent (or 0)
            irite       ibody of rite parent (or 0)
            vals        array of Argument values
            nnode       number of Nodes
            nedge       number of Edges
            nface       number of Faces
        """
        _ocsm.ocsmGetBody.argtypes = [ctypes.c_void_p,
                                      ctypes.c_int,
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_double),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int),
                                      ctypes.POINTER(ctypes.c_int)]
        _ocsm.ocsmGetBody.restype  =  ctypes.c_int

        type  = ctypes.c_int()
        ichld = ctypes.c_int()
        ileft = ctypes.c_int()
        irite = ctypes.c_int()
        vals  = (ctypes.c_double * 10)()
        nnode = ctypes.c_int()
        nedge = ctypes.c_int()
        nface = ctypes.c_int()

        status = _ocsm.ocsmGetBody(self._modl, ibody, type, ichld, ileft, irite, vals, nnode, nedge, nface)
        _processStatus(status, "GetBody")

        outvals = list(vals[0:10])

        return (type.value, ichld.value, ileft.value, irite.value, outvals, nnode.value, nedge.value, nface.value)

# ======================================================================

    def PrintBodys(self, filename):
        """
        ocsm.PrintBodys - print all Bodys to file

        inputs:
            filename    file to which output is appended (or "" for stdout)
        outputs:
            (None)
        """
        _ocsm.ocsmPrintBodys.argtypes = [ctypes.c_void_p,
                                         ctypes.c_char_p]
        _ocsm.ocsmPrintBodys.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmPrintBodys(self._modl, filename)
        _processStatus(status, "PrintBodys")

        return

# ======================================================================

    def PrintBrep(self, ibody, filename):
        """
        ocsm.PrintBrep - print the BRep associated with a specific Body

        inputs:
            ibody       Body index (1-nbody)
            filename    file to which output is appended (or "" for stdout)
        outputs:
            (None)
        """
        _ocsm.ocsmPrintBrep.argtypes = [ctypes.c_void_p,
                                        ctypes.c_int,
                                        ctypes.c_char_p]
        _ocsm.ocsmPrintBrep.restype  =  ctypes.c_int

        if (isinstance(filename, str)):
            filename = filename.encode()

        status = _ocsm.ocsmPrintBrep(self._modl, ibody, filename)
        _processStatus(status, "PrintBrep")

        return

# ======================================================================

    def EvalExpr(self, express):
        """
        ocsm.EvalExpr - evaluate an expression

        inputs:
            express     expression
        outputs:
            value       value
            dot         velocity
            string      value if string-valued
        """
        _ocsm.ocsmEvalExpr.argtypes = [ctypes.c_void_p,
                                       ctypes.c_char_p,
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.POINTER(ctypes.c_double),
                                       ctypes.c_char_p]
        _ocsm.ocsmEvalExpr.restype  =  ctypes.c_int

        if (isinstance(express, str)):
            express = express.encode()

        value  = ctypes.c_double()
        dot    = ctypes.c_double()
        string = ctypes.create_string_buffer(257)

        status = _ocsm.ocsmEvalExpr(self._modl, express, value, dot, string)
        _processStatus(status, "EvalExpr")

        return (value.value, dot.value, string.value.decode("utf-8"))

# ======================================================================

    def GetText(self, code):
        """
        ocsm.GetText - convert an OCSM numeric code to text

        inputs:
            code        code to look up
        outputs:
            text        text string
        """
        _ocsm.ocsmGetText.argtypes = [ctypes.c_int]
        _ocsm.ocsmGetText.restype  =  ctypes.c_char_p

        text = _ocsm.ocsmGetText(code)

        return text.decode("utf-8")

# ======================================================================

    def GetCode(self, text):
        """
        ocsm.GetCode - convert text to an OCSM numeric code

        inputs:
            text        text string
        outputs:
            code        numeric code
        """
        _ocsm.ocsmGetCode.argtypes = [ctypes.c_char_p]
        _ocsm.ocsmGetCode.restype  =  ctypes.c_int

        if (isinstance(text, str)):
            text = text.encode()

        code = _ocsm.ocsmGetCode(text)

        return code

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

    # this list comes from OpenCSM.h
    if   (status == -201):
        raise OcsmError("FILE_NOT_FOUND",              routine)
    elif (status == -202):
        raise OcsmError("ILLEGAL_STATEMENT",           routine)
    elif (status == -203):
        raise OcsmError("NOT_ENOUGH_ARGS",             routine)
    elif (status == -204):
        raise OcsmError("NAME_ALREADY_DEFINED",        routine)
    elif (status == -205):
        raise OcsmError("NESTED_TOO_DEEPLY",           routine)
    elif (status == -206):
        raise OcsmError("IMPROPER_NESTING",            routine)
    elif (status == -207):
        raise OcsmError("NESTING_NOT_CLOSED",          routine)
    elif (status == -208):
        raise OcsmError("NOT_MODL_STRUCTURE",          routine)
    elif (status == -209):
        raise OcsmError("PROBLEM_CREATING_PERTURB",    routine)
    elif (status == -211):
        raise OcsmError("MISSING_MARK",                routine)
    elif (status == -212):
        raise OcsmError("INSUFFICIENT_BODYS_ON_STACK", routine)
    elif (status == -213):
        raise OcsmError("WRONG_TYPES_ON_STACK",        routine)
    elif (status == -214):
        raise OcsmError("DID_NOT_CREATE_BODY",         routine)
    elif (status == -215):
        raise OcsmError("CREATED_TOO_MANY_BODYS",      routine)
    elif (status == -216):
        raise OcsmError("TOO_MANY_BODYS_ON_STACK",     routine)
    elif (status == -217):
        raise OcsmError("ERROR_IN_BODYS_ON_STACK",     routine)
    elif (status == -218):
        raise OcsmError("MODL_NOT_CHECKED",            routine)
    elif (status == -219):
        raise OcsmError("NEED_TESSELLATION",           routine)
    elif (status == -221):
        raise OcsmError("BODY_NOT_FOUND",              routine)
    elif (status == -222):
        raise OcsmError("FACE_NOT_FOUND",              routine)
    elif (status == -223):
        raise OcsmError("EDGE_NOT_FOUND",              routine)
    elif (status == -224):
        raise OcsmError("NODE_NOT_FOUND",              routine)
    elif (status == -225):
        raise OcsmError("ILLEGAL_VALUE",               routine)
    elif (status == -226):
        raise OcsmError("ILLEGAL_ATTRIBUTE",           routine)
    elif (status == -227):
        raise OcsmError("ILLEGAL_CSYSTEM",             routine)
    elif (status == -228):
        raise OcsmError("NO_SELECTION",                routine)
    elif (status == -231):
        raise OcsmError("SKETCH_IS_OPEN",              routine)
    elif (status == -232):
        raise OcsmError("SKETCH_IS_NOT_OPEN",          routine)
    elif (status == -233):
        raise OcsmError("COLINEAR_SKETCH_POINTS",      routine)
    elif (status == -234):
        raise OcsmError("NON_COPLANAR_SKETCH_POINTS",  routine)
    elif (status == -235):
        raise OcsmError("TOO_MANY_SKETCH_POINTS",      routine)
    elif (status == -236):
        raise OcsmError("TOO_FEW_SPLINE_POINTS",       routine)
    elif (status == -237):
        raise OcsmError("SKETCH_DOES_NOT_CLOSE",       routine)
    elif (status == -238):
        raise OcsmError("SELF_INTERSECTING",           routine)
    elif (status == -239):
        raise OcsmError("ASSERT_FAILED",               routine)
    elif (status == -241):
        raise OcsmError("ILLEGAL_CHAR_IN_EXPR",        routine)
    elif (status == -242):
        raise OcsmError("CLOSE_BEFORE_OPEN",           routine)
    elif (status == -243):
        raise OcsmError("MISSING_CLOSE",               routine)
    elif (status == -244):
        raise OcsmError("ILLEGAL_TOKEN_SEQUENCE",      routine)
    elif (status == -245):
        raise OcsmError("ILLEGAL_NUMBER",              routine)
    elif (status == -246):
        raise OcsmError("ILLEGAL_PMTR_NAME",           routine)
    elif (status == -247):
        raise OcsmError("ILLEGAL_FUNC_NAME",           routine)
    elif (status == -248):
        raise OcsmError("ILLEGAL_TYPE",                routine)
    elif (status == -249):
        raise OcsmError("ILLEGAL_NARG",                routine)
    elif (status == -251):
        raise OcsmError("NAME_NOT_FOUND",              routine)
    elif (status == -252):
        raise OcsmError("NAME_NOT_UNIQUE",             routine)
    elif (status == -253):
        raise OcsmError("PMTR_IS_DESPMTR",             routine)
    elif (status == -254):
        raise OcsmError("PMTR_IS_LOCALVAR",            routine)
    elif (status == -255):
        raise OcsmError("PMTR_IS_OUTPMTR",             routine)
    elif (status == -256):
        raise OcsmError("PMTR_IS_CONPMTR",             routine)
    elif (status == -257):
        raise OcsmError("WRONG_PMTR_TYPE",             routine)
    elif (status == -258):
        raise OcsmError("FUNC_ARG_OUT_OF_BOUNDS",      routine)
    elif (status == -259):
        raise OcsmError("VAL_STACK_UNDERFLOW",         routine)
    elif (status == -260):
        raise OcsmError("VAL_STACK_OVERFLOW",          routine)
    elif (status == -261):
        raise OcsmError("ILLEGAL_BRCH_INDEX",          routine)
    elif (status == -262):
        raise OcsmError("ILLEGAL_PMTR_INDEX",          routine)
    elif (status == -263):
        raise OcsmError("ILLEGAL_BODY_INDEX",          routine)
    elif (status == -264):
        raise OcsmError("ILLEGAL_ARG_INDEX",           routine)
    elif (status == -265):
        raise OcsmError("ILLEGAL_ACTIVITY",            routine)
    elif (status == -266):
        raise OcsmError("ILLEGAL_MACRO_INDEX",         routine)
    elif (status == -267):
        raise OcsmError("ILLEGAL_ARGUMENT",            routine)
    elif (status == -268):
        raise OcsmError("CANNOT_BE_SUPPRESSED",        routine)
    elif (status == -269):
        raise OcsmError("STORAGE_ALREADY_USED",        routine)
    elif (status == -270):
        raise OcsmError("NOTHING_PREVIOUSLY_STORED",   routine)
    elif (status == -271):
        raise OcsmError("SOLVER_IS_OPEN",              routine)
    elif (status == -272):
        raise OcsmError("SOLVER_IS_NOT_OPEN",          routine)
    elif (status == -273):
        raise OcsmError("TOO_MANY_SOLVER_VARS",        routine)
    elif (status == -274):
        raise OcsmError("UNDERCONSTRAINED",            routine)
    elif (status == -275):
        raise OcsmError("OVERCONSTRAINED",             routine)
    elif (status == -276):
        raise OcsmError("SINGULAR_MATRIX",             routine)
    elif (status == -277):
        raise OcsmError("NOT_CONVERGED",               routine)
    elif (status == -281):
        raise OcsmError("UDP_ERROR1",                  routine)
    elif (status == -282):
        raise OcsmError("UDP_ERROR2",                  routine)
    elif (status == -283):
        raise OcsmError("UDP_ERROR3",                  routine)
    elif (status == -284):
        raise OcsmError("UDP_ERROR4",                  routine)
    elif (status == -285):
        raise OcsmError("UDP_ERROR5",                  routine)
    elif (status == -286):
        raise OcsmError("UDP_ERROR6",                  routine)
    elif (status == -287):
        raise OcsmError("UDP_ERROR7",                  routine)
    elif (status == -288):
        raise OcsmError("UDP_ERROR8",                  routine)
    elif (status == -289):
        raise OcsmError("UDP_ERROR9",                  routine)
    elif (status == -291):
        raise OcsmError("OP_STACK_UNDERFLOW",          routine)
    elif (status == -292):
        raise OcsmError("OP_STACK_OVERFLOW",           routine)
    elif (status == -293):
        raise OcsmError("RPN_STACK_UNDERFLOW",         routine)
    elif (status == -294):
        raise OcsmError("RPN_STACK_OVERFLOW",          routine)
    elif (status == -295):
        raise OcsmError("TOKEN_STACK_UNDERFLOW",       routine)
    elif (status == -296):
        raise OcsmError("TOKEN_STACK_OVERFLOW",        routine)
    elif (status == -298):
        raise OcsmError("UNSUPPORTED",                 routine)
    elif (status == -299):
        raise OcsmError("INTERNAL_ERROR",              routine)

# ======================================================================

# cleanup udp storage by calling _ocsm.ocsmFree(None) at exit

atexit.register(_ocsmFree)

# ======================================================================
# OcsmError class
# ======================================================================

class OcsmError(Exception):

    # Constructor or Initializer
    def __init__(self, value, routine):
        self.value   = value
        self.routine = routine

    # __str__ is to print() the value
    def __str__(self):
        return("OcsmError: "+repr(self.value)+" detected during call to "+self.routine)

# ======================================================================

# print version number
#(imajor, iminor) = Version()
#print("\nocsm version {0:d}.{1:d} has been loaded\n".format(imajor,iminor))
#del imajor
#del iminor
