#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

## CAPS error exception class. 
# See \ref exception.py for a representative use case.
class CAPSError(Exception):
    ## \example exception.py
    # Example of raised error (pyCAPS.CAPSError) handling in pyCAPS.
    
    # Check with cythonDoxCleanup.py that the correct error codes are being replaced here!
    
    ## \showinitializer Dictionary of CAPS errors {errorCode : errorName}.
    capsError = { 
 cCAPS.CAPS_SUCCESS      : 'CAPS_SUCCESS',
 cCAPS.CAPS_BADRANK      : 'CAPS_BADRANK',
 cCAPS.CAPS_BADDSETNAME  : 'CAPS_BADDSETNAME',
 cCAPS.CAPS_NOTFOUND     : 'CAPS_NOTFOUND',
 cCAPS.CAPS_BADINDEX     : 'CAPS_BADINDEX',
 cCAPS.CAPS_NOTCHANGED   : 'CAPS_NOTCHANGED',
 cCAPS.CAPS_BADTYPE      : 'CAPS_BADTYPE',
 cCAPS.CAPS_NULLVALUE    : 'CAPS_NULLVALUE',
 cCAPS.CAPS_NULLNAME     : 'CAPS_NULLNAME',
 cCAPS.CAPS_NULLOBJ      : 'CAPS_NULLOBJ',
 cCAPS.CAPS_BADOBJECT    : 'CAPS_BADOBJECT',
 cCAPS.CAPS_BADVALUE     : 'CAPS_BADVALUE',
 cCAPS.CAPS_PARAMBNDERR  : 'CAPS_PARAMBNDERR',
 cCAPS.CAPS_NOTCONNECT   : 'CAPS_NOTCONNECT',
 cCAPS.CAPS_NOTPARMTRIC  : 'CAPS_NOTPARMTRIC',
 cCAPS.CAPS_READONLYERR  : 'CAPS_READONLYERR',
 cCAPS.CAPS_FIXEDLEN     : 'CAPS_FIXEDLEN',
 cCAPS.CAPS_BADNAME      : 'CAPS_BADNAME',
 cCAPS.CAPS_BADMETHOD    : 'CAPS_BADMETHOD',
 cCAPS.CAPS_CIRCULARLINK : 'CAPS_CIRCULARLINK',
 cCAPS.CAPS_UNITERR      : 'CAPS_UNITERR',
 cCAPS.CAPS_NULLBLIND    : 'CAPS_NULLBLIND',
 cCAPS.CAPS_SHAPEERR     : 'CAPS_SHAPEERR',
 cCAPS.CAPS_LINKERR      : 'CAPS_LINKERR',
 cCAPS.CAPS_MISMATCH     : 'CAPS_MISMATCH',
 cCAPS.CAPS_NOTPROBLEM   : 'CAPS_NOTPROBLEM',
 cCAPS.CAPS_RANGEERR     : 'CAPS_RANGEERR',
 cCAPS.CAPS_DIRTY        : 'CAPS_DIRTY',
 cCAPS.CAPS_HIERARCHERR  : 'CAPS_HIERARCHERR',
 cCAPS.CAPS_STATEERR     : 'CAPS_STATEERR',
 cCAPS.CAPS_SOURCEERR    : 'CAPS_SOURCEERR',
 cCAPS.CAPS_EXISTS       : 'CAPS_EXISTS',
 cCAPS.CAPS_IOERR        : 'CAPS_IOERR',
 cCAPS.CAPS_DIRERR       : 'CAPS_DIRERR',
 cCAPS.CAPS_NOTIMPLEMENT : 'CAPS_NOTIMPLEMENT',
 cCAPS.CAPS_EXECERR      : 'CAPS_EXECERR',
 cCAPS.CAPS_CLEAN        : 'CAPS_CLEAN',
 cCAPS.CAPS_BADINTENT    : 'CAPS_BADINTENT',
 cCAPS.CAPS_BADINIT      : 'CAPS_BADINIT',   
 cCAPS.CAPS_NOTNEEDED    : 'CAPS_NOTNEEDED',
 cCAPS.CAPS_NOSENSITVTY  : 'CAPS_NOSENSITVTY'}
        
    ## \showinitializer Dictionary of EGADS errors {errorCode : errorName}.
    egadsError = {
   cEGADS.EGADS_CNTXTHRD : 'EGADS_CNTXTHRD',
   cEGADS.EGADS_READERR  : 'EGADS_READERR',
   cEGADS.EGADS_TESSTATE : 'EGADS_TESSTATE',
   cEGADS.EGADS_EXISTS   : 'EGADS_EXISTS',
   cEGADS.EGADS_ATTRERR  : 'EGADS_ATTRERR',
   cEGADS.EGADS_TOPOCNT  : 'EGADS_TOPOCNT',
   cEGADS.EGADS_BADSCALE : 'EGADS_BADSCALE',
   cEGADS.EGADS_NOTORTHO : 'EGADS_NOTORTHO',
   cEGADS.EGADS_DEGEN    : 'EGADS_DEGEN',
   cEGADS.EGADS_CONSTERR : 'EGADS_CONSTERR',
   cEGADS.EGADS_TOPOERR  : 'EGADS_TOPOERR',
   cEGADS.EGADS_GEOMERR  : 'EGADS_GEOMERR',
   cEGADS.EGADS_NOTBODY  : 'EGADS_NOTBODY',
   cEGADS.EGADS_WRITERR  : 'EGADS_WRITERR',
   cEGADS.EGADS_NOTMODEL : 'EGADS_NOTMODEL',
   cEGADS.EGADS_NOLOAD   : 'EGADS_NOLOAD',
   cEGADS.EGADS_RANGERR  : 'EGADS_RANGERR',
   cEGADS.EGADS_NOTGEOM  : 'EGADS_NOTGEOM',
   cEGADS.EGADS_NOTTESS  : 'EGADS_NOTTESS',
   cEGADS.EGADS_EMPTY    : 'EGADS_EMPTY',
   cEGADS.EGADS_NOTTOPO  : 'EGADS_NOTTOPO',
   cEGADS.EGADS_REFERCE  : 'EGADS_REFERCE',
   cEGADS.EGADS_NOTXFORM : 'EGADS_NOTXFORM',
   cEGADS.EGADS_NOTCNTX  : 'EGADS_NOTCNTX',
   cEGADS.EGADS_MIXCNTX  : 'EGADS_MIXCNTX',
   cEGADS.EGADS_NODATA   : 'EGADS_NODATA',
   cEGADS.EGADS_NONAME   : 'EGADS_NONAME',
   cEGADS.EGADS_INDEXERR : 'EGADS_INDEXERR',
   cEGADS.EGADS_MALLOC   : 'EGADS_MALLOC',
   cEGADS.EGADS_NOTOBJ   : 'EGADS_NOTOBJ',
   cEGADS.EGADS_NULLOBJ  : 'EGADS_NULLOBJ',
   cEGADS.EGADS_NOTFOUND : 'EGADS_NOTFOUND', 
   cEGADS.EGADS_SUCCESS  : 'EGADS_SUCCESS', 
   cEGADS.EGADS_OUTSIDE  : 'EGADS_OUTSIDE',
   
   cEGADS.PRM_NOTCONVERGED   : 'PRM_NOTCONVERGED',
   cEGADS.PRM_BADNUMVERTICES : 'PRM_BADNUMVERTICES',
   cEGADS.PRM_ZEROPIVOT      : 'PRM_ZEROPIVOT',
   cEGADS.PRM_NEGATIVEAREAS  : 'PRM_NEGATIVEAREAS',
   cEGADS.PRM_TOLERANCEUNMET : 'PRM_TOLERANCEUNMET',
   cEGADS.PRM_NOGLOBALUV     : 'PRM_NOGLOBALUV',
   cEGADS.PRM_BADPARAM       : 'PRM_BADPARAM',
   cEGADS.PRM_BADDIVISION    : 'PRM_BADDIVISION',
   cEGADS.PRM_CANNOTFORMLOOP : 'PRM_CANNOTFORMLOOP',
   cEGADS.PRM_WIGGLEDETECTED : 'PRM_WIGGLEDETECTED',
   cEGADS.PRM_INTERNAL       : 'PRM_INTERNAL' }
    
    ## \showinitializer Dictionary of OCSM errors {errorCode : errorName}.
    ocsmError = {cOCSM.OCSM_FILE_NOT_FOUND : 'OCSM_FILE_NOT_FOUND',
    cOCSM.OCSM_ILLEGAL_STATEMENT     : 'OCSM_ILLEGAL_STATEMENT',
    cOCSM.OCSM_NOT_ENOUGH_ARGS       : 'OCSM_NOT_ENOUGH_ARGS',
    cOCSM.OCSM_NAME_ALREADY_DEFINED  : 'OCSM_NAME_ALREADY_DEFINED',
    cOCSM.OCSM_NESTED_TOO_DEEPLY     : 'OCSM_NESTED_TOO_DEEPLY',
    cOCSM.OCSM_IMPROPER_NESTING      : 'OCSM_IMPROPER_NESTING',
    cOCSM.OCSM_NESTING_NOT_CLOSED    : 'OCSM_NESTING_NOT_CLOSED',    
    cOCSM.OCSM_NOT_MODL_STRUCTURE    : 'OCSM_NOT_MODL_STRUCTURE',     
    cOCSM.OCSM_PROBLEM_CREATING_PERTURB  : 'OCSM_PROBLEM_CREATING_PERTURB',    
    
    cOCSM.OCSM_MISSING_MARK                : 'OCSM_MISSING_MARK',          
    cOCSM.OCSM_INSUFFICIENT_BODYS_ON_STACK : 'OCSM_INSUFFICIENT_BODYS_ON_STACK',
    cOCSM.OCSM_WRONG_TYPES_ON_STACK        : 'OCSM_WRONG_TYPES_ON_STACK',   
    cOCSM.OCSM_DID_NOT_CREATE_BODY         : 'OCSM_DID_NOT_CREATE_BODY',         
    cOCSM.OCSM_CREATED_TOO_MANY_BODYS      : 'OCSM_CREATED_TOO_MANY_BODYS',      
    cOCSM.OCSM_TOO_MANY_BODYS_ON_STACK     : 'OCSM_TOO_MANY_BODYS_ON_STACK',     
    cOCSM.OCSM_ERROR_IN_BODYS_ON_STACK     : 'OCSM_ERROR_IN_BODYS_ON_STACK',       
    cOCSM.OCSM_MODL_NOT_CHECKED            : 'OCSM_MODL_NOT_CHECKED',       
    cOCSM.OCSM_NEED_TESSELLATION           : 'OCSM_NEED_TESSELLATION',       
    
    cOCSM.OCSM_BODY_NOT_FOUND     : 'OCSM_BODY_NOT_FOUND',           
    cOCSM.OCSM_FACE_NOT_FOUND     : 'OCSM_FACE_NOT_FOUND',          
    cOCSM.OCSM_EDGE_NOT_FOUND     : 'OCSM_EDGE_NOT_FOUND',             
    cOCSM.OCSM_NODE_NOT_FOUND     : 'OCSM_NODE_NOT_FOUND',
    cOCSM.OCSM_ILLEGAL_VALUE      : 'OCSM_ILLEGAL_VALUE',             
    cOCSM.OCSM_ILLEGAL_ATTRIBUTE  : 'OCSM_ILLEGAL_ATTRIBUTE',            
    cOCSM.OCSM_ILLEGAL_CSYSTEM    : 'OCSM_ILLEGAL_CSYSTEM',             
    
    cOCSM.OCSM_SKETCH_IS_OPEN             : 'OCSM_SKETCH_IS_OPEN',             
    cOCSM.OCSM_SKETCH_IS_NOT_OPEN         : 'OCSM_SKETCH_IS_NOT_OPEN',        
    cOCSM.OCSM_COLINEAR_SKETCH_POINTS     : 'OCSM_COLINEAR_SKETCH_POINTS',    
    cOCSM.OCSM_NON_COPLANAR_SKETCH_POINTS : 'OCSM_NON_COPLANAR_SKETCH_POINTS',   
    cOCSM.OCSM_TOO_MANY_SKETCH_POINTS     : 'OCSM_TOO_MANY_SKETCH_POINTS',   
    cOCSM.OCSM_TOO_FEW_SPLINE_POINTS      : 'OCSM_TOO_FEW_SPLINE_POINTS',  
    cOCSM.OCSM_SKETCH_DOES_NOT_CLOSE      : 'OCSM_SKETCH_DOES_NOT_CLOSE',        
    cOCSM.OCSM_SELF_INTERSECTING          : 'OCSM_SELF_INTERSECTING',           
    cOCSM.OCSM_ASSERT_FAILED              : 'OCSM_ASSERT_FAILED',     
    
    cOCSM.OCSM_ILLEGAL_CHAR_IN_EXPR       : 'OCSM_ILLEGAL_CHAR_IN_EXPR',        
    cOCSM.OCSM_CLOSE_BEFORE_OPEN          : 'OCSM_CLOSE_BEFORE_OPEN',    
    cOCSM.OCSM_MISSING_CLOSE              : 'OCSM_MISSING_CLOSE',    
    cOCSM.OCSM_ILLEGAL_TOKEN_SEQUENCE     : 'OCSM_ILLEGAL_TOKEN_SEQUENCE',  
    cOCSM.OCSM_ILLEGAL_NUMBER             : 'OCSM_ILLEGAL_NUMBER',   
    cOCSM.OCSM_ILLEGAL_PMTR_NAME          : 'OCSM_ILLEGAL_PMTR_NAME',   
    cOCSM.OCSM_ILLEGAL_FUNC_NAME          : 'OCSM_ILLEGAL_FUNC_NAME',   
    cOCSM.OCSM_ILLEGAL_TYPE               : 'OCSM_ILLEGAL_TYPE',  
    cOCSM.OCSM_ILLEGAL_NARG               : 'OCSM_ILLEGAL_NARG',    
    
    cOCSM.OCSM_NAME_NOT_FOUND             : 'OCSM_NAME_NOT_FOUND',    
    cOCSM.OCSM_NAME_NOT_UNIQUE            : 'OCSM_NAME_NOT_UNIQUE',    
    #cOCSM.OCSM_PMTR_IS_EXTERNAL           : 'OCSM_PMTR_IS_EXTERNAL',    
    #cOCSM.OCSM_PMTR_IS_INTERNAL           : 'OCSM_PMTR_IS_INTERNAL',    
    #cOCSM.OCSM_PMTR_IS_CONSTANT           : 'OCSM_PMTR_IS_CONSTANT',    
    cOCSM.OCSM_WRONG_PMTR_TYPE            : 'OCSM_WRONG_PMTR_TYPE',    
    cOCSM.OCSM_FUNC_ARG_OUT_OF_BOUNDS     : 'OCSM_FUNC_ARG_OUT_OF_BOUNDS',    
    cOCSM.OCSM_VAL_STACK_UNDERFLOW        : 'OCSM_VAL_STACK_UNDERFLOW',     
    cOCSM.OCSM_VAL_STACK_OVERFLOW         : 'OCSM_VAL_STACK_OVERFLOW',     
    
    cOCSM.OCSM_ILLEGAL_BRCH_INDEX         : 'OCSM_ILLEGAL_BRCH_INDEX',     
    cOCSM.OCSM_ILLEGAL_PMTR_INDEX         : 'OCSM_ILLEGAL_PMTR_INDEX',    
    cOCSM.OCSM_ILLEGAL_BODY_INDEX         : 'OCSM_ILLEGAL_BODY_INDEX',     
    cOCSM.OCSM_ILLEGAL_ARG_INDEX          : 'OCSM_ILLEGAL_ARG_INDEX',   
    cOCSM.OCSM_ILLEGAL_ACTIVITY           : 'OCSM_ILLEGAL_ACTIVITY',     
    cOCSM.OCSM_ILLEGAL_MACRO_INDEX        : 'OCSM_ILLEGAL_MACRO_INDEX',     
    cOCSM.OCSM_ILLEGAL_ARGUMENT           : 'OCSM_ILLEGAL_ARGUMENT',     
    cOCSM.OCSM_CANNOT_BE_SUPPRESSED       : 'OCSM_CANNOT_BE_SUPPRESSED',    
    cOCSM.OCSM_STORAGE_ALREADY_USED       : 'OCSM_STORAGE_ALREADY_USED',    
    cOCSM.OCSM_NOTHING_PREVIOUSLY_STORED  : 'OCSM_NOTHING_PREVIOUSLY_STORED',    
    
    cOCSM.OCSM_SOLVER_IS_OPEN             : 'OCSM_SOLVER_IS_OPEN',    
    cOCSM.OCSM_SOLVER_IS_NOT_OPEN         : 'OCSM_SOLVER_IS_NOT_OPEN',    
    cOCSM.OCSM_TOO_MANY_SOLVER_VARS       : 'OCSM_TOO_MANY_SOLVER_VARS',     
    cOCSM.OCSM_UNDERCONSTRAINED           : 'OCSM_UNDERCONSTRAINED',     
    cOCSM.OCSM_OVERCONSTRAINED            : 'OCSM_OVERCONSTRAINED',     
    cOCSM.OCSM_SINGULAR_MATRIX            : 'OCSM_SINGULAR_MATRIX',     
    cOCSM.OCSM_NOT_CONVERGED              : 'OCSM_NOT_CONVERGED',     
    
    cOCSM.OCSM_UDP_ERROR1                 : 'OCSM_UDP_ERROR1',     
    cOCSM.OCSM_UDP_ERROR2                 : 'OCSM_UDP_ERROR2',     
    cOCSM.OCSM_UDP_ERROR3                 : 'OCSM_UDP_ERROR3',     
    cOCSM.OCSM_UDP_ERROR4                 : 'OCSM_UDP_ERROR4',     
    cOCSM.OCSM_UDP_ERROR5                 : 'OCSM_UDP_ERROR5',     
    cOCSM.OCSM_UDP_ERROR6                 : 'OCSM_UDP_ERROR6',     
    cOCSM.OCSM_UDP_ERROR7                 : 'OCSM_UDP_ERROR7',     
    cOCSM.OCSM_UDP_ERROR8                 : 'OCSM_UDP_ERROR8',     
    cOCSM.OCSM_UDP_ERROR9                 : 'OCSM_UDP_ERROR9',     
    
    cOCSM.OCSM_OP_STACK_UNDERFLOW         : 'OCSM_OP_STACK_UNDERFLOW',     
    cOCSM.OCSM_OP_STACK_OVERFLOW          : "OCSM_OP_STACK_OVERFLOW",     
    cOCSM.OCSM_RPN_STACK_UNDERFLOW        : 'OCSM_RPN_STACK_UNDERFLOW',     
    cOCSM.OCSM_RPN_STACK_OVERFLOW         : 'OCSM_RPN_STACK_OVERFLOW',     
    cOCSM.OCSM_TOKEN_STACK_UNDERFLOW      : 'OCSM_TOKEN_STACK_UNDERFLOW',    
    cOCSM.OCSM_TOKEN_STACK_OVERFLOW       : 'OCSM_TOKEN_STACK_OVERFLOW',    
    cOCSM.OCSM_UNSUPPORTED                : 'OCSM_UNSUPPORTED',    
    cOCSM.OCSM_INTERNAL_ERROR             : 'OCSM_INTERNAL_ERROR'}
    
    ## Initialize the CAPSError exception. 
    def __init__(self, code=None, msg=None):
        
        # Combine errors 
        self.errorDict = self.capsError.copy()
        self.errorDict.update(self.egadsError)
        self.errorDict.update(self.ocsmError)
        
        if code is None:
            code = -999
        if msg is None:
            msg = "CAPS error detected!!!"
        
        ## Error code encountered when running the CAPS
        self.errorCode = code
        
        ## Name of error code encountered
        self.errorName = "Unknown_Error"
        
        if code not in self.errorDict:
            msgPre = "CAPS error detected: error code " + str(code) + " = " + " undefined error, "
        else:
            msgPre = "CAPS error detected: error code " + str(code) + " = " + \
                                                          self.errorDict[code] + " received, "
            self.errorName = self.errorDict[code]
        
        
        super(CAPSError,self).__init__(msgPre + msg)
        
        
# Not currently used
# class UNKNOWN_ERROR(CAPSError):
#     def __init__(self, msg=None):
#         super(UNKNOWN_ERROR, self).__init__(msg)
# 
# class EGADS_SUCCESS(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_SUCCESS, self).__init__(msg)
# 
# class EGADS_OUTSIDE(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_OUTSIDE, self).__init__(msg)
# 
# class EGADS_NOTFOUND(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTFOUND, self).__init__(msg)
# 
# class CAPS_BADINTENT(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADINTENT, self).__init__(msg)
# 
# class CAPS_CLEAN(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_CLEAN, self).__init__(msg)
# 
# class CAPS_EXECERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_EXECERR, self).__init__(msg)
# 
# class CAPS_NOTIMPLEMENT(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTIMPLEMENT, self).__init__(msg)
# 
# class CAPS_DIRERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_DIRERR, self).__init__(msg)
# 
# class CAPS_IOERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_IOERR, self).__init__(msg)
# 
# class CAPS_EXISTS(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_EXISTS, self).__init__(msg)
# 
# class CAPS_SOURCEERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_SOURCEERR, self).__init__(msg)
# 
# class CAPS_STATEERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_STATEERR, self).__init__(msg)
# 
# class CAPS_HIERARCHERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_HIERARCHERR, self).__init__(msg)
# 
# class CAPS_DIRTY(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_DIRTY, self).__init__(msg)
# 
# class CAPS_RANGEERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_RANGEERR, self).__init__(msg)
# 
# class CAPS_NOTPROBLEM(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTPROBLEM, self).__init__(msg)
# 
# class CAPS_MISMATCH(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_MISMATCH, self).__init__(msg)
# 
# class CAPS_LINKERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_LINKERR, self).__init__(msg)
# 
# class CAPS_SHAPEERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_SHAPEERR, self).__init__(msg)
# 
# class CAPS_NULLBLIND(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NULLBLIND, self).__init__(msg)
# 
# class CAPS_UNITERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_UNITERR, self).__init__(msg)
# 
# class CAPS_CIRCULARLINK(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_CIRCULARLINK, self).__init__(msg)
# 
# class CAPS_BADMETHOD(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADMETHOD, self).__init__(msg)
# 
# class CAPS_BADNAME(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADNAME, self).__init__(msg)
# 
# class CAPS_FIXEDLEN(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_FIXEDLEN, self).__init__(msg)
# 
# class CAPS_READONLYERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_READONLYERR, self).__init__(msg)
# 
# class CAPS_NOTPARMTRIC(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTPARMTRIC, self).__init__(msg)
# 
# class CAPS_NOTCONNECT(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTCONNECT, self).__init__(msg)
# 
# class CAPS_PARAMBNDERR(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_PARAMBNDERR, self).__init__(msg)
# 
# class CAPS_BADVALUE(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADVALUE, self).__init__(msg)
# 
# class CAPS_BADOBJECT(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADOBJECT, self).__init__(msg)
# 
# class CAPS_NULLOBJ(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NULLOBJ, self).__init__(msg)
# 
# class CAPS_NULLNAME(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NULLNAME, self).__init__(msg)
# 
# class CAPS_NULLVALUE(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NULLVALUE, self).__init__(msg)
# 
# class CAPS_BADTYPE(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADTYPE, self).__init__(msg)
# 
# class CAPS_NOTCHANGED(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTCHANGED, self).__init__(msg)
# 
# class CAPS_BADINDEX(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADINDEX, self).__init__(msg)
# 
# class CAPS_NOTFOUND(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_NOTFOUND, self).__init__(msg)
# 
# class CAPS_BADDSETNAME(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADDSETNAME, self).__init__(msg)
# 
# class CAPS_BADRANK(CAPSError):
#     def __init__(self, msg=None):
#         super(CAPS_BADRANK, self).__init__(msg)
# 
# class EGADS_TESSTATE(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_TESSTATE, self).__init__(msg)
# 
# class EGADS_EXISTS(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_EXISTS, self).__init__(msg)
# 
# class EGADS_ATTRERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_ATTRERR, self).__init__(msg)
# 
# class EGADS_TOPOCNT(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_TOPOCNT, self).__init__(msg)
# 
# class EGADS_BADSCALE(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_BADSCALE, self).__init__(msg)
# 
# class EGADS_NOTORTHO(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTORTHO, self).__init__(msg)
# 
# class EGADS_DEGEN(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_DEGEN, self).__init__(msg)
# 
# class EGADS_CONSTERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_CONSTERR, self).__init__(msg)
# 
# class EGADS_TOPOERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_TOPOERR, self).__init__(msg)
# 
# class EGADS_GEOMERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_GEOMERR, self).__init__(msg)
# 
# class EGADS_NOTBODY(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTBODY, self).__init__(msg)
# 
# class EGADS_WRITERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_WRITERR, self).__init__(msg)
# 
# class EGADS_NOTMODEL(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTMODEL, self).__init__(msg)
# 
# class EGADS_NOLOAD(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOLOAD, self).__init__(msg)
# 
# class EGADS_RANGERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_RANGERR, self).__init__(msg)
# 
# class EGADS_NOTGEOM(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTGEOM, self).__init__(msg)
# 
# class EGADS_NOTTESS(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTTESS, self).__init__(msg)
# 
# class EGADS_EMPTY(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_EMPTY, self).__init__(msg)
# 
# class EGADS_NOTTOPO(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTTOPO, self).__init__(msg)
# 
# class EGADS_REFERCE(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_REFERCE, self).__init__(msg)
# 
# class EGADS_NOTXFORM(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTXFORM, self).__init__(msg)
# 
# class EGADS_NOTCNTX(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTCNTX, self).__init__(msg)
# 
# class EGADS_MIXCNTX(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_MIXCNTX, self).__init__(msg)
# 
# class EGADS_NODATA(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NODATA, self).__init__(msg)
# 
# class EGADS_NONAME(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NONAME, self).__init__(msg)
# 
# class EGADS_INDEXERR(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_INDEXERR, self).__init__(msg)
# 
# class EGADS_MALLOC(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_MALLOC, self).__init__(msg)
# 
# class EGADS_NOTOBJ(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NOTOBJ, self).__init__(msg)
# 
# class EGADS_NULLOBJ(CAPSError):
#     def __init__(self, msg=None):
#         super(EGADS_NULLOBJ, self).__init__(msg)
