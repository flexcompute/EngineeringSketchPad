# This software has been cleared for public release on 05 Nov 2020, case number 88ABW-2020-3462.

# Import Cython *.pxds 
cimport cEGADS

#cimport cAIMUtil

# Python version
# try:
#     from cpython.version cimport PY_MAJOR_VERSION
# except: 
#     print "Error: Unable to import cpython.version"
# #    raise ImportError
# 
# # Import Python modules
# try:
#     import sys
# except:
#     print "\tUnable to import sys module - Try ' pip install sys '"
# #    raise ImportError
# 
# try:
#     import os
# except:
#     print "\tUnable to import os module - Try ' pip install os '"
# #    raise ImportError
# 
# try: 
#     from pyNastran.op2.op2 import OP2
# except: 
#     print "\tUnable to import pyNastran module - Try ' pip install pyNastran '"
#     #    raise ImportError

from cpython.version cimport PY_MAJOR_VERSION

from pyNastran.op2.op2 import OP2

# Handles all byte to uni-code and and byte conversion. 
def _str(data,ignore_dicts = False):
        
    if PY_MAJOR_VERSION >= 3:
        if isinstance(data, bytes):
            return data.decode("utf-8")
   
        if isinstance(data,list):
            return [_str(item,ignore_dicts=True) for item in data]
        if isinstance(data, dict) and not ignore_dicts:
            return {
                    _str(key,ignore_dicts=True): _str(value, ignore_dicts=True)
                    for key, value in data.iteritems()
                    }
    return data

# Get objective results 
cdef public int nastran_getObjective(const char *filename, int *numObjective, double **objective) except -1:
    
    cdef:
        int numData
        double *data=NULL
  
    results = OP2(debug=False)
    results.read_op2(_str(filename))
    objFunc = results.op2_results.responses.convergence_data.obj_initial
    numData = int(len(objFunc))
    data = <double *> cEGADS.EG_alloc(numData*sizeof(double))
    if not data: 
        raise MemoryError()
    
    i = 0
    while i < numData:
        data[i] = objFunc[i]
        i += 1
    
    numObjective[0] = numData
    objective[0] = data
    
    return 0