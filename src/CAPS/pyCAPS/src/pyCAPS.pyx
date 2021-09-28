#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

## \namespace pyCAPS
# Python extension module for CAPS
#
#\mainpage Introduction
# \tableofcontents
# \section overviewpyCAPS Overview
# pyCAPS is a Python extension module to interact with Computational
# Aircraft Prototype Syntheses (CAPS) routines in the
# Python environment. Written in Cython, pyCAPS natively handles all type 
# conversions/casting, while logically grouping CAPS function calls together to simplify 
# a user's experience. Additional functionality not directly available through the 
# CAPS API (such has saving a geometric view) is also provided.  
#
# An overview of the basic pyCAPS functionality is provided in \ref gettingStarted. 
#
# \section differencespyCAPS Key differences between pyCAPS and CAPS
# - Manipulating the "owner" information for CAPS objects isn't currently supported
#
# - CAPS doesn't natively support an array of string values (an array of strings is viewed by CAPS as a single 
# concatenated string), however pyCAPS does. If a list of strings is provided this list is concatenated, separated by a ';' and 
# provided to CAPS as a single string. The number of rows and columns are correctly set to match the original list.
# If a string is received from CAPS by pyCAPS and the rows and columns are set correctly it will be unpacked correctly considering
# entries are separated by a ';'. (Important) If the rows and columns aren't set correctly and the string contains a ';',
#  the data will likely be unpacked incorrectly or raise an indexing error. (Note: not available when setting attributes on objects)
# 
# \section clearancepyCAPS Clearance Statement
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

# Cython declarations
#### cython: language_level=3

# C-Imports 
cimport cCAPS
cimport cEGADS
cimport cOCSM

# Python version
try:
    from cpython.version cimport PY_MAJOR_VERSION
except: 
    print("Error: Unable to import cpython.version")
    raise ImportError

# C-Memory 
try: 
    from libc.stdlib cimport malloc, realloc, free
except:
    print("Error: Unable to import libc.stdlib")
    raise ImportError

# C-string
#try: 
#    from libc.string cimport strcat # memcpy, strncasecmp strcat, strncat, memset, memchr, memcmp, memcpy, memmove
#except:
#    print("Error: Unable to import libc.string")
#    raise ImportError

# C-stdio
try:
    from libc.stdio cimport sprintf
except:
    print("Error: Unable to import libc.stdio")
    raise ImportError

# JSON string conversion       
try: 
    from json import dumps as jsonDumps
    from json import loads as jsonLoads
except:
    print("Error:  Unable to import json\n")
    raise ImportError

# OS module 
try:
    import os
except:
    print("Error: Unable to import os\n")
    raise ImportError

# Time module
try: 
    import time
except: 
    print("Error:  Unable to import time\n")
    raise ImportError

# Math module
try:
    from math import sqrt
except:
    print("Error: Unable to import math\n")
    raise ImportError

# atexit module
try:
    import atexit
except:
    print("Error: Unable to import atexit\n")
    raise ImportError

# Reg. expression module 
try: 
    from re import search as regSearch
    #from re import findall as regFindAll
    from re import M as regMultiLine
    from re import I as regIgnoreCase
    from re import compile as regCompile
     
except:
    print("Error: Unable to import re\n")
    raise ImportError

# Set version 
__version__ = "2.2.1"

# ----- Include files ----- #

# All functions involved with forming HTML tree are include file  
include "treeUtils.pxi"

# Functions to handle unit conversions
include "unitsUtils.pxi"

# Functions to deal with attributes 
include "attributeUtils.pxi"

# Class for exceptions 
include "exceptionUtils.pxi"

# Functions for web browser data viewer
include "webViewerUtils.pxi"

# All functions involved with colorMap selection
include "colorMapUtils.pxi"

# Functions to help when working with Tecplot files
include "tecplotUtils.pxi"

# ------------------------- #

# Print out any errors
cdef void report_Errors(int nErr, cCAPS.capsErrs * errors):
    
    cdef:
        int i = 0, j = 0 # Indexing 
        int status = 0 # Function return status
    
        int nLines 
        char **lines
        cCAPS.capsObj obj
    
    if errors != NULL:
        for i in range(nErr):

            status = cCAPS.caps_errorInfo(errors, i+1, &obj, &nLines, &lines)
            if status != 0: 
                print(" CAPS Error: caps_errorInfo[%d] = %d!\n"%(i+1, status))
                continue
            
            for j in range(nLines):
                 print(" %s\n"%(lines[j]))
        
        cCAPS.caps_freeError(errors)

# Handles all uni-code and and byte conversion. Originally intended to just convert 
#     JSON unicode returns to str - taken from StackOverFlow question - "How to get string 
#     objects instead of Unicode ones from JSON in Python?"
cpdef _byteify(data,ignore_dicts = False):
        
        if isinstance(data, bytes):
            return data
        
        if PY_MAJOR_VERSION < 3:
            if isinstance(data,unicode):
                return data.encode('utf-8')
        else:
            if isinstance(data,str):
                return data.encode('utf-8')
        
        if isinstance(data,list):
            return [_byteify(item,ignore_dicts=True) for item in data]
        if isinstance(data, dict) and not ignore_dicts:
            return {
                    _byteify(key,ignore_dicts=True): _byteify(value, ignore_dicts=True)
                    for key, value in data.iteritems()
                    }
        return data

# Handles all byte to uni-code conversion. 
cpdef _strify(data,ignore_dicts = False):
        
    if PY_MAJOR_VERSION >= 3:
        if isinstance(data, bytes):
            return data.decode("utf-8")
    
        if isinstance(data,list):
            return [_strify(item,ignore_dicts=True) for item in data]
        if isinstance(data, dict) and not ignore_dicts:
            return {
                    _strify(key,ignore_dicts=True): _strify(value, ignore_dicts=True)
                    for key, value in data.iteritems()
                    }
    return data

# Convert a returned CAPS value (void *) to an appropriate Python object
cdef object castValue2PythonObj(const void *data, 
                                cCAPS.capsvType valueType, 
                                int numRow, 
                                int numCol):
    
    cdef:
        int length = 0
        object valueOut, valueTemp
        cCAPS.capsTuple tupleTemp
        char *byteString
        object rowMajorSwitched = False
    
    #print("Data", <double> (<double *> data)[0])
    #print("Nrow = ", numRow, "Ncol = ", numCol)
   
    # what else to do with thing...
    if data == NULL: return None
 
    if numRow == 1 and numCol == 1:
        valueOut = 0
    else:
        
        if numRow == 1 or numCol == 1:
            if numRow != 1:
                valueOut = []*numRow
            if numCol != 1:
                valueOut = []*numCol
                numRow = numCol
                numCol = 1
        else:
            valueOut = [[] for i in range(numRow)]
            
    for i in range(numRow):
        for j in range(numCol):
            
            # Boolean
            if valueType == cCAPS.Boolean:
                
                if <int> ((<int *> data)[length]) == <int> cCAPS.true:
                    valueTemp = True
                else: 
                    valueTemp = False
                    
                length += 1
            
            # Integer 
            elif valueType == cCAPS.Integer:
                
                valueTemp = <object> ((<int *> data)[length])
                length += 1
            
            # Double          
            elif valueType == cCAPS.Double:
                
                valueTemp = <object> ((<double *> data)[length])
                length += 1
            
            # String          
            elif valueType == cCAPS.String:
                byteString = <char *> data
                
                if byteString:
                    valueTemp = _strify(_strify(byteString).split(';')[length])
                else:
                    valueTemp = None
                    
                length += 1
            
            # Tuple
            elif valueType == cCAPS.Tuple:
                #print("I'm a tuple")
                
                tupleTemp = (<cCAPS.capsTuple *> data)[length]
                length += 1
                
                if tupleTemp.name:
                    keyWord = _strify( <object>  tupleTemp.name)
                else:
                    keyword = None
                
                if tupleTemp.value:
                    try:
                        keyValue = _strify(jsonLoads( _strify(<object> tupleTemp.value), object_hook=_byteify))
                    except:
                        keyValue = <object> tupleTemp.value
                        
                        if isinstance(keyValue, bytes):
                            keyValue = _strify(keyValue)
                else:
                    keyValue = None
                            
                valueTemp = (keyWord, keyValue) 
                
            # Unknown type - Value type most likely 
            else:
                print("Can not convert type (", valueType, ") to Python Object!")
                return None
            
            if isinstance(valueOut, list):
                if numCol != 1: #if isinstance(valueOut[0], list):
                    
                    valueOut[i].append(valueTemp)
                else:
                    valueOut.append(valueTemp)
            else: 
                valueOut = valueTemp
    
    return valueOut

# Convert a Python object to a CAPS value that is then cast to a void * - returned value (if not a String)
# must be freed to avoid memory leaks
cdef void * castValue2VoidP(object data,
                            object byteData, # Make sure the byteified data stays in memory long enough to set it
                            cCAPS.capsvType valueType):
    cdef: 
        int length = 0
        
        int numRow
        int numCol
    
        int *valInt = NULL
        double *valDouble = NULL
        char *varChar = NULL
        cCAPS.capsTuple *valTuple = NULL
        object tempJSON = [] # List of temporary JSON strings 
        
    
    # Convert data to list 
    if not isinstance(data, list):
        data = [data] # Convert to list
    
    numRow, numCol = getPythonObjShape(data)
        
    length = numRow*numCol

    #print("castValue2VoidP Nrow = ", numRow, "Ncol = ", numCol)

    # Boolean
    if valueType == cCAPS.Boolean:
        #print("I'm a bool")
        
        valInt = <int *> malloc(length*sizeof(int))
        if not valInt:
            raise MemoryError()
    
        length = 0
        for i in range(numRow):
            for j in range(numCol):
            
                if numCol == 1:
                    valInt[length] = data[i]
                else:
                    valInt[length] = data[i][j]
                
                length += 1
                
        return <void *> valInt
    
    # Integer
    elif valueType == cCAPS.Integer:
        #print("I'm a integer")
        
        valInt = <int *> malloc(length*sizeof(int))
        if not valInt:
            raise MemoryError()
        
        length = 0
        for i in range(numRow):
            for j in range(numCol):
                
                if numCol == 1:
                    valInt[length] = data[i]
                else:
                    valInt[length] = data[i][j]
                
                length += 1
        
        return <void *> valInt
       
    # Double
    elif valueType == cCAPS.Double:
        #print("I'm a double")
        
        valDouble = <double *> malloc(length*sizeof(double))
        
        if not valDouble:
            raise MemoryError()
        
        length = 0
        for i in range(numRow):
            for j in range(numCol):
                
                if numCol == 1:
                    valDouble[length] = data[i]
                else:
                    valDouble[length] = data[i][j]
                  
                length += 1

        return <void *> valDouble
        
    # String
    elif valueType == cCAPS.String:
        #print("I'm a string")
        byteData.append(_byteify(';'.join( [ i if (isinstance(i,str)) else jsonDumps(i) for i in data] )))
        valChar = <char *> byteData[-1]
        
        return <void *> valChar
   
    # Tuple
    elif valueType == cCAPS.Tuple:
        #print("I'm a tuple")

        valTuple = <cCAPS.capsTuple *> malloc(length*sizeof(cCAPS.capsTuple))
        
        if not valTuple:
            raise MemoryError()
        
        for i in range(numRow):
           
            if (isinstance(data[i][1], bytes) or 
                isinstance(data[i][1], str)  or 
                isinstance(data[i][1], int)  or 
                isinstance(data[i][1], float)):
                
                # Need data to be a string for cast to char * below - otherwise we get a TypeError
                if isinstance(data[i][1], int)  or isinstance(data[i][1], float):
                    temp = str(data[i][1])
                else:
                    temp = data[i][1] 
                
                byteString =  _byteify(temp)
            else:
        
                tempJSON.append(jsonDumps(data[i][1]))
                byteString =  _byteify(tempJSON[-1])

            byteData.append((_byteify(data[i][0]), byteString))
                           
            valTuple[<int> i].name = <char *> byteData[-1][0]
            valTuple[<int> i].value = <char *> byteData[-1][1]
            
        return <void *> valTuple
        
    # Unknown type - Value type most likely 
    else:
        print("Can not convert Python object to type", valueType)
        raise TypeError

# Determine the equivalent CAPS value object type for a Python object     
cdef cCAPS.capsvType sortPythonObjType(object obj):

    # If we have a list loop back in to see what type of value the element is 
    if isinstance(obj, list):
        return sortPythonObjType(obj[0]) 
   
    if isinstance(obj, bool) and (str(obj) == 'True' or str(obj) == 'False'):
        return cCAPS.Boolean
    
    if isinstance(obj, int):
        return cCAPS.Integer
    
    if isinstance(obj, float):
        return cCAPS.Double
    
    if isinstance(obj, tuple):
        return cCAPS.Tuple
        
    return cCAPS.String

cdef object getPythonObjShape(object obj):
    
    dataValue = obj
    
    numRow = 0
    numCol = 0
    
    if not isinstance(dataValue, list):
        dataValue = [dataValue] # Convert to list
            
    numRow = len(dataValue)
        
    if isinstance(dataValue[0], list):
        numCol = len(dataValue[0]) 
    else:
        numCol = 1
        
    # Make sure shape is consistent
    if numCol != 1:
        for i in range(numRow):
            if len(dataValue[i]) != numCol:
                print("Inconsistent list sizes!", i, dataValue[i], numCol, dataValue[0])
                raise CAPSError(cCAPS.CAPS_BADVALUE, msg="Inconsistent list sizes!")
    
    return numRow, numCol        

# Save the an array of bodies to a file
cdef int saveBodies(int numBody, cEGADS.ego *bodies, object filename):

    cdef: 
        int status = 0 # Function return status
        int i = 0 # Indexing

        cEGADS.ego *bodiesCopy = NULL # Model bodies
        cEGADS.ego model
        cEGADS.ego context
        
    if numBody < 1:
        print("The number of bodies in is less than 1!!!\n")
        return cCAPS.CAPS_BADVALUE

    cEGADS.EG_getContext(bodies[0], &context)

    # Allocate body copy array
    bodiesCopy = <cEGADS.ego *> malloc(numBody*sizeof(cEGADS.ego))
    if not bodiesCopy:
        raise MemoryError()

    # Make copies of the bodies to hand to the model
    for i in range(numBody):
        status = cEGADS.EG_copyObject(bodies[i], NULL, &bodiesCopy[i])
        if status != cEGADS.EGADS_SUCCESS:
            free(bodiesCopy)
            return status

    # Make a model from the context and bodies
    status = cEGADS.EG_makeTopology(context, 
                                    <cEGADS.ego> NULL, 
                                    cEGADS.MODEL, 
                                    0, 
                                    <double *> NULL, 
                                    numBody, 
                                    bodiesCopy, 
                                    <int *> NULL,
                                    &model)
    if status != cEGADS.EGADS_SUCCESS:
        free(bodiesCopy)
        return status
        
    # Save the model
    filename = _byteify(filename)
    status = cEGADS.EG_saveModel(model, <char *>filename)
    if status != cEGADS.EGADS_SUCCESS:
        free(bodiesCopy)
        return status

    # Clean up - delete the model
    status = cEGADS.EG_deleteObject(model)
    if status != cEGADS.EGADS_SUCCESS:
        free(bodiesCopy)
        return status

    if bodiesCopy:
        free(bodiesCopy)

    return cCAPS.CAPS_SUCCESS

# Save the geometry of the current problem 
cdef int saveGeometry(cCAPS.capsObj pObj, object filename):

    cdef: 
        int status = 0 # Function return status
        int i = 0 # Indexing
        int numBody=0 
        
        cCAPS.capsProblem *problem = NULL
        
        cCAPS.capsErrs *errors = NULL
        int nErr
        
        cEGADS.ego body # Temporary body container
        cEGADS.ego *bodies = NULL # Model bodies
        cEGADS.ego model
        
        char *units = NULL

    # Check the problem object
    if pObj.type != cCAPS.PROBLEM:
        return cCAPS.CAPS_BADTYPE
     
    if pObj.blind == NULL:
        return cCAPS.CAPS_NULLBLIND
    
    # Get pointer to capsProblem from blind 
    problem = <cCAPS.capsProblem *> pObj.blind

    # Force regeneration of the geometry
    status = cCAPS.caps_preAnalysis(pObj, &nErr, &errors)

    # Check errors
    report_Errors(nErr, errors)

    if (status != cCAPS.CAPS_SUCCESS and 
        status != cCAPS.CAPS_CLEAN):
        return status
    
    # How many bodies do we have
    status = cCAPS.caps_size(pObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
    if status != cCAPS.CAPS_SUCCESS:
        return status

    if numBody < 1:
                
        print("The number of bodies in the problem is less than 1!!!\n")
        return cCAPS.CAPS_BADVALUE
    
    # Allocate body array
    bodies = <cEGADS.ego *> malloc(numBody*sizeof(cEGADS.ego))
    if not bodies:
        raise MemoryError()

    # Make copies of the bodies to hand to the model
    for i in range(numBody):
        status = cCAPS.caps_bodyByIndex(pObj, i+1, &bodies[i], &units)
        if status != cCAPS.CAPS_SUCCESS:
            free(bodies)
            return status
    
    status = saveBodies(numBody, bodies, filename);
    
    if bodies: free(bodies)
    
    return status

def _storeIteration(initial, analysisDir, noFileDir=True):
    
    def initialInstance(analysisDir):
        
        # Get current directory
        currentDir = os.getcwd()
        
        # Change to analysis directory
        os.chdir(analysisDir)
        
        # Determine instance number
        instance = 0
        key = "Instance_"
        
        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]
        
        for i in currentFiles:
            
            path = os.path.join(os.getcwd(), i)
    
            if os.path.isdir(path) and key in i:
                keepFiles.remove(i)
                
                if int(i.replace(key, "")) >= instance:
                    instance = int(i.replace(key, "")) + 1
            
        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 
    
        # Make new instance directory
        instance = os.path.join(os.getcwd(), key + str(instance))
        os.mkdir(instance)
        
        for i in keepFiles:
            # print("Change file - ", os.path.join(instance,i))
            os.rename(i, os.path.join(instance,i))
                
        # Change back to current directory
        os.chdir(currentDir)
    
    if initial:
        initialInstance(analysisDir)
        return 
    
    # Get current directory
    currentDir = os.getcwd()
    
    # Change to analysis directory
    os.chdir(analysisDir)
    
    # Determine instance number
    interation = 0
    key = "Iteration_"
    key2 = "Instance_"
    
    currentFiles = os.listdir(os.getcwd());
    keepFiles = currentFiles[:]
    
    for i in currentFiles:
        path = os.path.join(os.getcwd(), i)

        if os.path.isdir(path):
            
            if key in i or key2 in i:
                keepFiles.remove(i)
            
            if key in i:
                if int(i.replace(key, "")) >= interation:
                    interation = int(i.replace(key, "")) + 1
                    
    if noFileDir == False:
        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 
    
    # Make new interation directory
    interation = os.path.join(os.getcwd(), key + str(interation))
    if not os.path.isdir(interation):
        os.mkdir(interation)
    
    # Move files into new interation folder
    for i in keepFiles:
        # print"Change file - ", os.path.join(interation,i)
        os.rename(i, os.path.join(interation,i))
    
    # Change back to current directory
    os.chdir(currentDir)
        
cdef object createOpenMDAOComponent(analysisObj, 
                                    inputParam, 
                                    outputParam, 
                                    executeCommand = None, 
                                    inputFile = None, 
                                    outputFile = None, 
                                    sensitivity = {}, 
                                    changeDir = None, 
                                    saveIteration = False):
    
    try: 
        from openmdao import __version__ as version
    except:
        print("Error: Unable to determine OpenMDAO version!")
        return cCAPS.CAPS_NOTFOUND 
    try:
        from numpy import ndarray
    except:
        print("Error: Unable to import numpy!")
    
    print("OpenMDAO version - ", version)
    versionMajor = int(version.split(".")[0])
    versionMinor = int(version.split(".")[1])
    
    if versionMajor >=2:
        try: 
            from  openmdao.api import ExternalCodeComp, ExplicitComponent     
        except: 
            print("Error: Unable to import ExternalCodeComp or ExplicitComponent from openmdao.api!")
            return cCAPS.CAPS_NOTFOUND 
    
    elif versionMajor == 1 and versionMinor == 7:   
        try: 
            from  openmdao.api import ExternalCode, Component     
        except: 
            print("Error: Unable to import ExternalCode or Component from openmdao.api!")
            return cCAPS.CAPS_NOTFOUND 
    else:
        print("Error: Unsupported OpenMDAO version!")
        return cCAPS.CAPS_BADVALUE
    
    # Change to lists if not lists already 
    if not isinstance(inputParam, list):
        inputParam = [inputParam]
        
    if not isinstance(outputParam, list):
        outputParam = [outputParam]

    def filterInputs(varName, validInput, validGeometry):
        
        value = None
        
        if ":" in varName: # Is the input trying to modify a Tuple?
          
            analysisTuple = varName.split(":")[0] 
            if analysisTuple not in validInput.keys():
                if varName not in validGeometry.keys(): #Check geometry
                    print("It appears the input parameter", varName, "is trying to modify a Tuple input, but", analysisTuple,
                          "this is not a valid ANALYSISIN or GEOMETRYIN variable! It will not be added to the OpenMDAO component")
                    return None
                else: 
                    value = validGeometry[varName]
            else:  
                
                if ":" in varName:
                    tuples = validInput[analysisTuple]
                    if not isinstance(tuples,list):
                        tuples = [tuples]
                    
                    tupleValue = varName.split(":")[1] 
                    
                    for i in tuples:
                        if tupleValue == i[0]:
                            if varName.count(":") == 2:
                                dictValue = varName.split(":")[2] 
                                value = i[1][dictValue]
                            elif varName.count(":") == 1:
                                value = i[1]
        else:
            # Check to make sure variable is in analysis
            if varName not in validInput.keys():
                if varName not in validGeometry.keys():
                    print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                          "It will not be added to the OpenMDAO component.")
                    return None
                else:
                    value = validGeometry[varName]
            else:
                value = validInput[varName]
                
        # Check for None
        if value is None:
            value = 0.0
        
        return value
    
    def filterOutputs(varName, validOutput):
        
        value = None
        
        # Check to make sure variable is in analysis
        #if varName not in validOutput.keys():
        if varName not in validOutput:
            print("Output parameter", varName, "is not a valid ANALYSISOUT variable!",
                  "It will not be added to the OpenMDAO component.")
            return None
        #else:
        #    value = validOutput[varName]
        
        # Check for None
        if value is None:
            value = 0.0
        
        return value
    
    def setInputs(varName, varValue, analysisObj, validInput, validGeometry):
       
        if isinstance(varValue, ndarray):
            varValue = varValue.tolist()
        
        if isinstance(varValue,list):
            numRow, numCol = getPythonObjShape(varValue)
            if numRow == 1 and numCol == 1:
                    varValue = varValue[0]
            elif numRow == 1 and numCol > 1:
                varValue = varValue[0]
                
        if ":" in varName: # Is the input trying to modify a Analysis Tuple?
            
            analysisTuple = varName.split(":")[0] 
            
            if analysisTuple not in validInput.keys():
                if varName in validGeometry.keys(): #Check geometry
                    analysisObj.capsProblem.geometry.setGeometryVal(varName, varValue)
                else:
                    print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                          "It will not be modified to the OpenMDAO component. This should have already been reported!")
            else:
                if varName.count(":") == 2:
                   
                    value = setElementofCAPSTuple(analysisObj.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1], # Tuple key
                                                  varValue, # Value to set 
                                                  varName.split(":")[2]) # Dictionary key
                else: 
                    
                    value = setElementofCAPSTuple(analysisObj.getAnalysisVal(analysisTuple), # Validity has already been checked
                                                  varName.split(":")[1],  # Tuple key
                                                  varValue) # Value to set 
                
                analysisObj.setAnalysisVal(analysisTuple, value)
            
        else:
            
            if varName in validInput.keys():
                
                analysisObj.setAnalysisVal(varName, varValue)
            
            elif varName in validGeometry.keys():
                
                analysisObj.capsProblem.geometry.setGeometryVal(varName, varValue)
            else:
            
                print("Input parameter", varName, "is neither a valid ANALYSISIN nor GEOMETRYIN variable!", 
                      "It will not be modified to the OpenMDAO component. This should have already been reported!")    
    
#     def gatherParents(analysis, analysisParents):
#         
#         analysisParents += analysis.parents
#         
#         if analysis.parents:
#             for parent in analysis.parents:
#                 analysisParents = gatherParents(analysis.capsProblem.analysis[parent], analysisParents)
#            
#         return set(analysisParents)
#     
    def checkParentExecution(analysis, checkInfo = False):
        
        def infoDictCheck(analysis):
            
            infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)
        
            if infoDict["executionFlag"] == True:
                return True
            else:
                print(analysis.aimName, "can't self execute!")
                return False
        
        if analysis.parents: # If we have parents 
            
            # Check analysis parents
            for parent in analysis.parents:
        
                if not checkParentExecution(analysis.capsProblem.analysis[parent], True):
                    return False
                
                if checkInfo:
                    return infoDictCheck(analysis)
            
        if checkInfo:
            return infoDictCheck(analysis)
    
    def executeParent(analysis, runAnalysis = False):
        
        def executeAnalysis(analysis):
            infoDict = analysis.getAnalysisInfo(printInfo = False, infoDict = True)
                
            # Only execute analsysis if it isn't up to date
            if infoDict["status"] != 0: # 0 == Up to date
                analysis.preAnalysis()
                analysis.postAnalysis()
            
            return
            
        if analysis.parents:
            for parent in analysis.parents:
                executeParent(analysis.capsProblem.analysis[parent], runAnalysis = True)
        
            if runAnalysis:
                executeAnalysis(analysis)
            
            return 
        
        if runAnalysis:
            executeAnalysis(analysis)
        
        return 
    
    def initialInstance(analysisDir):
        
         # Get current directory
        currentDir = os.getcwd()
        
        # Change to analysis directory
        os.chdir(analysisDir)
        
        # Determine instance number
        instance = 0
        key = "Instance_"
        
        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]
        
        for i in currentFiles:
            
            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path) and key in i:
                keepFiles.remove(i)
                
                if int(i.replace(key, "")) >= instance:
                    instance = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 

        # Make new instance directory
        instance = os.path.join(os.getcwd(), key + str(instance))
        os.mkdir(instance)
        
        for i in keepFiles:
            # print"Change file - ", os.path.join(instance,i)
            os.rename(i, os.path.join(instance,i))
                
        # Change back to current directory
        os.chdir(currentDir)
        
    def storeIteration(analysisDir):
        
        # Get current directory
        currentDir = os.getcwd()
        
        # Change to analysis directory
        os.chdir(analysisDir)
        
        # Determine instance number
        interation = 0
        key = "Iteration_"
        key2 = "Instance_"
        
        currentFiles = os.listdir(os.getcwd());
        keepFiles = currentFiles[:]
        
        for i in currentFiles:
            path = os.path.join(os.getcwd(), i)

            if os.path.isdir(path):
                
                if key in i or key2 in i:
                    keepFiles.remove(i)
                
                if key in i:
                    if int(i.replace(key, "")) >= interation:
                        interation = int(i.replace(key, "")) + 1

        # If we don't have any files nothing needs to be done
        if not keepFiles:
            # Change back to current directory
            os.chdir(currentDir)
            return 
        
        # Make new interation directory
        interation = os.path.join(os.getcwd(), key + str(interation))
        if not os.path.isdir(interation):
            os.mkdir(interation)
        
        # Move files into new interation folder
        for i in keepFiles:
            # print"Change file - ", os.path.join(interation,i)
            os.rename(i, os.path.join(interation,i))
        
        # Change back to current directory
        os.chdir(currentDir)
    
    if versionMajor >=2:
        class openMDAOExternal_V2(ExternalCodeComp):      
            def __init__(self, 
                         analysisObj, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal_V2, self).__init__()
                
                self.analysisObj = analysisObj
                self.inputParam = inputParam
                self.outputParam = outputParam
                self.executeCommand = executeCommand
                self.inputFile = inputFile
                self.outputFile = outputFile
                self.sensitivity = sensitivity
                self.changeDir = changeDir
                self.saveIteration = saveIteration
                
                # Get check to make sure         
                self.validInput    = self.analysisObj.getAnalysisVal()
                self.validOutput   = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()
            
            def setup(self):
                
                for i in self.inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
             
                    self.add_input(i, val=value)
                
                for i in self.outputParam:
                    value = filterOutputs(i, self.validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if self.inputFile is not None:
                    if isinstance(self.inputFile, list):
                        self.options['external_input_files'] = self.inputFile
                    else:
                        self.options['external_input_files'] = [self.inputFile]
                
                if self.outputFile is not None:
                    if isinstance(self.outputFile, list):
                        self.options['external_output_files'] = self.outputFile
                    else:
                        self.options['external_output_files'] = [self.outputFile]
               
                if isinstance(self.executeCommand, list):
                    self.options['command'] = self.executeCommand
                else:
                    self.options['command'] = [self.executeCommand]
                    
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in self.sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = self.sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = self.sensitivity[i]
                    elif "form" in i: 
                        form = self.sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = self.sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                  
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
#                  self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)
                
                #self.changeDir = changeDir
                
                #self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)           
                     
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)
                
                # Run parents 
                executeParent(self.analysisObj, runAnalysis = False)
                
                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)
                        
                self.analysisObj.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)
                    
                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysisObj.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)
            
        class openMDAOExternal_V2ORIGINAL(ExternalCodeComp):      
            def __init__(self, 
                         analysisObj, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal_V2, self).__init__()
                
                self.analysisObj = analysisObj
                
                # Get check to make sure         
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()
                
                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
             
                    self.add_input(i, val=value)
                
                for i in outputParam:
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]
                
                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]
               
                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]
                    
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i: 
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                  
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
#                 self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                self.declare_partials('*', '*', method=method, step=step, form=form, step_calc=step_calc)
                
                self.changeDir = changeDir
                
                self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)           
                     
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)
                
                # Run parents 
                executeParent(self.analysisObj, runAnalysis = False)
                
                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)
                        
                self.analysisObj.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)
                    
                #Parent compute function actually runs the external code
                super(openMDAOExternal_V2, self).compute(params, unknowns)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysisObj.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
    #         def compute_partials(self, inputs, partials):
    #             outputs = {}
    #      the parent compute function actually runs the external code
    #         super(ParaboloidExternalCodeCompDerivs, self).compute(inputs, outputs)
    # 
    #         # parse the derivs file from the external code and set partials
    #         with open(self.derivs_file, 'r') as derivs_file:
    #             partials['f_xy', 'x'] = float(derivs_file.readline())
    #             partials['f_xy', 'y'] = float(derivs_file.readline())
                
        class openMDAOComponent_V2(ExplicitComponent):
                
            def __init__(self, 
                         analysisObj, 
                         inputParam, outputParam, 
                         sensitivity = {}, 
                         saveIteration = False):
                    
                super(openMDAOComponent_V2, self).__init__()
    
                self.analysisObj = analysisObj
                    
                # Get check to make sure         
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()
                
                for i in inputParam:
                    
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_input(i, val=value)
                
                for i in outputParam:
                    
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                   
                method = "fd"
                step= None
                form = None
                step_calc = None
                
                for i in sensitivity.keys():
                    if "method" in i or "type" in i:
                        method = sensitivity[i]
                    elif "step" == i or "step_size" in i:
                        step = sensitivity[i]
                    elif "form" in i: 
                        form = sensitivity[i]
                    elif "step_calc" in i:
                        step_calc = sensitivity[i]
                    else:
                        print("Warning: Unreconized setSensitivity key - ", i)
                
                if method != "fd":
                    print("Only finite difference is currently supported!")
                    return cCAPS.CAPS_BADVALUE
                
                self.declare_partials(of='*', wrt='*', method=method, step=step, form=form, step_calc=step_calc)
                
                self.saveIteration = saveIteration
                
                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)  
               
            def compute(self, params, unknowns):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)
               
                # Run parents 
                executeParent(self.analysisObj, runAnalysis = False)
                
                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)
                    
                self.analysisObj.preAnalysis()
                self.analysisObj.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
    #         def compute_partials(self, inputs, partials):
    #             partials['area', 'length'] = inputs['width']
    #             partials['area', 'width'] = inputs['length']
    
    else:
        class openMDAOExternal(ExternalCode):      
            def __init__(self, 
                         analysisObj, 
                         inputParam, outputParam, executeCommand, 
                         inputFile = None, outputFile = None, 
                         sensitivity = {}, 
                         changeDir = True,
                         saveIteration = False):
                
                super(openMDAOExternal, self).__init__()
                
                self.analysisObj = analysisObj
                
                # Get check to make sure         
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()
                
                for i in inputParam:
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_param(i, val=value)
                
                for i in outputParam:
                    value = filterOutputs(i, validOutput)
                    
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                    
                if inputFile is not None:
                    if isinstance(inputFile, list):
                        self.options['external_input_files'] = inputFile
                    else:
                        self.options['external_input_files'] = [inputFile]
                
                if outputFile is not None:
                    if isinstance(outputFile, list):
                        self.options['external_output_files'] = outputFile
                    else:
                        self.options['external_output_files'] = [outputFile]
               
                if isinstance(executeCommand, list):
                    self.options['command'] = executeCommand
                else:
                    self.options['command'] = [executeCommand]
                    
                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]
                    
                self.changeDir = changeDir
                
                self.saveIteration = saveIteration     
                
                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)           
                     
            def solve_nonlinear(self, params, unknowns, resids):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)
               
                # Run parents 
                executeParent(self.analysisObj, runAnalysis = False)
                
                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)
                        
                self.analysisObj.preAnalysis()
                    
                currentDir = None
                # Change the directory for the analysis - unless instructed otherwise 
                if self.changeDir:
                    currentDir = os.getcwd()
                    os.chdir(self.analysisObj.analysisDir)
                    
                #Parent solve_nonlinear function actually runs the external code
                super(openMDAOExternal, self).solve_nonlinear(params, unknowns, resids)
                
                # Change the directory back after we are done with the analysis - if it was changed
                if currentDir is not None: 
                    os.chdir(currentDir)
                    
                self.analysisObj.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J
        
        class openMDAOComponent(Component):
                
            def __init__(self, 
                         analysisObj, 
                         inputParam, outputParam, 
                         sensitivity = {}, 
                         saveIteration = False):
                    
                super(openMDAOComponent, self).__init__()
    
                self.analysisObj = analysisObj
                    
                # Get check to make sure         
                self.validInput  = self.analysisObj.getAnalysisVal()
                validOutput      = self.analysisObj.getAnalysisOutVal(namesOnly = True)
                self.validGeometry = self.analysisObj.capsProblem.geometry.getGeometryVal()
                
                for i in inputParam:
                    
                    value = filterInputs(i, self.validInput, self.validGeometry)
                    if value is None:
                        continue 
                    
                    self.add_param(i, val=value)
                
                for i in outputParam:
                    
                    value = filterOutputs(i, validOutput)
                    if value is None:
                        continue 
                    
                    self.add_output(i, val=value)
                   
                for i in sensitivity.keys():
                    self.deriv_options[i] = sensitivity[i]
                    
                self.saveIteration = saveIteration
                
                if self.saveIteration:
                    initialInstance(self.analysisObj.analysisDir)  
        
            def solve_nonlinear(self, params, unknowns, resids):
                
                # Set inputs
                for i in params.keys():
                    setInputs(i, params[i], self.analysisObj, self.validInput, self.validGeometry)
               
                # Run parents 
                executeParent(self.analysisObj, runAnalysis = False)
                
                if self.saveIteration:
                    storeIteration(self.analysisObj.analysisDir)
                    
                self.analysisObj.preAnalysis()
                self.analysisObj.postAnalysis()
                
                for i in unknowns.keys():
                    unknowns[i] = self.analysisObj.getAnalysisOutVal(varname = i)
            
            # TO add sensitivities       
            #def linearize(self, params, unknowns, resids):
            #    J = {}
            #    return J
    
    if checkParentExecution(analysisObj) == False: #Parents can't all self execute
        return cCAPS.CAPS_BADVALUE 

    infoDict = analysisObj.getAnalysisInfo(printInfo = False, infoDict = True)
    
    if executeCommand is not None and infoDict["executionFlag"] == True:

        if versionMajor >=2:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCodeComp to ExplicitComponent")
        else:
            print("An execution command was provided, but the AIM says it can run itself! Switching ExternalCode to Component")
    
    if executeCommand is not None and infoDict["executionFlag"] == False:
        
        if versionMajor >=2:
            return openMDAOExternal_V2(analysisObj, inputParam, outputParam, 
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir, 
                                    saveIteration)
    
        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOExternal(analysisObj, inputParam, outputParam, 
                                    executeCommand, inputFile, outputFile, sensitivity, changeDir, 
                                    saveIteration)
    else:
        
        if versionMajor >=2:
           return openMDAOComponent_V2(analysisObj, inputParam, outputParam, sensitivity, saveIteration)
    
        elif versionMajor == 1 and versionMinor == 7:
            return openMDAOComponent(analysisObj, inputParam, outputParam, sensitivity, saveIteration)


## Defines a CAPS problem object.
# A capsProblem is the top-level object for a single mission/problem. It maintains a single set
# of interrelated geometric models (see \ref pyCAPS.capsGeometry), 
# analyses to be executed (see \ref pyCAPS.capsAnalysis), 
# connectivity and data (see \ref pyCAPS.capsBound)
# associated with the run(s), which can be both multi-fidelity and multi-disciplinary.
cdef class capsProblem:
    
    cdef:
        # CAPS problem object 
        cCAPS.capsObj problemObj
        
        # Current function return status.
        readonly int status
        
        # CAPS geometry class for the problem
        readonly capsGeometry geometry
        
        # CAPS analyses classes for the problem
        readonly object analysis
        #readonly object aimGlobalCount
        #readonly object analysisDir
        
        # Data Bound
        readonly object dataBound
        
        # Intent (analysis+geometry) in which to load the AIM
        public object capsIntent
        
        # Raise an exception after an error
        public object raiseException
    
        # CAPS value object in the problem
        readonly object value
        
#     ## Analysis intent dictionary.
#     capsIntent = {
#                         'ALL'           : cCAPS.ALL,
#                         'WAKE'          : cCAPS.WAKE,
#                         'STRUCTURE'     : cCAPS.STRUCTURE,
#                         'LINEARAERO'    : cCAPS.LINEARAERO,
#                         'FULLPOTENTIAL' : cCAPS.FULLPOTENTIAL,
#                         'CFD'           : cCAPS.CFD
#                         }
    
#     ## Geometry representation dictionary.
#     geometricRep = {
#                         'ALL'      : cCAPS.ALL,
#                         'NODE'     : cEGADS.NODE,
#                         'WIREBODY' : cEGADS.WIREBODY,
#                         'FACEBODY' : cEGADS.FACEBODY,
#                         'SHEETBODY': cEGADS.SHEETBODY,
#                         'SOLIDBODY': cEGADS.SOLIDBODY
#                         }
    
    ## Initialize the problem.
    # See \ref problem.py for a representative use case.
    # \param libDir Deprecated option, no longer required. 
    # \param raiseException Raise an exception after a CAPS error is encountered (default - True). See \ref raiseException.
    def __cinit__(self, libDir=None, raiseException=True):
        ## \example problem.py 
        # Basic example use case of the initiation (pyCAPS.capsProblem.__init__) and 
        # termination (pyCAPS.capsProblem.closeCAPS)
        
        if libDir is not None:
            print("\n\nWarning: libDir (during __init__ of capsProblem) is a deprecated option as it\n" +
                   "is no longer needed, please remove this argument from your script.\n" +
                   "Further use will result in an error in future releases!!!!!!\n")
        
        ## Current CAPS status code 
        self.status =  cCAPS.CAPS_SUCCESS
        
        ## Raise an exception after a CAPS error is found (default - True). Disabling (i.e. setting to False) may have unexpected 
        # consequences; in general the value should be set to True. 
        self.raiseException = raiseException
         
        # Geometry 
        ## Geometry object loaded into the problem. Set via \ref loadCAPS.
        self.geometry = None

        # Analysis
        # Number of AIMs loaded into the problem.
        #self.aimGlobalCount = 0

        # Current analysis directory which was used to load the latest AIM.
        #self.analysisDir = None
        
        # Default intent used to load the AIM if capsIntent is not provided.
        #self.capsIntent = None
        
        ## Dictionary of analysis objects loaded into the problem. Set via \ref loadAIM. 
        self.analysis = {}
        
        ## Dictionary of data transfer/bound objects loaded into the problem. Set 
        # via \ref createDataBound or \ref createDataTransfer.
        self.dataBound = {}
        
        ## Dictionary of value objects loaded into the problem. Set via \ref createValue.
        self.value = {}

        ## Register cleanup function to guarantee cleanup of the capsProblem
        ## See: https://stackoverflow.com/questions/16333054/what-are-the-implications-of-registering-an-instance-method-with-atexit-in-pytho
        def cleanup(): 
            if <void*>self.problemObj != NULL:
                self.closeCAPS()

        atexit.register(cleanup)
            
    
    # During garbage college a call to closeCAPS is made
    # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html
    def __dealloc__(self):
        # Cannot call python functions in __dealloc__
        if <void*>self.problemObj != NULL:
            cCAPS.caps_close(self.problemObj)

        self.problemObj = NULL 
    
    # __del__ is not used in cython extension
    # https://cython.readthedocs.io/en/latest/src/userguide/special_methods.html
    #def __del__(self):
    #   self.closeCAPS()

    ## Loads a *.csm, *.caps, or *.egads file into the problem. 
    # See \ref problem1.py, \ref problem2.py, and \ref problem8.py for example use cases.
    #
    # \param capsFile CAPS file to load. Options: *.csm, *.caps, or *.egads. If the filename has a *.caps extension 
    # the pyCAPS analysis, bound, and value objects will be re-populated (see remarks).
    #  
    # \param projectName CAPS project name. projectName=capsFile if not provided.
    # \param verbosity Level of output verbosity. See \ref setVerbosity . 
    #
    # \return Optionally returns the reference to the \ref geometry class object (\ref pyCAPS.capsGeometry). 
    #
    # \remark 
    # Caveats of loading an existing CAPS file:
    # - Can currently only load *.caps files generated from pyCAPS originally.
    # - OpenMDAO objects won't be re-populated for analysis objects 
    def loadCAPS(self, capsFile, projectName=None, verbosity=None):
        ## \example problem1.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() function in which are multiple 
        # problems with geometry are created.
        
        ## \example problem2.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() function - compare "geometry" attribute with 
        # returned geometry object.
        
        ## \example problem8.py
        # Example use case for the pyCAPS.capsProblem.loadCAPS() and pyCAPS.capsProblem.saveCAPS() functions - using a CAPS file
        # to initiate a new problem.

        capsGeometry(self, capsFile, projectName)
        
        if verbosity is not None:
            self.setVerbosity(verbosity)
            
        return self.geometry
    
    ## Set the verbosity level of the CAPS output. 
    # See \ref problem5.py for a representative use case.
    #
    # \param verbosityLevel Level of output verbosity. Options: 0 (or "minimal"), 
    # 1 (or "standard") [default], and 2 (or "debug").
    def setVerbosity(self, verbosityLevel):
        ## \example problem5.py 
        # Basic example for setting the verbosity of a problem using pyCAPS.capsProblem.setVerbosity() function. 
        
        verbosity = {"minimal"  : 0, 
                     "standard" : 1, 
                     "debug"    : 2}
        
        if not isinstance(verbosityLevel, (int,float)):
            if verbosityLevel.lower() in verbosity.keys():
                verbosityLevel = verbosity[verbosityLevel.lower()]
            
        if int(verbosityLevel) not in verbosity.values():
            self.status = cCAPS.CAPS_BADVALUE
            self.checkStatus(msg = "while setting verbosity, invalid verbosity level!")
            
        self.status = cCAPS.caps_outLevel(self.problemObj, 
                                          <int > verbosityLevel)
        if self.status < 0:
            self.checkStatus(msg = "while setting verbosity" 
                                 + " (during a call to caps_outLevel)")
            
    ## Save a CAPS problem.
    # See \ref problem8.py for example use case. 
    # \param filename File name to use when saving CAPS problem.  
    def saveCAPS(self, filename="saveCAPS.caps"):
        
        if filename is "saveCAPS.caps":
            print("Using default file name - " + str(filename))

        file_path, file_extension = os.path.splitext(filename)
        if ".caps" not in file_extension:
            filename += ".caps"
        
        if os.path.isfile(filename):
            print("Warning: "+ str(filename) + " will be overwritten!")
            os.remove(filename)
            
        filename = _byteify(filename)
        self.status = cCAPS.caps_save(self.problemObj, 
                                      <const char *> filename)
        self.checkStatus(msg = "while saving CAPS")

    ## Close a CAPS problem.
    # See \ref problem1.py for a representative use case.
    def closeCAPS(self):
        
        self.status = cCAPS.caps_close(self.problemObj)
        self.checkStatus(msg = "while closing CAPS")
        
        self.problemObj = NULL 
    
    ## Report what analyses loaded into the problem are dirty.
    #
    # \return Optionally returns a list of names of the dirty analyses. An empty list is returned if no
    # analyses are dirty.
    def dirtyAnalysis(self):
      
        numDirty = 0
        dirtyAnalysis = []
        
        for i in self.analysis:
            if self.analysis[i].getAnalysisInfo(printInfo=False) > 0:
                numDirty = numDirty + 1
                dirtyAnalysis.append(self.analysis[i].aimName)
        
        print("Number of dirty analyses = ", numDirty)
        if len(dirtyAnalysis) > 0:
            print("Dirty analyses = ", dirtyAnalysis)
        
        return dirtyAnalysis

    ## Load an AIM (Analysis Interface Module) into the problem. See 
    # examples \ref problem3.py and  \ref problem4.py for typical representative use cases. 
    # 
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param aim Name of the requested AIM. 
    #
    # \param altName Alternative name to use when referencing AIM inside 
    # the problem (dictionary key in \ref analysis). The name of the AIM, aim, will be used if no 
    # \ref capsProblem.loadAIM$altName is provided (see remarks).
    #
    # \param analysisDir  Directory for AIM analysis. If none is provided the directory of the last loaded AIM
    # will be used; if no AIMs have been load the current working directory is used. CAPS requires that an unique 
    # directory be specified for each instance of AIM.
    #
    # \param capsIntent Analysis intention in which to invoke the AIM. 
    #
    # \param parents Single or list of parent AIM names to initilize the AIM with.
    #
    # \param copyAIM Name of AIM to copy. Creates the new AIM instance by duplicating an existing AIM. Analysis directory 
    # (\ref loadAIM$analysisDir) and \ref loadAIM$altName should be provided and different 
    # from the AIM being copied. See example \ref analysis2.py for a representative use case.
    #
    # \param copyParents When copying an AIM, should the same parents also be used (default - True).  
    #
    # \return Optionally returns the reference to the analysis dictionary (\ref analysis) entry created 
    # for the analysis class object (\ref pyCAPS.capsAnalysis).
    #
    # \remark If no \ref capsProblem.loadAIM$altName is provided and an AIM with the name, \ref capsProblem.loadAIM$aim, 
    # has already been loaded, an alternative 
    # name will be automatically specified with the syntax \ref capsProblem.loadAIM$aim_[instance number]. If an  
    # \ref capsProblem.loadAIM$altName is provided it must be unique compared to other instances of the loaded AIM. 
    def loadAIM(self, **kwargs):
        ## \example problem3.py
        # Example use cases for the pyCAPS.capsProblem.loadAIM() function.
        
        ## \example problem4.py
        # Example use cases for the pyCAPS.capsProblem.loadAIM() function - compare analyis dictionary entry
        # with returned analyis object.
        
        ## \example analysis2.py 
        # Duplicate an AIM using the \ref pyCAPS.capsProblem.loadAIM()$copyAIM keyword argument of the pyCAPS.capsProblem.loadAIM() function.
        
        # Load the aim
        return capsAnalysis(self, **kwargs)

    ## Alteranative to \ref createDataBound.
    # Enforces that at least 2 AIMs must be already loaded into the problem. See \ref createDataBound 
    # for details.
    def createDataTransfer(self, **kwargs):
        
        if len(self.analysis.keys()) < 2:
            print("At least two AIMs need to be loaded before a data transfer can be setup")
            self.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while creating a data transfer")
        
        return self.createDataBound(**kwargs)

    ## Create a CAPS data bound/transfer into the problem.
    #
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param capsBound Name of capsBound to use for the data bound. 
    # \param variableName Single or list of variables names to add. 
    # \param aimSrc Single or list of AIM names that will be the data sources for the bound.
    # \param aimDest Single or list of AIM names that will be the data destinations during the transfer. 
    # \param transferMethod Single or list of transfer methods to use during the transfer.
    # \param initValueDest Single or list of initial values for the destication data.
    #
    # \return Optionally returns the reference to the data bound dictionary (\ref dataBound) entry created 
    # for the bound class object (\ref pyCAPS.capsBound).
    def createDataBound(self, **kwargs):
        
        return capsBound(self, **kwargs)

    
    ## Create a CAPS value object. Only a value of subtype in PARAMETER is currently supported. 
    # See \ref value.py for a representative use case. Value objects are stored the \ref value dictionary.
    #
    # \param name Name used to define the value. This will be used as the keyword entry in the \ref value 
    # dictionary.
    # 
    # \param data Data value(s) for the variable. Note that type casting in done automatically based on the 
    # determined type of the Python object. Be careful with the Integers and Floats/Reals, for example 10 would be 
    # type casted as an Integer, while 10.0 would be a float - this small discrepancy may lead to type errors when 
    # linking values to analysis inputs. 
    # 
    # \param units Units associated with the value (default - None). 
    #
    # \param limits Valid/acceptable range for the value (default - None).
    # 
    # \param fixedLength Should the length of the value object be fixed (default - True)? For example if 
    # the object is initialized with a value of [1, 2] it can not be changed to [1, 2, 3]. 
    #
    # \param fixedShape Should the shape of the value object be fixed (default - True)? For example if 
    # the object is initialized with a value of 1 it can not be changed to [1, 2] or [[1, 2, 3], [4, 5, 6]] can 
    # not be changed to [[1, 2], [4, 5]]. 
    #
    # \return Optionally returns the reference to the value dictionary (\ref value) entry created 
    # for the value (\ref pyCAPS._capsValue).
    def createValue(self, name, data, units=None, limits=None, fixedLength=True, fixedShape=True):
        ## \example value.py
        # Example use cases for pyCAPS.capsProblem.createValue() function.
        
        if name in self.value.keys():
            print("A value with the name,", name, "already exists!" )
            self.status = cCAPS.CAPS_BADNAME
            self.checkStatus()
            
        self.value[name] = _capsValue(self, name, data, 
                                      units=units, limits=limits, 
                                      subType="Parameter", 
                                      fixedLength=fixedLength, fixedShape=fixedShape)
        
        return self.value[name]
    
    ## Create a link between a created CAPS value parameter and analyis inputs of <b>all</b> loaded AIMs, automatically. 
    # Valid CAPS value, parameter objects must be created with \ref createValue(). Note, only links to 
    # ANALYSISIN inputs are currently made at this time. 
    # See \ref value6.py for a representative use case.
    #
    # \param value Value to use when creating the link (default - None). A combination (i.e. a single or list)
    # of \ref value dictionary entries and/or value object instances (returned from a 
    # call to \ref createValue()) can be used. If no value is provided, all entries in the 
    # value dictionary (\ref value) will be used.   
    def autoLinkValue(self, value=None):
        ## \example value6.py
        # Example use cases for pyCAPS.capsProblem.autoLinkValue() function.
        
        cdef: 
            cCAPS.capsObj tempObj
            cCAPS.capsObj aObj, vObj
        
        # If none we want to try to link all values
        if value is None:
            value = list(self.value.keys())
        
        # If not a list convert it to one 
        if not isinstance(value, list):
            value = [value]
        
        # Loop through the values 
        for i in value:
            
            # If it is a capsValue instance grab the name 
            if isinstance(i, _capsValue):
               i = i.name 
           
            # Check to make sure name is in value dictionary 
            if i not in self.value.keys():
                print("Unable to find", i, "value dictionary!")
                self.status = cCAPS.CAPS_NOTFOUND
                self.checkStatus(msg = "while autolinking values in problem!")
            
            varname = _byteify(self.value[i].name) # Byteify the value name 
            
            found = False
            # Loop through analysis dictionary  
            for analysisName in self.analysis.keys():
                
                aObj = (<capsAnalysis> self.analysis[analysisName]).__analysisObj()
                # Look to see if the value name is present as an analysis input variable
                self.status = cCAPS.caps_childByName(aObj,
                                                     cCAPS.VALUE,
                                                     cCAPS.ANALYSISIN,
                                                     <char *> varname,
                                                     &tempObj)
                
                # If not found skip to the next 
                if self.status == cCAPS.CAPS_NOTFOUND:
                    continue
                elif self.status != cCAPS.CAPS_SUCCESS:
                    self.checkStatus(msg = "while getting analysis variable - " + str(varname) 
                                     + " (during a call to caps_childByName)")
        
                vObj = (<_capsValue> self.value[i]).__valueObj()
                # Try to make the link if we made it this far
                self.status = cCAPS.caps_linkValue(vObj, cCAPS.Copy, tempObj)
                self.checkStatus(msg = "while autolink value - " + str(varname) + " in " + str(analysisName)
                                 + " (during a call to caps_linkValue)")
                
                print("Linked", str(self.value[i].name), "to analysis", str(analysisName), "input", _strify(varname))
                found = True
                
            if not found:
                print("No linkable data found for", str(self.value[i].name))
                
    ## Add an attribute (that is meta-data) to the problem object. See example 
    # \ref problem7.py for a representative use case.
    # 
    # \param name Name used to define the attribute.
    # 
    # \param data Data value(s) for the attribute. Note that type casting in done 
    # automatically based on the determined type of the Python object.
    def addAttribute(self, name, data):
        ## \example problem7.py
        # Example use cases for interacting the 
        # pyCAPS.capsProblem.addAttribute() and 
        # pyCAPS.capsProblem.getAttribute() functions.
        
        valueObj = _capsValue(self, name, data, subType="User")
            
        self.status = cCAPS.caps_setAttr(self.problemObj, 
                                         NULL, 
                                         valueObj.valueObj)
        self.checkStatus(msg = "while adding attribute to problem!")
        
        # valueObj should be automatically deleted after the function loses scopes
        # del valueObj # not neccessary
    
    ## Get an attribute (that is meta-data) that exists on the problem object. See example 
    # \ref problem7.py for a representative use case.
    # \param name Name of attribute to retrieve. 
    # \return Value of attribute "name" 
    def getAttribute(self, name): 
        cdef:
            cCAPS.capsObj valueObj
            
        name = _byteify(name)
        self.status = cCAPS.caps_attrByName(self.problemObj,
                                            <char *> name, 
                                            &valueObj)
        self.checkStatus(msg = "while getting attribute for problem!")
        
        return _capsValue(self, None, None).__setupObj(valueObj).value
    
    def __getAIMParents(self, parentList):
        """
        Check the requested parent list against the analysis dictionary
        """
        
        for i in parentList:
            if i in self.analysis.keys():
                pass
            else:
                print("Requested parent", i, "has not been initialized/loaded")
                return False

        return True
    
    def __repopulateProblem(self):
        
        cdef: 
            int numObj
            cCAPS.capsObj tempObj 
        
        # Repopulate analysis class    
        self.status = cCAPS.caps_size(self.problemObj, 
                                      cCAPS.ANALYSIS, cCAPS.NONE, 
                                      &numObj)
        self.checkStatus(msg = "while repopulating analysis in capsProblem ")
        
        for i in range(numObj):
            self.status = cCAPS.caps_childByIndex(self.problemObj, 
                                                  cCAPS.ANALYSIS, cCAPS.NONE, 
                                                  i+1, 
                                                  &tempObj)
            self.checkStatus(msg = "while repopulating analysis in capsProblem")
            
            tempAnalysis = capsAnalysis(self, __fromSave=True).__setupObj(tempObj) 
          
            self.analysis[tempAnalysis.aimName] = tempAnalysis
        
            # Set other analysis related variable in the problem 
            #self.aimGlobalCount += 1
            #self.capsIntent  = self.analysis[tempAnalysis.aimName].capsIntent
            #self.analysisDir = self.analysis[tempAnalysis.aimName].analysisDir

        # Repopulate value class 
        self.status = cCAPS.caps_size(self.problemObj, 
                                      cCAPS.VALUE, cCAPS.PARAMETER, 
                                      &numObj)
        self.checkStatus(msg = "while repopulating values in capsProblem")
        
        for i in range(numObj):
            self.status = cCAPS.caps_childByIndex(self.problemObj, 
                                                  cCAPS.VALUE, cCAPS.PARAMETER,
                                                  i+1, 
                                                  &tempObj)
            self.checkStatus(msg = "while repopulating values in capsProblem")
            
            tempValue =  _capsValue(self, None, None).__setupObj(tempObj)
        
            self.value[tempValue.name] = tempValue
        
        
        # Repopulate databound class 
        self.status = cCAPS.caps_size(self.problemObj, 
                                      cCAPS.BOUND, cCAPS.NONE, 
                                      &numObj)
        self.checkStatus(msg = "while repopulating dataBound in capsProblem")
        
        for i in range(numObj):
            print("Unable to currently repopulate dataBound from a *.caps file")
            self.status = cCAPS.CAPS_BADVALUE
            self.checkStatus(msg = "while repopulating dataBound in capsProblem")
   
            #     BOUND, VERTEXSET, DATASET
            #self.status = cCAPS.caps_childByIndex(self.problemObj, 
            #                                       cCAPS.BOUND, cCAPS.NONE, 
            #                                      i+1, 
            #                                      &tempObj)
            #self.checkStatus(msg = "while repopulating  values in capsProblem")
            
            #tempValue =  _capsValue(self, None, None).__setupObj(tempObj)
        
            #self.dataBound[tempValue.name] = tempValue
        
    ## Create a HTML dendrogram/tree of the current state of the problem. See 
    # example \ref problem6.py for a representative use case. The HTML file relies on 
    # the open-source JavaScript library, D3, to visualize the data. 
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file. 
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to 
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet 
    # access will be necessary to view the tree.  
    #
    # \param filename Filename to use when saving the tree (default - "myProblem"). Note an ".html" is 
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False 
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access 
    # will be necessary to view the tree.   
    #
    # \param analysisGeom Show the geometry for each analysis entity (default - False).
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    #
    # \param reverseMap Reverse the geometry attribute map (default - False).
    def createTree(self, filename = "myProblem", **kwargs):
        ## \example problem6.py
        # Example use case for the pyCAPS.capsProblem.createTree() function.
        
        embedJSON = kwargs.pop("embedJSON", True)
        
        if embedJSON is False:
            try:
                # Write JSON data file 
                fp = open(filename + ".json",'w')
                fp.write( createTree(self, 
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                    )
                         )
                
                fp.close()
                
                # Write HTML file
                self.status = writeTreeHTML(filename +".html", 
                                            filename + ".json", 
                                            None, 
                                            kwargs.pop("internetAccess", True))
            except:
                self.status = cCAPS.CAPS_IOERR
        else: 
            self.status = writeTreeHTML(filename+".html", 
                                        None, 
                                        createTree(self, 
                                                   kwargs.pop("analysisGeom", False),
                                                   kwargs.pop("internalGeomAttr", False),
                                                   kwargs.pop("reverseMap", False)
                                                   ),
                                        kwargs.pop("internetAccess", True)
                                        )
        
        self.checkStatus(msg = "while creating HTML tree")
        
    def checkStatus(self, msg=None):
        # Check status flag returned by CAPS
        
        #def statusContinueTimer():
        #    print".....press Ctrl-C if you wish to proceed (you have 5 seconds)."
        #    try:
        #       time.sleep(5)
        #    except KeyboardInterrupt:
        #        return False
        #    
        #    return True
        
        if self.status == cCAPS.CAPS_SUCCESS:
            return 
        
        if self.raiseException is False: # If we have an error and have raiseException set to False
            
            print("Warning: Error detected, but raiseException is disabled - this may have unexepected consequences!")
            return 
        
        raise CAPSError(self.status, msg=msg)
                
    # Sanitize a filename
    def _createFilename(self, filename, directory, extension):
    
        if filename is None:
            self.status = cCAPS.CAPS_BADVALUE
            self.checkStatus(msg = "while creating a filename. Filename is None!")
         
        # Check to see if directory exists
        if not os.path.isdir(directory):
            print("Directory does not currently exist while creating a filename - it will" +
                    " be made automatically")
            os.makedirs(directory)
         
        filename = os.path.join(directory, filename)
        
        if "." not in extension:
            extension = "." + extension
        
        # only add an extension if the filename is missing one
        file_path, file_extension = os.path.splitext(filename)
        if file_extension == "":
            filename = filename+extension
    
        # Check to see if file exists
        try:
            os.remove(filename)
            print("File already exists - file " + filename + " will be deleted")
        except OSError:
           pass
       
        return filename
    
    cdef __getBodyName(self, object ibody, object nameList, cEGADS.ego body):
        cdef:
            # Attributes
            int atype
            int alen
            const int *ints
            const double *reals 
            const char *string
            
        name = "Body_" + str(ibody+1)
        
        attr = _byteify("_name")
        # Get body _name - if available 
        self.status = cEGADS.EG_attributeRet(body, <char *>attr, &atype, &alen, &ints, &reals, &string)
        if self.status == cEGADS.EGADS_NOTFOUND:
            pass
        elif self.status != cEGADS.EGADS_SUCCESS:
            self.checkStatus(msg= "while getting trying to retrieve body name (during a call to EG_attributeRet)")
            
        else:
            if atype == cEGADS.ATTRSTRING:
                name = _strify(string)
        if name:
            if name in nameList:
                name = "Body_" + str(ibody+1)
                print("Body", ibody, ": Name '", name, "' found more than once. Changing name to: '", name, "'")
                
            nameList.append(name)
            
        return name, nameList

## Defines a CAPS analysis object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.loadAIM .
cdef class capsAnalysis:

    cdef:
        readonly capsProblem capsProblem
        
        cCAPS.capsObj analysisObj
    
        readonly object aimName
        readonly object officialName
        readonly object analysisDir
        readonly object parents
        readonly object capsIntent
        readonly object unitSystem
        
        public object openMDAOComponent

    ## Initializes an analysis object. See \ref capsProblem.loadAIM for additional information and caveats. 
    #
    # \param problemObject capsProblem object in which the analysis object will be added to. 
    #
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param aim Name of the requested AIM. 
    #
    # \param altName Alternative name to use when referencing AIM inside 
    # the problem (dictionary key in \ref capsProblem.analysis). The name of the AIM, aim, will be used if no 
    # alternative name is provided (see remarks).
    #
    # \param analysisDir  Directory for AIM analysis. If none is provided the directory of the last loaded AIM
    # will be used; if no AIMs have been load the current working directory is used. CAPS requires that an unique 
    # directory be specified for each instance of AIM.
    #
    # \param capsIntent Analysis intention in which to invoke the AIM. 
    #
    # \param parents Single or list of parent AIM names to initilize the AIM with.
    #
    # \param copyAIM Name of AIM to copy. Creates the new AIM instance by duplicating an existing AIM. Analysis directory 
    # (\ref capsProblem.loadAIM$analysisDir) and \ref capsProblem.loadAIM$altName should be provided and different 
    # from the AIM being copied. See example \ref analysis2.py for a representative use case.
    #
    # \param copyParents When copying an AIM, should the same parents also be used (default - True).  
    #
    # \remark If no $altName is provided and an AIM with the name, aim, 
    # has already been loaded into the problem, an alternative 
    # name will be automatically specified with the syntax $aim_[instance number]. If an  
    # $altName is provided it must be unique compared to other instances of the loaded AIMs in the problem. 
    def __cinit__(self, problemObject, **kwargs):
        
        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError("Provided problemObject is not a capsProblem or CapsProblem class")
        
        ## Reference to the problem object that loaded the AIM during a call to \ref capsProblem.loadAIM .
        self.capsProblem = problemObject
        
        # Check the requested parent list against the analysis dictionary
        def checkAIMParents(parentList):
            for i in parentList:
                if i in self.capsProblem.analysis.keys():
                    pass
                else:
                    print("Requested parent", i, "has not been initialized/loaded")
                    return False

            return True
    
        # Check the analysis dictionary for redundant names - need unique names
        def checkAnalysisDict(name, printInfo=True):
            if name in self.capsProblem.analysis.keys():
                if printInfo:
                    print("\nAnalysis object with the name, " + str(name) +  ", has already been initilized")
                return True
            
            else:
                return False
        # Check the analysis dictionary for redundant directories if official aim names are the same
        def checkAnalysisDir(name, analysidDir):
            
            for aim in self.capsProblem.analysis.values():
                
                if name == aim.officialName and analysidDir == aim.analysisDir:
                    print("\nAnalysis object with the name, ", str(name),  ", and directory, ", str(analysidDir),
                          ", has already been initilized - Please use a different directory!!\n")
                    return True
            
            return False
        
        # Check to see if an AIM has been fully specified
        def checkAIMStatus(analysisDir, currentAIM, copyAIM):
            if copyAIM is not None:
                if analysisDir is None:
                    print("\nError: Unable to copy AIM, analysisDir (analysis directory) have not been set!\n")
                    return False
                else: 
                    return True
                
            if analysisDir is None or currentAIM is None:
                print("\nError: Unable to load AIM, either the analysisDir (analysis directory) and/or aim variables have not been set!\n")
                return False
            else:
                return True  
        
        if "capsFidelity" in kwargs.keys():
            print("Warning: IMPORTANT - The term capsFidelity is no longer supported please change to capsIntent in both your *.csm and *.py files")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while loading aim")

        if self.capsProblem.geometry == None:
            print("\nError: A *.csm, *.caps, or *.egads file has not been loaded!")
            self.capsProblem.status = cCAPS.CAPS_NULLOBJ
            self.capsProblem.checkStatus(msg = "while loading aim")

            
        currentAIM = kwargs.pop("aim", None)
        capsIntent = kwargs.pop("capsIntent", "")
        analysisDir = kwargs.pop("analysisDir", None)
        parents = kwargs.pop("parents", [])
        copyAIM = kwargs.pop("copyAIM", None)
        copyParents = kwargs.pop("copyParent", True)
        
        __fromSave = kwargs.pop("__fromSave", False)
        
        if capsIntent:            
            if not isinstance(capsIntent, list):
                capsIntent = [capsIntent] 
            capsIntent = ';'.join(capsIntent)

        if not analysisDir:
            print("No analysis directory provided - defaulting to current working directory")
            analysisDir = os.getcwd()

        if parents:
            if not isinstance(parents, list):
                parents = [parents] 
            
            temp = []
            for i in parents:
                
                if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                    temp.append(i.aimName)
                elif isinstance(i, str):
                    temp.append(i)
                else:    
                    raise TypeError("Provided parent is not a capsAnalysis or CapsAnalysis class or a string name")   
            
            # Update parents 
            parents = temp[:]
           
            # Get parent index from anaylsis dictionary
            if not checkAIMParents(parents):
                print("Unable to find parents")
                self.capsProblem.status = cCAPS.CAPS_NULLOBJ
                self.capsProblem.checkStatus(msg = "while loading aim")

    
        if copyAIM:
            if copyAIM not in self.capsProblem.analysis.keys():
                print("\n Error: AIM " + str(copyAIM) + " to duplicate has not been loaded into the problem!")
                
                self.capsProblem.status = cCAPS.CAPS_BADNAME
                self.capsProblem.checkStatus(msg = "while loading aim")
            else: 
                currentAIM = self.capsProblem.analysis[copyAIM].officialName 
                capsIntent = self.capsProblem.analysis[copyAIM].capsIntent
            
            if copyParents == True:
                parents = self.capsProblem.analysis[copyAIM].parents 
                 
        altName = kwargs.pop("altName", currentAIM)
        
        if not __fromSave:
            # Check aim status
            if not checkAIMStatus(analysisDir, currentAIM, copyAIM):
                self.capsProblem.status = cCAPS.CAPS_BADVALUE
                self.capsProblem.checkStatus(msg = "while loading aim - " + currentAIM)
            
            # Check analysis dictionary for redunant names 
            if checkAnalysisDict(altName):
                
                if altName == currentAIM: 
                    i = 0
                    while i < 1000:
                        i +=1
                        
                        altName = currentAIM + "_" + str(i)
                        if not checkAnalysisDict(altName, printInfo=False):
                            break
                    else:
                        print("Unable to auto-specify altName, are there more than 1000 AIM instances loaded?")
                        self.capsProblem.status = cCAPS.CAPS_BADNAME
                        self.capsProblem.checkStatus(msg = "while loading aim - " + currentAIM)
                    
                    print("An altName of", altName, "will be used!")
                    
                else:
                    print("The altName,", altName, ", has already been used!")
                    self.capsProblem.status = cCAPS.CAPS_BADNAME
                    self.capsProblem.checkStatus(msg = "while loading aim - " + currentAIM)
            
            # Check analysis for redunant directories if aim has already been loaded
            if checkAnalysisDir(currentAIM, analysisDir):
                self.capsProblem.status = cCAPS.CAPS_DIRERR
                self.capsProblem.checkStatus(msg = "while loading aim - " + currentAIM)
        
        self.__createObj(altName, currentAIM, analysisDir, capsIntent,
                         unitSystem=kwargs.pop("unitSystem", None),
                         parentList=parents, 
                         copyAIM=copyAIM)
        
        # Load the aim
        self.capsProblem.analysis[altName] = self
    
    cdef __createObj(self, 
                     name, officialName, analysisDir, capsIntent, # Use if creating an analysis from scratch, else set to None and call __setupObj
                     unitSystem=None, 
                     parentList=None, 
                     copyAIM=None):
        
        cdef:
            cCAPS.capsObj *aobjTemp = NULL
            cCAPS.capsObj aobjCopy = NULL
            char *intent = NULL
            
        ## Reference name of AIM loaded for the analysis object (see \ref capsProblem.loadAIM$altName).  
        self.aimName = name
        
        ## Name of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$aim).
        self.officialName = officialName
        
        ## Analysis directory of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$analysis). 
        # If the directory does not exist it will be made automatically.
        self.analysisDir = analysisDir
        
        ## Unit system the AIM was loaded with (if applicable). 
        self.unitSystem = unitSystem
        
        ## List of parents of the AIM loaded for the anlaysis object (see \ref capsProblem.loadAIM$parents).
        self.parents = parentList
        
        ## Analysis intent of the AIM loaded for the analysis object (see \ref capsProblem.loadAIM$capsIntent
        self.capsIntent = capsIntent 
        
        ## OpenMDAO "component" object (see \ref createOpenMDAOComponent for instantiation).    
        self.openMDAOComponent = None
        
        # It is assumed we are going to call __setupObj
        if self.aimName is None:
            return 
        
        # Check to see if directory exists
        if not os.path.isdir(self.analysisDir):
            print("Analysis directory does not currently exist - it will be made automatically")
            os.makedirs(self.analysisDir)
        
        # Allocate a temporary list of analysis object for the parents
        if len(parentList) != 0:

            aobjTemp = <cCAPS.capsObj *> malloc(len(self.parents)*sizeof(cCAPS.capsObj))
            
            if not aobjTemp:
                raise MemoryError()
            else:
                for i, parent in enumerate(self.parents):
                    aobjTemp[i] = (<capsAnalysis> self.capsProblem.analysis[parent]).__analysisObj()
                
        officialName = _byteify(officialName)
        analysisDir = _byteify(self.analysisDir)
        unitSys = _byteify(self.unitSystem)
        
        if capsIntent:
            capsIntent = _byteify(capsIntent)
            intent = capsIntent
        
        
        if copyAIM is not None:
            aobjCopy = (<capsAnalysis> self.capsProblem.analysis[copyAIM]).__analysisObj()
            self.capsProblem.status = cCAPS.caps_dupAnalysis(aobjCopy,
                                                             <char *> analysisDir,
                                                             <int>len(self.parents),
                                                             aobjTemp,
                                                             &self.analysisObj)
        else:
            if unitSys:
                self.capsProblem.status = cCAPS.caps_load(self.capsProblem.problemObj,
                                                          <char *> officialName,
                                                          <char *> analysisDir,
                                                          <char *> unitSys,
                                                          intent,
                                                          <int>len(self.parents),
                                                          aobjTemp,
                                                          &self.analysisObj)
            else:

                self.capsProblem.status = cCAPS.caps_load(self.capsProblem.problemObj,
                                                          <char *> officialName,
                                                          <char *> analysisDir,
                                                          NULL,
                                                          intent,
                                                          <int>len(self.parents),
                                                          aobjTemp,
                                                          &self.analysisObj)
        if aobjTemp != NULL:
            free(aobjTemp)
            
        if copyAIM is not None:
            self.capsProblem.checkStatus(msg = "while loading aim - " + str(self.aimName) 
                                             + " (during a call to caps_dupAnalysis)")
        else:
            self.capsProblem.checkStatus(msg = "while loading aim - " + str(self.aimName) 
                                             + " (during a call to caps_load)")
        
        # Add attribute information for saving the problem
        self.addAttribute("__aimName"     , self.aimName)
        self.addAttribute("__officialName", self.officialName)
        if self.parents:
            self.addAttribute("__parents", self.parents)
    
    
    cdef __setupObj(self, cCAPS.capsObj analysisObj):
        
        self.analysisObj = analysisObj
            
        self.aimName = self.getAttribute("__aimName")
        self.officialName = self.getAttribute("__officialName")
        try:
            self.parents = self.getAttribute("__parents")
            
        except CAPSError as e:
            
            if e.errorName != "CAPS_NOTFOUND":
                self.capsProblem.status = e.errorName
                self.capsProblem.checkStatus(msg = "while repopulating aim - " + str(self.aimName))
            
            self.parents = None
                    
        analysisInfo = self.getAnalysisInfo(printInfo=False, infoDict=True)
        
        self.analysisDir = analysisInfo["analysisDir"]
        self.capsIntent = analysisInfo["intent"]
        self.unitSystem = analysisInfo["unitSystem"]
        
        return self
    
    cdef cCAPS.capsObj __analysisObj(self):
        return self.analysisObj
    
    ## Sets an ANALYSISIN variable for the analysis object.
    # \param varname Name of CAPS value to set in the AIM.
    # \param value Value to set. Type casting is automatically done based on the CAPS value type 
    # \param units Applicable units of the current variable (default - None). Only applies to real values.    
    # specified in the AIM.
    def setAnalysisVal(self, varname, value, units=None):
       
        cdef: 
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int dim
            cCAPS.capsFixed lfixed
            cCAPS.capsFixed sfixed
            cCAPS.capsNull nval
            int nrow
            int ncol

            int vlen
            int nErr
            const char *valunits
        
            void *data
             
        if varname is None:
            print("No variable name provided")
            return
        
        varname = _byteify(varname)
       
        self.capsProblem.status = cCAPS.caps_childByName(self.analysisObj,
                                                         cCAPS.VALUE,
                                                         cCAPS.ANALYSISIN,
                                                         <char *> varname,
                                                         &tempObj)
        self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(varname) 
                                         + " (during a call to caps_childByName)")
    
        self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                      &vtype,
                                                      &vlen,
                                                      <const void ** > &data,
                                                      &valunits,
                                                      &nErr,
                                                      &errors)
        # Check errors
        report_Errors(nErr, errors)
        self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(varname) 
                                         + " (during a call to caps_getValue)")
        
        # Special treatement for nullifying a value
        if value is None:
            self.capsProblem.status = cCAPS.caps_getValueShape(tempObj, &dim, &lfixed, &sfixed, &nval, &nrow, &ncol)
            self.capsProblem.checkStatus(msg = "while getting analysis variable shape - " + str(varname) 
                                             + " (during a call to caps_getValueShape)")

            self.capsProblem.status = cCAPS.caps_setValue(tempObj, 
                                                          nrow, 
                                                          ncol, 
                                                          <const void *> NULL)
            self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(varname) 
                                             + " to None (during a call to caps_setValue)")
            return


        self.capsProblem.status, value = convertUnitPythonObj(tempObj,
                                                              vtype,
                                                              _byteify(units),
                                                              value)
        self.capsProblem.checkStatus(msg = "while converting analysis variable - " + str(varname) 
                                         + " units (during a call to caps_convert)")
       
        numRow, numCol = getPythonObjShape(value)
            
        byteValue  = []
        data = castValue2VoidP(value, byteValue, vtype)
        
        try:
            self.capsProblem.status = cCAPS.caps_setValue(tempObj, 
                                                          <int> numRow, 
                                                          <int> numCol, 
                                                          <const void *> data)
            
            if self.capsProblem.status == cCAPS.CAPS_LINKERR:
                print("CAPS_LINKERR detected while trying to set analysis variable ",
                        str(varname), " - link will be broken in order to set value!!\n")

                # Break the link
                self.capsProblem.status = cCAPS.caps_linkValue(NULL, cCAPS.Copy, tempObj)
                #self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(varname) 
                #                         + " (during a call to caps_linkValue)")
                
                # Try to set the value again 
                self.capsProblem.status = cCAPS.caps_setValue(tempObj, 
                                                              <int> numRow, 
                                                              <int> numCol, 
                                                              <const void *> data)
        finally:
            if isinstance(byteValue,list):
                del byteValue[:]
            
            if data and vtype != cCAPS.String:
                free(data)
            
        self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(varname) 
                                         + " (during a call to caps_setValue)")

    ## Gets an ANALYSISIN variable for the analysis object.
    # \param varname Name of CAPS value to retrieve from the AIM. If no name is provided a dictionary
    # containing all ANALYSISIN values is returned. See example \ref analysis4.py for a representative 
    # use case.
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param units When set to True returns the units along with the value as specified by the analysis. When 
    # set to a string (e.g. units="ft") the returned value is converted into the specified units. (default - False).    
    #
    # \param namesOnly Return only a list of variable names (no values) if creating a dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values. Units are also returned if applicable 
    # based on the "units" keyword (does not apply to dictionary returns).
    def getAnalysisVal(self, varname = None, **kwargs):
        ## \example analysis4.py
        # Example use cases for interacting the 
        # pyCAPS.capsAnalysis.getAnalysisVal() and 
        # pyCAPS.capsAnalysis.getAnalysisOutVal() functions.
       
        cdef:
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int vlen
            int nErr
            const char *valunits
            const void *data

            int vdim, vnrow, vncol
            cCAPS.capsFixed vlfixed, vsfixed
            cCAPS.capsNull  vntype
        
            int index, numAnalysisInput # returning a dictionary variables 
            char * name
            
        analysisIn = {}
        
        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", False)
         
        if varname is not None: # Get just the 1 variable
            
            varname = _byteify(varname)
            self.capsProblem.status = cCAPS.caps_childByName(self.analysisObj,
                                                             cCAPS.VALUE,
                                                             cCAPS.ANALYSISIN,
                                                             <char *> varname,
                                                             &tempObj)
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname) 
                                             + " (during a call to caps_childByName)")
            
            self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                          &vtype,
                                                          &vlen,
                                                          &data,
                                                          &valunits,
                                                          &nErr,
                                                          &errors)
            # Check errors
            report_Errors(nErr, errors)
            
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname)  
                                             + " (during a call to caps_getValue)")
            
            self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                               &vdim,
                                                               &vlfixed,
                                                               &vsfixed,
                                                               &vntype,
                                                               &vnrow,
                                                               &vncol)
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname) 
                                             + " (during a call to caps_getValueShape)")
    
            
            if vntype == cCAPS.IsNull:
                
                if units:
                    
                    if isinstance(units, str):
                        print("Analysis variable, ", _strify(varname) , ", is None; unable to convert units to", units )
                        return None, units 
                    
                    if valunits:
                        return None, <object> _strify(<char *> valunits)
                    else:
                        return None, None
                        
                else:
                    return None
            else:
                value = castValue2PythonObj(data, vtype, vnrow, vncol)
                
                if units:
                    
                    if isinstance(units, str):
                        
                        if not valunits:
                            print("No units assigned to analysis variable, ", _strify(varname) , ", unable to convert units to", units )
                            return value, None
                        
                        value = capsConvert(value, <object> _strify(<char *> valunits), units)
                        if value == None:
                            self.capsProblem.status = cCAPS.CAPS_UNITERR
                            self.capsProblem.checkStatus(msg = "while converting analysis variable - " + _strify(varname)  
                                                         + " units (during a call to pyCAPS.capsConvert)")
                        
                        return value, units
                        
                    else:
                        if valunits:
                            return value, <object> _strify(<char *> valunits)
                        else:
                            return value, None
                else: 
                    return value
                   
        else: # Build a dictionary 
             
            self.capsProblem.status = cCAPS.caps_size(self.analysisObj,
                                                      cCAPS.VALUE,
                                                      cCAPS.ANALYSISIN,
                                                      &numAnalysisInput)
            self.capsProblem.checkStatus(msg = "while getting number of analysis inputs" 
                                             + " (during a call to caps_size)")
            
            for i in range(1, numAnalysisInput+1):
                self.capsProblem.status = cCAPS.caps_childByIndex(self.analysisObj,
                                                                  cCAPS.VALUE,
                                                                  cCAPS.ANALYSISIN,
                                                                  i,
                                                                  &tempObj)
                self.capsProblem.checkStatus(msg = "while getting analysis input variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_childByIndex)")
                 
                name = tempObj.name
                 
                self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                              &vtype,
                                                              &vlen,
                                                              &data,
                                                              &valunits,
                                                              &nErr,
                                                              &errors)
                # Check errors
                report_Errors(nErr, errors)
             
                self.capsProblem.checkStatus(msg = "while getting analysis input variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValue)")
             
                self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                                   &vdim,
                                                                   &vlfixed,
                                                                   &vsfixed,
                                                                   &vntype,
                                                                   &vnrow,
                                                                   &vncol)
                self.capsProblem.checkStatus(msg = "while getting analysis input variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValueShape)")
     
             
                if vntype == cCAPS.IsNull:
                    analysisIn[<object> _strify(name)] = None
                else:
                    analysisIn[<object> _strify(name)] = castValue2PythonObj(data, vtype, vnrow, vncol)
             
            if namesOnly is True:
                return analysisIn.keys()
            else:
                return analysisIn

    ## Gets an ANALYSISOUT variable for the analysis object.
    # \param varname Name of CAPS value to retrieve from the AIM. If no name is provided a dictionary
    # containing all ANALYSISOUT values is returned. See example \ref analysis4.py for a 
    # representative use case.
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param units When set to True returns the units along with the value as specified by the analysis. When 
    # set to a string (e.g. units="ft") the returned value is converted into the specified units. (default - False).    
    #
    # \param namesOnly Return only a list of variable names (no values) if creating a dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values. Units are also returned if applicable 
    # based on the "units" keyword (does not apply to dictionary returns).
    def getAnalysisOutVal(self, varname = None, **kwargs):
        ## \example analysis4.py
        # Example use cases for interacting the 
        # pyCAPS.capsAnalysis.getAnalysisVal() and 
        # pyCAPS.capsAnalysis.getAnalysisOutVal() functions.
       
        cdef:
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int vlen
            int nErr
            const char *valunits
            const void *data

            int vdim, vnrow, vncol
            cCAPS.capsFixed vlfixed, vsfixed
            cCAPS.capsNull  vntype
            int index, numAnalysisOutput # returning a dictionary variables 
            char * name
            
        analysisOut = {}
        
        namesOnly = kwargs.pop("namesOnly", False)
        units = kwargs.pop("units", False)
        
        stateOfAnalysis = self.getAnalysisInfo(printInfo=False)
        
        if stateOfAnalysis != 0 and namesOnly == False:
            print("Analysis state isn't clean (state = ", stateOfAnalysis, ") can only return names of ouput variables; set namesOnly keyword to True")
            self.capsProblem.status = cCAPS.CAPS_DIRTY
            self.capsProblem.checkStatus(msg = "while getting analysis out variables")
             
        if varname is not None:
            
            varname = _byteify(varname)
            self.capsProblem.status = cCAPS.caps_childByName(self.analysisObj,
                                                             cCAPS.VALUE,
                                                             cCAPS.ANALYSISOUT,
                                                             <char *> varname,
                                                             &tempObj)
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname) 
                                             + " (during a call to caps_childByName)")
            
            self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                          &vtype,
                                                          &vlen,
                                                          &data,
                                                          &valunits,
                                                          &nErr,
                                                          &errors)
            
            # Check errors
            report_Errors(nErr, errors)
            
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname) 
                                             + " (during a call to caps_getValue)")
            
            self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                               &vdim,
                                                               &vlfixed,
                                                               &vsfixed,
                                                               &vntype,
                                                               &vnrow,
                                                               &vncol)
            self.capsProblem.checkStatus(msg = "while getting analysis variable - " + _strify(varname) 
                                             + " (during a call to caps_getValueShape)")
    
            
            #self.capsProblem.status = castValue2PythonObj(data, vtype, vnrow, vncol, value)
            #self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
            #                                 + "(during value casting)")
            
            if vntype == cCAPS.IsNull:
                if units:
                    
                    if isinstance(units, str):
                        print("Analysis variable, ", _strify(varname) , ", is None; unable to convert units to", units )
                        return None, units 
                    
                    if valunits:
                        return None, <object> _strify(<char *> valunits)
                    else:
                        return None, None
                        
                else:
                    return None
            else:
                value =  castValue2PythonObj(data, vtype, vnrow, vncol)
                
                if units:
                    
                    if isinstance(units, str):
                        
                        if not valunits:
                            print("No units assigned to analysis variable, ", _strify(varname), ", unable to convert units to", units )
                            return value, None
                        
                        value = capsConvert(value, <object> _strify(<char *> valunits), units)
                        if value == None:
                            self.capsProblem.status = cCAPS.CAPS_UNITERR
                            self.capsProblem.checkStatus(msg = "while converting analysis variable - " + _strify(varname)  
                                                         + " units (during a call to pyCAPS.capsConvert)")
                        
                        return value, units
                        
                    else:
                        if valunits:
                            return value, <object> _strify(<char *> valunits)
                        else:
                            return value, None
                else: 
                    return value
                
        else: # Build a dictionary
            self.capsProblem.status = cCAPS.caps_size(self.analysisObj,
                                                      cCAPS.VALUE,
                                                      cCAPS.ANALYSISOUT,
                                                      &numAnalysisOutput)
            self.capsProblem.checkStatus(msg = "while getting number of analysis outputs" 
                                             + " (during a call to caps_size)")
            
            for i in range(1, numAnalysisOutput+1):
                self.capsProblem.status = cCAPS.caps_childByIndex(self.analysisObj,
                                                                  cCAPS.VALUE,
                                                                  cCAPS.ANALYSISOUT,
                                                                  i,
                                                                  &tempObj)
                self.capsProblem.checkStatus(msg = "while getting analysis output variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_childByIndex)")
                 
                name = tempObj.name
                
                if namesOnly is True:
                    analysisOut[<object> _strify(name)] = None
                else:
                    self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                                  &vtype,
                                                                  &vlen,
                                                                  &data,
                                                                  &valunits,
                                                                  &nErr,
                                                                  &errors)
                    # Check errors
                    report_Errors(nErr, errors)
                 
                    self.capsProblem.checkStatus(msg = "while getting analysis output variable (index) - " + str(i+1) 
                                                     + " (during a call to caps_getValue)")
                 
                    self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                                       &vdim,
                                                                       &vlfixed,
                                                                       &vsfixed,
                                                                       &vntype,
                                                                       &vnrow,
                                                                       &vncol)
                    self.capsProblem.checkStatus(msg = "while getting analysis output variable (index) - " + str(i+1) 
                                                     + " (during a call to caps_getValueShape)")
         
                 
                    if vntype == cCAPS.IsNull:
                        analysisOut[<object> _strify(name)] = None
                    else:
                        analysisOut[<object> _strify(name)] = castValue2PythonObj(data, vtype, vnrow, vncol)
             
            if namesOnly is True:
                return analysisOut.keys()
            else:
                return analysisOut
    
    ## Gets analysis information for the analysis object. See example \ref analysis6.py for 
    # a representative use case.
    #
    # \param printInfo Print information to sceen if True.
    # \param **kwargs See below.
    # 
    # \return Cleanliness state of analysis object or a dictionary containing analysis 
    # information (infoDict must be set to True)
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing analysis information instead 
    # of just the cleanliness state (default - False) 
    def getAnalysisInfo(self, printInfo=True, **kwargs):
        ## \example analysis6.py
        # Example use cases for interacting the 
        # pyCAPS.capsAnalysis.getAnalysisInfo() function.
         
        cdef:
            int i = 0
            int execution = 0
            char *intent = NULL
            int numParent = 0
            int numField = 0
            int *ranks = NULL
            int cleanliness = 0
            char *analysisPath = NULL
            char *unitSys = NULL
            char **fnames = NULL

            cCAPS.capsObj *aobjTemp = NULL
            
            object dirtyInfo = {0 : "Up to date", 
                                1 : "Dirty analysis inputs", 
                                2 : "Dirty geometry inputs", 
                                3 : "Both analysis and geometry inputs are dirty", 
                                4 : "New geometry", 
                                5 : "Post analysis required",
                                6 : "Execution and Post analysis required"}
            
        self.capsProblem.status = cCAPS.caps_analysisInfo(self.analysisObj, 
                                                          &analysisPath, 
                                                          &unitSys, 
                                                          &intent, 
                                                          &numParent, 
                                                          &aobjTemp, 
                                                          &numField, 
                                                          &fnames,
                                                          &ranks, 
                                                          &execution, 
                                                          &cleanliness)

        self.capsProblem.checkStatus(msg = "while getting analysis information - " + str(self.aimName)
                                         + " (during a call to caps_analysisInfo)")
        
        if intent:
            capsIntent = _strify(intent)
            if ";" in capsIntent: capsIntent = capsIntent.split(";")
        else:
            capsIntent = None
            
        if printInfo:
            print("Analysis Info:")
            print("\tName           = ", self.aimName)
            print("\tIntent         = ", capsIntent)
            print("\tAnalysisDir    = ", _strify(analysisPath))
            print("\tNumber of Parents = ", numParent)
            print("\tExecution Flag = ", execution)
            print("\tNumber of Fields = ", numField)
            
            if numField > 0:
                print("\tField name (Rank):")

            for i in range(numField):
                print("\t ", _strify(fnames[i]), "(", ranks[i], ")")
    
            if cleanliness in dirtyInfo:
                print("\tDirty state    = ", cleanliness, " - ", dirtyInfo[cleanliness] )
            else: 
                print("\tDirty state    = ", cleanliness, " ",)
        
        if execution == 1:
            executionFlag = True
        else:
            executionFlag = False 
        
        unitSystem = None
        if unitSys:
            unitSystem = _strify(<char *> unitSys)
             
        if kwargs.pop("infoDict", False):
            return {"name"     : _strify(self.aimName), 
                    "intent"   : capsIntent, 
                    "analysisDir"   : _strify(analysisPath),
                    "numParents"    : numParent, 
                    "executionFlag" : executionFlag, 
                    "fieldName"     : [_strify(fnames[i]) for i in range(numField)],
                    "fieldRank"     : [ranks[i] for i in range(numField)], 
                    "unitSystem"    : unitSystem,
                    "status"        : cleanliness
                    }
        
        return <object> cleanliness
    
    ## Retrieve a list of geometric attribute values of a given name 
    # ("attributeName") for the bodies loaded 
    # into the analysis. Level in which to search the bodies is 
    # determined by the attrLevel keyword argument. 
    # See \ref analysis3.py for a representative use case.
    #
    # \param attributeName Name of attribute to retrieve values for. 
    # \param **kwargs See below.
    # 
    # \return A list of attribute values.
    #
    # Valid keywords:
    # \param bodyIndex Specific body in which to retrieve attribute information from.
    # \param attrLevel Level to which to search the body(ies). Options: \n
    #                  0 (or "Body") - search just body attributes\n
    #                  1 (or "Face") - search the body and all the faces [default]\n
    #                  2 (or "Edge") - search the body, faces, and all the edges\n
    #                  3 (or "Node") - search the body, faces, edges, and all the nodes
    # 
    def getAttributeVal(self, attributeName, **kwargs):
        ## \example analysis3.py
        # Example use case of the pyCAPS.capsAnalysis.getAttributeVal() function
        cdef:
           int numBody = 0
            
           cEGADS.ego body # Temporary body container
            
           char *units = NULL
        
        # Initiate attribute value list     
        attributeList = []
        
        attributeName = _byteify(attributeName)
        
        bodyIndex = kwargs.pop("bodyIndex", None)
        attrLevel = kwargs.pop("attrLevel", "face")
        
        attribute = {"body" : 0, 
                     "face" : 1, 
                     "edge" : 2, 
                     "node" : 3}
        
        if not isinstance(attrLevel, (int,float)):
            if attrLevel.lower() in attribute.keys():
                attrLevel = attribute[attrLevel.lower()]
            
        if int(attrLevel) not in attribute.values():
            print("Error: Invalid attribute level! Defaulting to 'Face'")
            attrLevel = attribute["face"]
        
        # Force regeneration of the geometry
        self.capsProblem.geometry.buildGeometry()
        
        # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.analysisObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values for aim " 
                                         + str(self.aimName) + "(during a call to caps_size)")
        
    
        if numBody < 1:
            print("The number of bodies in the problem is less than 1!!!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values for aim " 
                                         + str(self.aimName))
        
        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i
            
            if bodyIndex is not None:
                if int(bodyIndex) != i+1:
                    continue
            
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.analysisObj, i+1, &body, &units)
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values (body " 
                                         + str(i+1) + " )"
                                         + " (during a call to caps_bodyByIndex)")
            
            self.capsProblem.status, returnList = createAttributeValList(body,
                                                                         <int> attrLevel, 
                                                                         <char *> attributeName)
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values (body " 
                                         + str(i+1) + " ) for aim " + str(self.aimName)
                                         + " (during a call to createAttributeValList)")
            attributeList = attributeList + returnList
        
        try:
            return list(set(attributeList)) 
        except TypeError:
            print("Unhashable type in attributelist")
            return attributeList
                    
    ## Create geometric attribution map (embeded dictionaries) for the bodies loaded 
    # into the analysis.
    # Dictionary layout: 
    # - Body 1 
    #   - Body : Body level attributes
    #   - Faces 
    #     - 1 : Attributes on the first face of the body
    #     - 2 : Attributes on the second face of the body
    #     - " : ... 
    #   - Edges 
    #     - 1 : Attributes on the first edge of the body
    #     - 2 : Attributes on the second edge of the body
    #     - " : ...
    #   - Nodes : 
    #     - 1 : Attributes on the first node of the body
    #     - 2 : Attributes on the second node of the body
    #     - " : ...
    # - Body 2
    #   - Body : Body level attributes 
    #   - Faces 
    #     - 1 : Attributes on the first face of the body
    #     - " : ...
    #   - ...
    # - ...
    #
    # Dictionary layout (reverseMap = True): 
    # - Body 1 
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not 
    #       - Faces : Face numbers at which the attribute exist
    #       - Edges : Edge numbers at which the attribute exist
    #       - Nodes : Node numbers at which the attribute exist
    #     -  Value : Next value of attribute with the same name
    #       - Body : True if value exist at body level, None if not
    #       -  "   : ...
    #     - ...
    #   - Atribute : Attribute name
    #     - Value : Value of attribute
    #       -   "  : ... 
    #     - ...
    # - Body 2
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not 
    #       -   "   : ...
    #     - ...
    #   - ...
    # - ...
    #
    # \param getInternal Get internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param reverseMap Reverse the attribute map (default - False). See above table for details.
    #
    # \return Dictionary containing attribution map
    def getAttributeMap(self, getInternal = False, **kwargs):
        
        cdef:
            int numBody = 0
           
            cCAPS.capsErrs *errors = NULL
            int nErr
            
            cEGADS.ego body # Temporary body container
            
            char *units = NULL
        
        attributeMap = {}
        
        # Force regeneration of the geometry
        self.capsProblem.geometry.buildGeometry()
        
        # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.analysisObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while creating attribute map"
                                          + " for aim " + str(self.aimName) 
                                          + " (during a call to caps_size)")
        
    
        if numBody < 1:
            print("The number of bodies in the problem is less than 1!!!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg =  "while creating attribute map"
                                          + " for aim " + str(self.aimName))
        
        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i
       
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.analysisObj, i+1, &body, &units)
            self.capsProblem.checkStatus(msg = "while creating attribute map on body " + str(i+1) 
                                        + " for aim " + str(self.aimName)
                                        + " (during a call to caps_bodyByIndex)")
            
            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}
            self.capsProblem.status, attributeMap["Body " + str(i + 1)] = createAttributeMap(body, getInternal)
            self.capsProblem.checkStatus(msg = "while creating attribute map on body " + str(i+1) 
                                                + " for aim " + str(self.aimName)
                                                + " (during a call to createAttributeMap)")
        if kwargs.pop("reverseMap", False):
            return reverseAttributeMap(attributeMap)
        else:
            return attributeMap 
     

    ## Alternative to \ref preAnalysis(). Warning: May be deprecated in future versions.
    def aimPreAnalysis(self):
        print("\n\nWarning: aimPreAnalysis is a deprecated. Please use preAnalysis.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.preAnalysis()
        
    ## Run the pre-analysis function for the AIM. If the specified analysis directory doesn't exist it 
    # will be made automatically. 
    def preAnalysis(self):
        
        cdef: 
            cCAPS.capsErrs *errors
            int nErr
        
        # Check to see if directory exists
        if not os.path.isdir(self.analysisDir):
            print("Analysis directory does not currently exist - it will be made automatically")
            os.makedirs(self.analysisDir)
            
        self.capsProblem.status = cCAPS.caps_preAnalysis(self.analysisObj,
                                                         &nErr,
                                                         &errors)
        # Check errors
        report_Errors(nErr, errors)
                                                         
        self.capsProblem.checkStatus(msg = "while getting aim pre-analysis - " + str(self.aimName)
                                         + " (during a call to caps_preAnalysis)")
    
    ## Alternative to \ref postAnalysis(). Warning: May be deprecated in future versions.
    def aimPostAnalysis(self):
        print("\n\nWarning: aimPostAnalysis is a deprecated. Please use postAnalysis.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.postAnalysis()
        
    ## Run the post-analysis function for the AIM.
    def postAnalysis(self):
         
        cdef: 
            char *name = NULL
            
            cCAPS.capsoType type
            cCAPS.capssType subtype 
            
            cCAPS.capsObj link
            cCAPS.capsObj parent
            cCAPS.capsOwn owner
            
            cCAPS.capsErrs *errors
            int nErr
        
        self.capsProblem.status = cCAPS.caps_info(self.capsProblem.problemObj, #self.analysisObj,
                                                  &name,
                                                  &type,
                                                  &subtype,
                                                  &link,
                                                  &parent,
                                                  &owner)

        self.capsProblem.checkStatus(msg = "while getting aim post-analysis - " + str(self.aimName)
                                         + " (during a call to caps_info)")

        self.capsProblem.status = cCAPS.caps_postAnalysis(self.analysisObj,
                                                          owner,
                                                          &nErr,
                                                          &errors)
        # Check errors
        report_Errors(nErr, errors)
                                                         
        self.capsProblem.checkStatus(msg = "while getting aim post-analysis - " + str(self.aimName)
                                         + " (during a call to caps_postAnalysis)")
    
    ## Create an OpenMDAO Component[1.7.3]/ExplicitComponent[2.8+] object; an external code component 
    # (ExternalCode[1.7.2]/ExternalCodeComp[2.8+]) is created if the 
    # executeCommand keyword arguement is provided. This functionality should work with either
    # verison 1.7.3 or >=2.8 of OpenMDAO.
    #
    # \param inputVariable Input variable(s)/parameter(s) to add to the OpenMDAO component. 
    # Variables may be either 
    # analysis input variables or geometry design parameters. Note, that the setting of analysis inputs 
    # supersedes the setting of geometry design parameters; issues may arise if analysis input and 
    # geometry design variables have the same name. 
    # If the analysis parameter wanting to be added to the OpenMDAO component is part of a
    # capsTuple the following notation should be used:
    #    "AnalysisInput:TupleKey:DictionaryKey", for example "AVL_Control:ControlSurfaceA:deflectionAngle" 
    # would correspond to the AVL_Control input variable, the ControlSurfaceA element of the input values
    #  (that is the name of the control surface being created) and finally deflectionAngle corresponds 
    # to the name of the dictionary entry that is to be used as the component parameter. 
    # If the tuple's value isn't a dictionary just "AnalysisInput:TupleKey" is needed.    
    #
    # \param outputVariable Output variable(s)/parameter(s) to add to the OpenMDAO component. Only scalar output variables
    # are currently supported 
    # 
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param changeDir Automatically switch into the analysis directory set for the AIM when executing an 
    # external code (default - True). 
    #
    # \param saveIteration If the generated OpenMDAO component is going to be called multiple times, 
    # the inputs and outputs from the analysis and the AIM will be
    # automatically bookkept ( = True) by moving the files to a folder within the AIM's analysis directory 
    # ( \ref analysisDir ) named "Iteration_#" were # represents the iteration number (default - False). By default ( = False)
    # input and output files will be continously overwritten. 
    # Notes:
    # - If the AIM has 'parents' their genertated files will not be bookkept. 
    # - If previous iteration folders already exist, the iteration folders and any other files in the directory will be moved to a folder named
    # "Instance_#".
    # - This bookkeeping method will likely fail if the iterations are run concurrently!
    #
    # \param executeCommand Command to be executed when running an external code. Command must be a list of command 
    # line arguements (see OpenMDAO documentation). If provided an ExternalCode[1.7.2]/ExternalCodeComp[2.8+]
    # object is created; if not provided or set to None a Component[1.7.3]/ExplicitComponent[2.8+] object is created (default - None).
    #
    # \param inputFile Optional list of input file names for OpenMDAO to check the existence of before 
    # OpenMDAO excutes the "solve_nonlinear"[1.7.3]/"compute"[2.8+] (default - None). This is redundant as the AIM automatically 
    # does this already.
    #
    # \param outputFile Optional list of output names for OpenMDAO to check the existence of before 
    # OpenMDAO excutes the "solve_nonlinear"[1.7.3]/"compute"[2.8+]  (default - None). This is redundant as the AIM automatically 
    # does this already.
    # 
    # \param stdin Set I/O connection for the standard input of an ExternalCode[1.7.2]/ExternalCodeComp[2.8+] component. 
    # The use of this depends on the expected AIM execution. 
    # 
    # \param stdout Set I/O connection for the standard ouput of an ExternalCode[1.7.2]/ExternalCodeComp[2.8+] component. 
    # The use of this depends on the expected AIM execution. 
    # 
    # \param setSensitivity Optional dictionary containing sensitivity/derivative settings/parameters. Currently only 
    # Finite difference is is supported!. See OpenMDAO 
    # documentation for addtional details of "deriv_options"(version 1.7) or "declare_partials"(version 2.8). 
    # Common values for a finite difference calculation would 
    # be setSensitivity['type'] = "fd" (Note in the version 2.8 documentation this varibale has been changed to "method" 
    # both variations will work when using version 2.8+), setSensitivity['form'] =  "forward" or "backward" or "central", and  
    # setSensitivity['step_size'] = 1.0E-6 (Note in the version 2.8 documentation this varibale has been changed to "step" 
    # both variations will work when using version 2.8+).
    #
    # \return Optionally returns the reference to the OpenMDAO component object created; a reference to the component is stored 
    # in \ref openMDAOComponent . "None" is returned if a failure occurred during object creation.
    def createOpenMDAOComponent(self, inputVariable, outputVariable, **kwargs):

        print("Creating an OpenMDAO component for analysis, " + str(self.aimName))
        
        self.openMDAOComponent = createOpenMDAOComponent(self, 
                                                         inputVariable, 
                                                         outputVariable, 
                                                         kwargs.pop("executeCommand", None), 
                                                         kwargs.pop("inputFile", None), 
                                                         kwargs.pop("outputFile", None), 
                                                         kwargs.pop("setSensitivity", {}), 
                                                         kwargs.pop("changeDir", True), 
                                                         kwargs.pop("saveIteration", False))
        if isinstance(self.openMDAOComponent, int):
            self.capsProblem.status = self.openMDAOComponent
            self.openMDAOComponent = None
        
            self.capsProblem.checkStatus(msg = "while getting creating OpenMDAO component for analysis - " 
                                         + str(self.aimName))
            return None
        else:
            
            if kwargs.get("stdin", None) is not None:
                self.openMDAOComponent.stdin = kwargs.get("stdin", None)
                
            if kwargs.get("stdout", None) is not None:
                self.openMDAOComponent.stdout = kwargs.get("stdout", None)
            
            return self.openMDAOComponent
    
    ## Create a HTML dendrogram/tree of the current state of the analysis. See 
    # example \ref analysis5.py for a representative use case. The HTML file relies on 
    # the open-source JavaScript library, D3, to visualize the data. 
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file. 
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to 
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet 
    # access will be necessary to view the tree.  
    #
    # \param filename Filename to use when saving the tree (default - "aimName"). Note an ".html" is 
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False 
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access 
    # will be necessary to view the tree.  
    #
    # \param analysisGeom Show the geometry currently load into the analysis in the tree (default - False).
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False). Note: "analysisGeom" must also be set to True.
    #
    # \param reverseMap Reverse the attribute map (default - False). See \ref getAttributeMap for details.
    def createTree(self, filename = "aimName", **kwargs):
        ## \example analysis5.py
        # Example use cases for the pyCAPS.capsAnalysis.createTree() function.
      
        if filename == "aimName":
            filename = self.aimName
            
        embedJSON = kwargs.pop("embedJSON", True)
        
        if embedJSON is False:
            try:
                # Write the JSON data file 
                fp = open(filename + ".json",'w')
                fp.write( createTree(self, 
                                     kwargs.pop("analysisGeom", False),
                                     kwargs.pop("internalGeomAttr", False),
                                     kwargs.pop("reverseMap", False)
                                    )
                         )
                fp.close()
                
                # Write the HTML file
                self.capsProblem.status = writeTreeHTML(filename +".html", 
                                                        filename + ".json",
                                                        None,
                                                        kwargs.pop("internetAccess", True)
                                                        )
            except:
                self.capsProblem.status = cCAPS.CAPS_IOERR
        else: 
            self.capsProblem.status = writeTreeHTML(filename+".html", 
                                                    None, 
                                                    createTree(self, 
                                                               kwargs.pop("analysisGeom", False),
                                                               kwargs.pop("internalGeomAttr", False),
                                                               kwargs.pop("reverseMap", False)
                                                               ),
                                                    kwargs.pop("internetAccess", True)
                                                    )
        
        self.capsProblem.checkStatus(msg = "while creating HTML tree")
    
    ## Add an attribute (that is meta-data) to the analysis object. See example 
    # \ref analysis7.py for a representative use case.
    # 
    # \param name Name used to define the attribute.
    # 
    # \param data Data value(s) for the attribute. Note that type casting in done 
    # automatically based on the determined type of the Python object.
    def addAttribute(self, name, data):
        ## \example analysis7.py
        # Example use cases for interacting the 
        # pyCAPS.capsAnalysis.addAttribute() and 
        # pyCAPS.capsAnalysis.getAttribute() functions.
        
        valueObj = _capsValue(self.capsProblem, name, data, subType="User")
            
        self.capsProblem.status = cCAPS.caps_setAttr(self.analysisObj, 
                                                     NULL, 
                                                     valueObj.valueObj)
        self.capsProblem.checkStatus(msg = "while adding attribute for analysis - " 
                                     + str(self.aimName))
        
        # valueObj should be automatically deleted after the function loses scopes
        # del valueObj # not neccessary
    
    ## Get an attribute (that is meta-data) that exists on the analysis object. See example 
    # \ref analysis7.py for a representative use case.
    # \param name Name of attribute to retrieve. 
    # \return Value of attribute "name".
    def getAttribute(self, name): 
        cdef:
            cCAPS.capsObj valueObj
            
        name = _byteify(name)
        self.capsProblem.status = cCAPS.caps_attrByName(self.analysisObj,
                                                        <char *> name, 
                                                        &valueObj)
        self.capsProblem.checkStatus(msg = "while getting attribute for analysis - " 
                                     + str(self.aimName))
        
        return _capsValue(self.capsProblem, None, None).__setupObj(valueObj).value
    
    ## Get sensitivity values. Note the AIM must have aimBackdoor function to interact with. 
    # \param inputVar Input variable to retrieve the sensitivity with respect to. 
    # \param outputVar Output variable to retrieve the sensitivity value for. 
    # \return Sensitivity value. 
    # 
    # Note for AIM developers:
    # 
    # This function makes use of the caps_AIMbackdoor function. The JSONin variable provided is a JSON dictionary with 
    # the following form - {"mode": "Sensitivity", "outputVar": outputVar, "inputVar": inputVar}. Similarly, a JSON dictionary 
    # is expected to return (JSONout) with the following form - {"sensitivity": sensitivityVal}; furthermore it is 
    # assumed that JSON string returned (JSONout) is freeable.
    def getSensitivity(self, inputVar, outputVar):
        
        cdef:
            char *JSONout=NULL # Free able 
                                           
        JSONin = _byteify(jsonDumps({"mode": "Sensitivity", 
                                     "outputVar": str(outputVar), 
                                     "inputVar" : str(inputVar)}))
        
        self.capsProblem.status = cCAPS.caps_AIMbackdoor(self.analysisObj, 
                                                         <char *> JSONin, 
                                                         &JSONout)
        
        self.capsProblem.checkStatus(msg = "while getting sensitivity value for analysis - " + str(self.aimName))
        
        if JSONout:
            temp =  _strify(jsonLoads( _strify(<object> JSONout), object_hook=_byteify))
             
            cEGADS.EG_free(JSONout)
            return float(temp["sensitivity"])
        
        return None
    
    ## Alternative to \ref backDoor(). Warning: May be deprecated in future versions.
    def aimBackDoor(self, JSONin):
        print("\n\nWarning: aimBackDoor is a deprecated. Please use backDoor.\n",
                   "Further use will result in an error in future releases!!!!!!\n")
        self.backDoor(JSONin)
        
    ## Make a call to the AIM backdoor function. Important: it is assumed that the JSON string returned (JSONout) by the AIM
    # is freeable.
    #
    # \param JSONin JSON string input to the backdoor function. If the value isn't a string (e.g. a Python dictionary)
    # it will automatically be converted to JSON string. 
    #
    # \return Returns a Python object of the converted JSON out string from the back door function if any.
    def backDoor(self, JSONin):
        cdef:
           
            char *JSONout = NULL
            object output = None
            
        if not isinstance(JSONin, str):
            JSONin = jsonDumps(JSONin)
        
        JSONin = _byteify(JSONin)
        
        self.capsProblem.status = cCAPS.caps_AIMbackdoor(self.analysisObj, 
                                                         <char *> JSONin, 
                                                         &JSONout)
        
        self.capsProblem.checkStatus(msg = "while calling AIMbackdoor for analysis - " + str(self.aimName))
        
        if JSONout:
            output = _strify(jsonLoads( _strify(<object> JSONout), object_hook=_byteify))
            
            cEGADS.EG_free(JSONout)
            return output
        
        return None

    ## View the geometry associated with the analysis. 
    #  If the analysis produces a surface tessellation, then that is shown. 
    #  Otherwise the bodies are shown with default tessellation parameters.
    #  Note that the geometry must be built and will not autoamtically be built by this function.
    #  
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param portNumber Port number to start the server listening on (default - 7681).
    def viewGeometry(self, **kwargs):
        cdef: 
            int numBody = 0
            cEGADS.ego *bodies = NULL # Array of bodies to display
            
        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))
        
        # Force regeneration of the geometry
        self.capsProblem.geometry.buildGeometry()
        
        # Get the bodies with possible tessellations from the analysis
        self.capsProblem.status = cCAPS.caps_getBodies(self.analysisObj, &numBody, &bodies)
        self.capsProblem.checkStatus(msg = "while generating view of geometry"
                                     + " (during a call to caps_getBodies)")

        if numBody == 0 or bodies == <cEGADS.ego*> NULL:
            print("No bodies available for display.")
            return
        
        #viewTess = False
        #for i in range(numBody):
        #    if bodies[i+numBody] != <cEGADS.ego> NULL:
        #        viewTess = True
        #        break
        
        #if viewTess:
        #    print("Viewing analysis tessellation...")
        #    # Show the tessellation used by the aim
        #    viewer.addEgadsTess(numBody, bodies+numBody)
        #else:
        # Show the bodies
        print("Viewing analysis bodies...")
        viewer.addEgadsBody(numBody, bodies)

        # Start up server 
        viewer.startServer()
        
        del viewer

    ## Save the current geometry used by the AIM to a file.
    # \param filename File name to use when saving geometry file.   
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    def saveGeometry(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):
        cdef: 
            int numBody = 0
            cEGADS.ego *bodies = NULL # Array of bodies to display

        filename = self.capsProblem._createFilename(filename, directory, extension)
        
        # Force regeneration of the geometry
        self.capsProblem.geometry.buildGeometry()

         # Get the bodies with possible tessellations from the analysis
        self.capsProblem.status = cCAPS.caps_getBodies(self.analysisObj, &numBody, &bodies)
        self.capsProblem.checkStatus(msg = "while saving analysis geometry"
                                     + " (during a call to caps_getBodies)")

        if numBody == 0 or bodies == <cEGADS.ego*> NULL:

            self.capsProblem.status = cCAPS.CAPS_SOURCEERR
            self.capsProblem.checkStatus(msg = "No bodies available to save!")

        self.capsProblem.status = saveBodies(numBody, bodies, filename)
        self.capsProblem.checkStatus(msg = "while saving geometry.")
    
    ## Get bounding boxes for the geometry used by the AIM.
    # 
    # \return Returns a dictionary of the bounding boxes, {"BodyName", [ x_min, y_min, z_min, x_max, y_max, z_max]}.
    def getBoundingBox(self):
        cdef: 
            int numBody = 0
            cEGADS.ego *bodies = NULL # Array of bodies
            
            double box[6]
            
        boundingBox = {}
        nameList = []
        
        # Force regeneration of the geometry
        self.capsProblem.geometry.buildGeometry()
        
        self.capsProblem.status = cCAPS.caps_getBodies(self.analysisObj, &numBody, &bodies)
        self.capsProblem.checkStatus(msg = "while geting bounding box"
                                     + " (during a call to caps_getBodies)")
        
        if numBody == 0 or bodies == <cEGADS.ego*> NULL:
            self.capsProblem.status = cCAPS.CAPS_SOURCEERR
            self.capsProblem.checkStatus(msg = "No bodies available to retrieve the bounding box for!")

        for i in range(numBody):
            self.capsProblem.status = cEGADS.EG_getBoundingBox(bodies[i], box)
            self.capsProblem.checkStatus(msg = "while geting bounding box"
                                         + " (during a call to EG_getBoundingBox)")
        
            
            name, nameList = self.capsProblem.__getBodyName(i, nameList, bodies[i])
            
            boundingBox[name] = [box[0], box[1], box[2], 
                                 box[3], box[4], box[5]]
        
        return boundingBox
    
## Defines a CAPS geometry object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.loadCAPS .
cdef class capsGeometry:
    cdef:
        capsProblem capsProblem
        
        readonly object geomName
    
    ## Initializes a geometry object. See \ref capsProblem.loadCAPS for additional information and caveats. 
    #
    # \param problemObject capsProblem object in which the geometry object will be added to. 
    #
    # \param capsFile CAPS file to load. Options: *.csm, *.caps, or *.egads. If the filename has a *.caps extension 
    # the pyCAPS analysis, bound, and value objects will be re-populated (see remarks in \ref capsProblem.loadCAPS ).
    #  
    # \param projectName CAPS project name. projectName=capsFile if not provided.
    def __cinit__(self, problemObject, capsFile, projectName=None):
        
        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError
        
        self.capsProblem = problemObject
    
        ## Geometry file loaded into problem. Note that the directory path has been removed.  
        self.geomName = os.path.basename(capsFile)
        
        if projectName is None: 
            projectName = capsFile 
        
        if self.capsProblem.geometry is not None:
            print("\nWarning: Can NOT load multiple files into a problem. A CAPS file has already been loaded!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while loading file - " + str(capsFile) + ". Can NOT load multiple files into a problem.")
        
        # Check to see if file is coming from a url 
        #if capsFile.startswith("http"):
        #    import urllib2
        #    capsFile.split
        #    with open('test.mp3','wb') as f:
        #       f.write(urllib2.urlopen(capsFile).read())
         
        capsFile = _byteify(capsFile)
        projectName = _byteify(projectName)
        self.capsProblem.status = cCAPS.caps_open(<const char *> capsFile, 
                                                  <const char *> projectName, 
                                                  &self.capsProblem.problemObj)
        capsFile = _strify(capsFile)
        self.capsProblem.checkStatus(msg = "while loading file - " + str(capsFile))
        
        self.capsProblem.geometry = self
        
        if ".caps" in capsFile: 
            self.capsProblem.__repopulateProblem()
            
    ## Manually force a build of the geometry. The geometry will only be rebuilt if a
    # design parameter (see \ref setGeometryVal) has been changed.
    def buildGeometry(self):
        cdef:
            cCAPS.capsErrs *errors
            int nErr

        # Check the problem object    
        if self.capsProblem.problemObj.blind == NULL:
            self.capsProblem.status = cCAPS.CAPS_NULLBLIND
            self.capsProblem.checkStatus(msg = "while building geometry.")
        
        # Force regeneration of the geometry
        self.capsProblem.status = cCAPS.caps_preAnalysis(self.capsProblem.problemObj, &nErr, &errors)
        
        # Check errors
        report_Errors(nErr, errors)
        
        if (self.capsProblem.status != cCAPS.CAPS_SUCCESS and 
            self.capsProblem.status != cCAPS.CAPS_CLEAN):
            
            self.capsProblem.checkStatus(msg = "while building geometry.")

    ## Save the current geometry to a file.
    # \param filename File name to use when saving geometry file.   
    # \param directory Directory where to save file. Default current working directory.
    # \param extension Extension type for file if filename does not contain an extension.
    def saveGeometry(self, filename="myGeometry", directory=os.getcwd(), extension=".egads"):
        
        filename = self.capsProblem._createFilename(filename, directory, extension)

        # Force regeneration of the geometry
        self.buildGeometry()

        self.capsProblem.status = saveGeometry(self.capsProblem.problemObj,
                                               filename)
        self.capsProblem.checkStatus(msg = "while saving geometry.")
    
    ## Sets a GEOMETRYIN variable for the geometry object.
    # \param varname Name of geometry design parameter to set in the geometry.
    # \param value Value of geometry design parameter. Type casting is automatically done based 
    # on value indicated by CAPS.
    def setGeometryVal(self, varname=None, value=None):
    
        cdef: 
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int vlen
            int nErr
            const char *valunits
        
            void *data
             
        if varname is None or value is None:
            print("No variable name and/or value provided")
            return

        #const char *valunits
        varname = _byteify(varname)
        self.capsProblem.status = cCAPS.caps_childByName(self.capsProblem.problemObj,
                                                         cCAPS.VALUE,
                                                         cCAPS.GEOMETRYIN,
                                                         <char *> varname,
                                                         &tempObj)
        self.capsProblem.checkStatus(msg = "while setting geometry variable - " + str(varname) 
                                         + " (during a call to caps_childByName)")
    
        self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                      &vtype,
                                                      &vlen,
                                                      <const void ** > &data,
                                                      &valunits,
                                                      &nErr,
                                                      &errors)
        # Check errors
        report_Errors(nErr, errors)
        
        self.capsProblem.checkStatus(msg = "while setting geometry variable - " + str(varname) 
                                         + " (during a call to caps_getValue)")
        
        numRow, numCol = getPythonObjShape(value)
            
        byteValue  = []
        data = castValue2VoidP(value, byteValue, vtype)
        
        try:
            self.capsProblem.status = cCAPS.caps_setValue(tempObj, 
                                                          <int> numRow, 
                                                          <int> numCol, 
                                                          <const void *> data)
        finally:
            if isinstance(byteValue,list):
                del byteValue[:]
            
            if data and vtype != cCAPS.String:
                free(data)
            
            
            
        self.capsProblem.checkStatus(msg = "while setting geometry variable - " + str(varname) 
                                         + " (during a call to caps_setValue)")

    ## Gets a GEOMETRYIN variable for the geometry object.
    # \param varname Name of geometry design parameter to retrieve from the geometry. 
    # If no name is provided a dictionary containing all design parameters is returned.
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param namesOnly Return only a list of parameter names (no values) if creating a 
    # dictionary (default - False).
    #
    # \return Value of "varname" or dictionary of all values. 
    def getGeometryVal(self, varname=None, **kwargs):

        cdef: 
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int vlen
            int nErr
            const char *valunits
            const void *data

            int vdim, vnrow, vncol
            cCAPS.capsFixed vlfixed, vsfixed
            cCAPS.capsNull  vntype
            
            int index, numGeometryInput # returning a dictionary variables 
            char * name
        
        geometryIn = {}
        
        namesOnly = kwargs.pop("namesOnly", False)
        
        if varname is not None: # Get just the 1 variable
        
            varname = _byteify(varname)
            self.capsProblem.status = cCAPS.caps_childByName(self.capsProblem.problemObj,
                                                             cCAPS.VALUE,
                                                             cCAPS.GEOMETRYIN,
                                                             <char *> varname,
                                                             &tempObj)
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_childByName)")
            
            self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                          &vtype,
                                                          &vlen,
                                                          &data,
                                                          &valunits,
                                                          &nErr,
                                                          &errors)
            
            # Check errors
            report_Errors(nErr, errors)
            
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_getValue)")
            
            self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                               &vdim,
                                                               &vlfixed,
                                                               &vsfixed,
                                                               &vntype,
                                                               &vnrow,
                                                               &vncol)
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_getValueShape)")

            if vntype == cCAPS.IsNull:
                return None
            else:
                return castValue2PythonObj(data, vtype, vnrow, vncol)
        
        else: # Build a dictionary 
             
            self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj,
                                                      cCAPS.VALUE,
                                                      cCAPS.GEOMETRYIN,
                                                      &numGeometryInput)
            self.capsProblem.checkStatus(msg = "while getting number of geometry design paramaters" 
                                             + " (during a call to caps_size)")
            
            for i in range(1, numGeometryInput+1):
                self.capsProblem.status = cCAPS.caps_childByIndex(self.capsProblem.problemObj,
                                                                  cCAPS.VALUE,
                                                                  cCAPS.GEOMETRYIN,
                                                                  i,
                                                                  &tempObj)
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_childByIndex)")
                 
                name = tempObj.name
                 
                self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                              &vtype,
                                                              &vlen,
                                                              &data,
                                                              &valunits,
                                                              &nErr,
                                                              &errors)
                # Check errors
                report_Errors(nErr, errors)
             
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValue)")
             
                self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                                   &vdim,
                                                                   &vlfixed,
                                                                   &vsfixed,
                                                                   &vntype,
                                                                   &vnrow,
                                                                   &vncol)
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValueShape)")
     
             
                if vntype == cCAPS.IsNull:
                    geometryIn[<object> _strify(name)] = None
                else:
                    geometryIn[<object> _strify(name)] = castValue2PythonObj(data, vtype, vnrow, vncol)
             
            if namesOnly is True:
                return geometryIn.keys()
            else:
                return geometryIn
            
    ## Gets a GEOMETRYOUT variable for the geometry object.
    # \param varname Name of geometry (local) parameter to retrieve from the geometry. 
    # If no name is provided a dictionary containing all local variables is returned.
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param namesOnly Return only a list of parameter names (no values) if creating a dictionary (default - False).
    # \param ignoreAt Ignore @ geometry variables when creating a dictionary (default - True).
    #
    # \return Value of "varname" or dictionary of all values. 
    def getGeometryOutVal(self, varname=None, **kwargs):

        cdef: 
            cCAPS.capsObj tempObj
            cCAPS.capsoType type
            cCAPS.capsvType vtype
            cCAPS.capsErrs *errors
        
            int vlen
            int nErr
            const char *valunits
            const void *data

            int vdim, vnrow, vncol
            cCAPS.capsFixed vlfixed, vsfixed
            cCAPS.capsNull  vntype
            
            int index, numGeometryOutput # returning a dictionary variables 
            char * name
        
        # Make sure the geometry is built
        self.buildGeometry()
        
        geometryOut = {}
        
        namesOnly = kwargs.pop("namesOnly", False)
        ignoreAt  = kwargs.pop("ignoreAt", True)
        
        if varname is not None: # Get just the 1 variable
        
            varname = _byteify(varname)
            self.capsProblem.status = cCAPS.caps_childByName(self.capsProblem.problemObj,
                                                             cCAPS.VALUE,
                                                             cCAPS.GEOMETRYOUT,
                                                             <char *> varname,
                                                             &tempObj)
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_childByName)")
            
            self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                          &vtype,
                                                          &vlen,
                                                          &data,
                                                          &valunits,
                                                          &nErr,
                                                          &errors)
            
            # Check errors
            report_Errors(nErr, errors)
            
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_getValue)")
            
            self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                               &vdim,
                                                               &vlfixed,
                                                               &vsfixed,
                                                               &vntype,
                                                               &vnrow,
                                                               &vncol)
            self.capsProblem.checkStatus(msg = "while getting geometry variable - " + str(varname) 
                                             + " (during a call to caps_getValueShape)")

            if vntype == cCAPS.IsNull:
                return None
            else:
                return castValue2PythonObj(data, vtype, vnrow, vncol)
        
        else: # Build a dictionary 
             
            self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj,
                                                      cCAPS.VALUE,
                                                      cCAPS.GEOMETRYOUT,
                                                      &numGeometryOutput)
            self.capsProblem.checkStatus(msg = "while getting number of geometry local variables" 
                                             + " (during a call to caps_size)")
            
            for i in range(1, numGeometryOutput+1):
                self.capsProblem.status = cCAPS.caps_childByIndex(self.capsProblem.problemObj,
                                                                  cCAPS.VALUE,
                                                                  cCAPS.GEOMETRYOUT,
                                                                  i,
                                                                  &tempObj)
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_childByIndex)")
                 
                name = tempObj.name
                varname = _strify(name)
                if ignoreAt: 
                    if "@" in varname:
                        continue
                     
                self.capsProblem.status = cCAPS.caps_getValue(tempObj,
                                                              &vtype,
                                                              &vlen,
                                                              &data,
                                                              &valunits,
                                                              &nErr,
                                                              &errors)
                # Check errors
                report_Errors(nErr, errors)
             
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValue)")
             
                self.capsProblem.status = cCAPS.caps_getValueShape(tempObj,
                                                                   &vdim,
                                                                   &vlfixed,
                                                                   &vsfixed,
                                                                   &vntype,
                                                                   &vnrow,
                                                                   &vncol)
                self.capsProblem.checkStatus(msg = "while getting geometry variable (index) - " + str(i+1) 
                                                 + " (during a call to caps_getValueShape)")
     
             
                if vntype == cCAPS.IsNull:
                    geometryOut[varname] = None
                else:
                    geometryOut[varname] = castValue2PythonObj(data, vtype, vnrow, vncol)
             
            if namesOnly is True:
                return geometryOut.keys()
            else:
                return geometryOut
    
    ## View or take a screen shot of the geometry configuration. The use of this function to save geometry requires the \b matplotlib module. 
    # \b <em>Important</em>: If both showImage = True and filename is not None, any manual view changes made by the user
    # in the displayed image will be reflected in the saved image.
    #  
    # \param **kwargs See below.
    # 
    # Valid keywords:
    # \param viewerType What viewer should be used (default - "capsViewer"). 
    #                   Options: "capsViewer" or "matplotlib" (options are case insensitive). 
    #                   Important: if $filename is not None, the viewer is changed to matplotlib.
    # \param portNumber Port number to start the server listening on (default - 7681).
    # \param title Title to add to each figure (default - None). 
    # \param filename Save image(s) to file specified (default - None). Note filename should 
    # not contain '.' other than to indicate file type extension (default type = *.png). 
    # 'file' - OK, 'file2.0Test' - BAD, 'file2_0Test.png' - OK, 'file2.0Test.jpg' - BAD. 
    # \param directory Directory path were to save file. If the directory doesn't exist 
    # it will be made. (default - current directory). 
    # \param viewType  Type of view for the image(s). Options: "isometric" (default),
    # "fourview", "top" (or "-zaxis"), "bottom" (or "+zaxis"), "right" (or "+yaxis"), 
    # "left" (or "-yaxis"), "front" (or "+xaxis"), "back" (or "-xaxis").
    # \param combineBodies Combine all bodies into a single image (default - False).
    # \param ignoreBndBox Ignore the largest body (default - False).
    # \param showImage Show image(s) (default - False).
    # \param showAxes Show the xyz axes in the image(s) (default - False).
    # \param showTess Show the edges of the tessellation (default - False). 
    # \param dpi Resolution in dots-per-inch for the figure (default - None).
    # \param tessParam Custom tessellation paremeters, see EGADS documentation for makeTessBody function. 
    # values will be scaled by the norm of the bounding box for the body (default - [0.0250, 0.0010, 15.0]).
    def viewGeometry(self, **kwargs):
        
        viewDefault = "capsviewer"

        viewer = kwargs.pop("viewerType", viewDefault).lower()
        filename = kwargs.get("filename", None)
        
        # Are we going to use the viewer?
        if filename is None and viewer.lower() == "capsviewer":
            self._viewGeometryCAPSViewer(**kwargs)
        else:
            self._viewGeometryMatplotlib(**kwargs)

    def _viewGeometryCAPSViewer(self, **kwargs):
        
        cdef: 
            
            int bodyIndex = 0 # Indexing
            int numBody = 0
            cEGADS.ego *bodies  = NULL # Array of bodies to display
            char *units = NULL

        def cleanup():
            cEGADS.EG_free(bodies)
        
        # Force regeneration of the geometry   
        self.buildGeometry()

        # Initialize the viewer class
        viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))
        
         # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while generating view of geometry"
                                     + " (during a call to caps_size)")

        try:
            # Allocate the array
            bodies = <cEGADS.ego *> cEGADS.EG_alloc(numBody*sizeof(cEGADS.ego))

            # Loop through and get all bodies
            for bodyIndex in range(numBody):
                self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, bodyIndex+1, &bodies[bodyIndex], &units)
                self.capsProblem.checkStatus(msg = "while generating view of geometry"
                                             + " (during a call to caps_bodyByIndex)")

            # Load the bodies into the viewer
            viewer.addEgadsBody(numBody, bodies)
        except:
            cleanup()
            raise

        # Start up server 
        viewer.startServer()
        
        del viewer
        cleanup()

    def _viewGeometryMatplotlib(self, **kwargs):

        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
        except:
            print("Error: Unable to import matplotlib - viewing the geometry with matplotlib is not possible")
            raise ImportError 
        
        def createFigure(figure, view, xArray, yArray, zArray, triArray, linewidth, edgecolor, transparent, facecolor):
            figure.append(plt.figure())  
            
            if view == "fourview": 
                for axIndex in range(4):
                    figure[-1].add_subplot(2, 2, axIndex+1, projection='3d')
                    
            else:
                figure[-1].gca(projection='3d')
                
            for ax in figure[-1].axes:
                ax.plot_trisurf(xArray, 
                                yArray, 
                                zArray, 
                                triangles=triArray,
                                linewidth=linewidth,
                                edgecolor=edgecolor,
                                alpha=transparent,
                                color=facecolor)
               
        cdef: 
            
            int i = 0, bodyIndex = 0, faceIndex = 0, bndBoxIndex = 0# Indexing
            
            int numBody = 0, numFace = 0
            int numPointTot = 0
            
            cEGADS.ego body # Temporary body container
            cEGADS.ego tess # Body tessellation 
            cEGADS.ego *faces = NULL # Faces on the body 
            
            double box[6]
            double size = 0, sizeBndBox = 0
            double params[3]
                     
            char *units = NULL
            
            #EGADS function return variables
            int plen = 0, tlen = 0
            const int *tris = NULL
            const int *triNeighbor = NULL
            const int *ptype = NULL
            const int *pindex = NULL
            
            const double *points = NULL
            const double *uv = NULL

            int     *triConn = NULL
            double  *xyzs = NULL
            
            int globalID = 0
            
            # Data arrays 
            object xArray = [], yArray = [], zArray = [], triArray = [], gIDArray = []
            
            # Figure 
            object axIndex = 0, figIndex =0 # Indexing 
            object figure = []
            
        # Plot parameters from keyword Args
        facecolor = "#D3D3D3"
        edgecolor = facecolor
        transparent = 1.0
        linewidth = 0
        
        viewType = ["isometric", "fourview", 
                    "top"   , "-zaxis", 
                    "bottom", "+zaxis", 
                    "right",  "+yaxis", 
                    "left",   "-yaxis", 
                    "front",  "+xaxis", 
                    "back",   "-xaxis"]
        
        title = kwargs.pop("title", None)
        combineBodies = kwargs.pop("combineBodies", False)
        ignoreBndBox = kwargs.pop("ignoreBndBox", False)
        showImage =  kwargs.pop("showImage", False)
        showAxes =  kwargs.pop("showAxes", False)
        
        showTess = kwargs.pop("showTess", False)
        dpi = kwargs.pop("dpi", None)
        
        tessParam = kwargs.pop("tessParam", [0.0250, 0.0010, 15.0])
                    
        if showTess:
            edgecolor = "black"
            linewidth = 0.2
        
        filename = kwargs.pop("filename", None)
        directory = kwargs.pop("directory", None)
        
        if (filename is not None and "." not in filename):
            filename += ".png"
        
        view = str(kwargs.pop("viewType", "isometric")).lower()
        if view not in viewType:
            view = viewType[0]
            print("Unrecongized viewType, defaulting to " + str(view))
                
        # Force regeneration of the geometry   
        self.buildGeometry()  
            
        # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                         + " (during a call to caps_size)")
    
        if numBody < 1:
            print("The number of bodies in the problem is less than 1!!!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while saving screen shot of geometry")
        
        # If ignoring the bounding box body, determine which one it is 
        if numBody > 1 and ignoreBndBox == True:
            for bodyIndex in range(numBody):
                self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, bodyIndex+1, &body, &units)
                self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                             + " (during a call to caps_bodyByIndex)")
                
                self.capsProblem.status = cEGADS.EG_getBoundingBox(body, box)
                self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                             + " (during a call to EG_getBoundingBox)")
    
                size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                            (box[1]-box[4])*(box[1]-box[4]) +
                            (box[2]-box[5])*(box[2]-box[5]))
                if size > sizeBndBox:
                    sizeBndBox = size
                    bndBoxIndex = bodyIndex
        else: 
            ignoreBndBox = False
                    
        # Loop through bodies 1 by 1 
        for bodyIndex in range(numBody):
            
            # If ignoring the bounding box - skip this body 
            if ignoreBndBox and bodyIndex == bndBoxIndex:
                continue
                
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, bodyIndex+1, &body, &units)
            self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                         + " (during a call to caps_bodyByIndex)")

            # Get tessellation             
            self.capsProblem.status = cEGADS.EG_getBoundingBox(body, box)
            self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                         + " (during a call to EG_getBoundingBox)")

            size = sqrt((box[0]-box[3])*(box[0]-box[3]) +
                        (box[1]-box[4])*(box[1]-box[4]) +
                        (box[2]-box[5])*(box[2]-box[5]))


            params[0] = tessParam[0]*size
            params[1] = tessParam[1]*size
            params[2] = tessParam[2]

            self.capsProblem.status =  cEGADS.EG_makeTessBody(body, params, &tess)
            self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                         + " (during a call to EG_makeTessBody)")
            
            # How many faces do we have? 
            self.capsProblem.status = cEGADS.EG_getBodyTopos(body, 
                                                             <cEGADS.ego> NULL, 
                                                             cEGADS.FACE, 
                                                             &numFace, 
                                                             &faces)
            self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                         + " (during a call to EG_makeTessBody)")
       
            # If not combining bodies erase arrays 
            if not combineBodies:
                del xArray[:]
                del yArray[:]
                del zArray[:]
                del triArray[:]
                
            del gIDArray[:]
            
            # Loop through face by face and get tessellation information 
            for faceIndex in range(numFace):
                self.capsProblem.status = cEGADS.EG_getTessFace(tess, faceIndex+1, &plen, &points, &uv, &ptype, 
                                                                &pindex, &tlen, &tris, &triNeighbor)
                self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                             + " (during a call to EG_getTessFace)")

                for i in range(plen):
                     self.capsProblem.status =  cEGADS.EG_localToGlobal(tess, faceIndex+1, i+1, &globalID)
                     self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                                  + " (during a call to EG_localToGlobal)")
                     
                     tri = []
                     if globalID not in gIDArray:
                         xArray.append(points[3*i  ])
                         yArray.append(points[3*i+1])
                         zArray.append(points[3*i+2])
                         gIDArray.append(globalID)
                
                for i in range(tlen):
                    tri = []
                    self.capsProblem.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*i + 0], &globalID)
                    self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                                  + " (during a call to EG_localToGlobal)")
                    
                    tri.append(globalID-1 + numPointTot)
                    
                    self.capsProblem.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*i + 1], &globalID)
                    self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                                  + " (during a call to EG_localToGlobal)")
                    
                    tri.append(globalID-1 + numPointTot)
                    
                    self.capsProblem.status = cEGADS.EG_localToGlobal(tess, faceIndex+1, tris[3*i + 2], &globalID)
                    self.capsProblem.checkStatus(msg = "while saving screen shot of geometry"
                                                  + " (during a call to EG_localToGlobal)")
                    
                    tri.append(globalID-1 + numPointTot)
                    
                    triArray.append(tri)
            
            if combineBodies:
                numPointTot += len(gIDArray)
            else: 
                createFigure(figure, view, 
                             xArray, yArray, zArray, triArray, 
                             linewidth, edgecolor, transparent, facecolor)

        if combineBodies:
            createFigure(figure, view, 
                         xArray, yArray, zArray, triArray, 
                         linewidth, edgecolor, transparent, facecolor)
        
        for fig in figure:
            
            if len(fig.axes) > 1:
                fig.subplots_adjust(left=2.5, bottom=0.0, right=97.5, top=0.95,
                                    wspace=0.0, hspace=0.0)
            if title is not None:
                fig.suptitle(title)
            
            if dpi is not None:
                fig.set_dpi(dpi)
                
            for axIndex, ax in enumerate(fig.axes):
                xLimits = ax.get_xlim()
                yLimits = ax.get_ylim()
                zLimits = ax.get_zlim()
            
                center = []
                center.append((xLimits[1] + xLimits[0])/2.0)
                center.append((yLimits[1] + yLimits[0])/2.0)
                center.append((zLimits[1] + zLimits[0])/2.0)
            
                scale = max(abs(xLimits[1] - xLimits[0]),
                            abs(yLimits[1] - yLimits[0]),
                            abs(zLimits[1] - zLimits[0]))      
            
                scale = 0.6*(scale/2.0)
                
                ax.set_xlim(-scale + center[0], scale + center[0])
                ax.set_ylim(-scale + center[1], scale + center[1])
                ax.set_zlim(-scale + center[2], scale + center[2])
               
                if showAxes:
                    ax.set_xlabel("x", size="medium")
                    ax.set_ylabel("y", size="medium")
                    ax.set_zlabel("z", size="medium")
                else:
                    ax.set_axis_off() # Turn the axis off
                
                #ax.invert_xaxis() # Invert x- axis
                if axIndex == 0:
                    
                    if view in ["top", "-zaxis"]:
                        ax.view_init(azim=-90, elev=90) # Top
                        
                    elif view in ["bottom", "+zaxis"]:
                        ax.view_init(azim=-90, elev=-90) # Bottom
                    
                    elif view in ["right", "+yaxis"]:
                        ax.view_init(azim=-90, elev=0) # Right
                    
                    elif view in ["left", "-yaxis"]:
                        ax.view_init(azim=90, elev=0) # Left
                    
                    elif view in ["front", "+xaxis"]:
                        ax.view_init(azim=180, elev=0) # Front
                    
                    elif view in ["back", "-xaxis"]:
                        ax.view_init(azim=0, elev=0) # Back
                    
                    else:
                        ax.view_init(azim=120, elev=30) # Isometric
                        
                if axIndex == 1:
                    ax.view_init(azim=180, elev=0) # Front 
                
                elif axIndex == 2:
                    ax.view_init(azim=-90, elev=90) # Top
                
                elif axIndex == 3:
                    ax.view_init(azim=-90, elev=0) # Right
                
        if showImage:
            plt.show() 
                
        if filename is not None:
            
            for figIndex, fig in enumerate(figure):
                if len(figure) > 1:
                    fname = filename.replace(".", "_" + str(figIndex) + ".")
                else:
                    fname = filename
                    
                if directory is not None:
                    if not os.path.isdir(directory):
                        print("Directory ( " + directory + " ) does not currently exist - it will be made automatically")
                        os.makedirs(directory)
                        
                    fname = os.path.join(directory, fname)
                
                print("Saving figure - ", fname)
                fig.savefig(fname, dpi=dpi)
    
    ## Retrieve a list of attribute values of a given name ("attributeName") for the bodies in the current geometry. 
    # Level in which to search the bodies is 
    # determined by the attrLevel keyword argument.
    #
    # \param attributeName Name of attribute to retrieve values for. 
    # \param **kwargs See below.
    # 
    # \return A list of attribute values.
    #
    # Valid keywords:
    # \param bodyIndex Specific body in which to retrieve attribute information from.
    # \param attrLevel Level to which to search the body(ies). Options: \n
    #                  0 (or "Body") - search just body attributes\n
    #                  1 (or "Face") - search the body and all the faces [default]\n
    #                  2 (or "Edge") - search the body, faces, and all the edges\n
    #                  3 (or "Node") - search the body, faces, edges, and all the nodes
    # 
    def getAttributeVal(self, attributeName, **kwargs):
        
        cdef: 
            int i = 0 # Indexing
            
            int numBody = 0
            
            cEGADS.ego body # Temporary body container
            
            char *units = NULL
            
        # Initiate attribute value list     
        attributeList = []
        
        attributeName = _byteify(attributeName)
        
        bodyIndex = kwargs.pop("bodyIndex", None)
        attrLevel = kwargs.pop("attrLevel", "face")
        
        attribute = {"body" : 0, 
                     "face" : 1, 
                     "edge" : 2, 
                     "node" : 3}
        
        if not isinstance(attrLevel, (int,float)):
            if attrLevel.lower() in attribute.keys():
                attrLevel = attribute[attrLevel.lower()]
            
        if int(attrLevel) not in attribute.values():
            print("Error: Invalid attribute level! Defaulting to 'Face'")
            attrLevel = attribute["face"]
                
        # Force regeneration of the geometry
        self.buildGeometry()

        # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) 
                                         + " values on geometry  (during a call to caps_size)")
        
        if numBody < 1:
            print("The number of bodies in the problem is less than 1!!!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values on geometry")
            
        # Loop through bodies
        for i in range(numBody):
            #print"Looking through body", i
            
            if bodyIndex is not None:
                if int(bodyIndex) != i+1:
                    continue
            
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, i+1, &body, &units)
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values (body " 
                                         + str(i+1) + " )"
                                         + " (during a call to caps_bodyByIndex)")
            
            self.capsProblem.status, returnList = createAttributeValList(body,
                                                                         <int> attrLevel, 
                                                                         <char *> attributeName)
            self.capsProblem.checkStatus(msg = "while getting attribute " + str(attributeName) + " values (body " 
                                         + str(i+1) + " )"
                                         + " (during a call to createAttributeValList)")
            
            attributeList = attributeList + returnList
        
#         if isinstance(attributeList, list):
#             return attributeList[0] if len(attributeList) == 1 else attributeList
#         else:
#             return list(set(attributeList)) 
    
        try:
            return list(set(attributeList)) 
        except TypeError:
            print("Unhashable type in attributelist")
            return attributeList
        
    ## Create attribution map (embeded dictionaries) of each body in the current geometry.
    # Dictionary layout: 
    # - Body 1 
    #   - Body : Body level attributes
    #   - Faces 
    #     - 1 : Attributes on the first face of the body
    #     - 2 : Attributes on the second face of the body
    #     - " : ... 
    #   - Edges 
    #     - 1 : Attributes on the first edge of the body
    #     - 2 : Attributes on the second edge of the body
    #     - " : ...
    #   - Nodes : 
    #     - 1 : Attributes on the first node of the body
    #     - 2 : Attributes on the second node of the body
    #     - " : ...
    # - Body 2
    #   - Body : Body level attributes 
    #   - Faces 
    #     - 1 : Attributes on the first face of the body
    #     - " : ...
    #   - ...
    # - ...
    #
    # Dictionary layout (reverseMap = True): 
    # - Body 1 
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not 
    #       - Faces : Face numbers at which the attribute exist
    #       - Edges : Edge numbers at which the attribute exist
    #       - Nodes : Node numbers at which the attribute exist
    #     -  Value : Next value of attribute with the same name
    #       - Body : True if value exist at body level, None if not
    #       -  "   : ...
    #     - ...
    #   - Atribute : Attribute name
    #     - Value : Value of attribute
    #       -   "  : ... 
    #     - ...
    # - Body 2
    #   - Attribute : Attribute name
    #     - Value : Value of attribute
    #       - Body  : True if value exist at body level, None if not 
    #       -   "   : ...
    #     - ...
    #   - ...
    # - ...
    #
    # \param getInternal Get internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param reverseMap Reverse the attribute map (default - False). See above table for details.
    #         
    # \return Dictionary containing attribution map
    def getAttributeMap(self, getInternal = False, **kwargs):
        
        cdef: 
            int bodyIndex = 0# Indexing
            
            int numBody = 0
                        
            cEGADS.ego body # Temporary body container
            
            char *units = NULL
        
        attributeMap = {}
        
        # Force regeneration of the geometry
        self.buildGeometry()

        # How many bodies do we have
        self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while creating attribute map (during a call to caps_size)")
        
        if numBody < 1:
            print("The number of bodies in the problem is less than 1!!!\n")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "while creating attribute map")
            
        for bodyIndex in range(numBody, 0, -1):
            
            # No plus +1 on bodyIndex since we are going backwards
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, bodyIndex, &body, &units)
            self.capsProblem.checkStatus(msg = "while creating attribute map (on body " + str(bodyIndex) + ")"
                                         + " (during a call to caps_bodyByIndex)")    
            
            #{"Body" : None, "Faces" : None, "Edges" : None, "Nodes" : None}

            self.capsProblem.status, attributeMap["Body " + str(bodyIndex)] = createAttributeMap(body, getInternal)
            self.capsProblem.checkStatus(msg = "while creating attribute map (on body " + str(bodyIndex) + ")"
                                         + " (during a call to createAttributeMap)")  
        
        if kwargs.pop("reverseMap", False):
            return reverseAttributeMap(attributeMap)
        else:
            return attributeMap 
    
    ## Create a HTML dendrogram/tree of the current state of the geometry. The HTML file relies on 
    # the open-source JavaScript library, D3, to visualize the data. 
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file. 
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to 
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet 
    # access will be necessary to view the tree.  
    #
    # \param filename Filename to use when saving the tree (default - "myGeometry"). Note an ".html" is 
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False 
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default - True)? If set to True internet access 
    # will be necessary to view the tree.  
    #
    # \param internalGeomAttr Show the internal attributes (denoted by starting with an underscore, for example
    # "_AttrName") that exist on the geometry (default - False). 
    #
    # \param reverseMap Reverse the attribute map (default - False). See \ref getAttributeMap for details.
    def createTree(self, filename = "myGeometry", **kwargs):
            
        embedJSON = kwargs.pop("embedJSON", True)
        
        if embedJSON is False:
            try:
                # Write the JSON data file 
                fp = open(filename + ".json",'w')
                fp.write(createTree(self, 
                                    False, 
                                    kwargs.pop("internalGeomAttr", False), 
                                    kwargs.pop("reverseMap", False)
                                    )
                         )
                fp.close()
                
                # Write the HTML file
                self.capsProblem.status = writeTreeHTML(filename +".html", 
                                                        filename + ".json", 
                                                        None, 
                                                        kwargs.pop("internetAccess", True))
            except:
                self.capsProblem.status = cCAPS.CAPS_IOERR
        else: 
            self.capsProblem.status = writeTreeHTML(filename+".html", 
                                                    None, 
                                                    createTree(self, 
                                                               False, 
                                                               kwargs.pop("internalGeomAttr", False), 
                                                               kwargs.pop("reverseMap", False)
                                                               ),
                                                    kwargs.pop("internetAccess", True)
                                                    )
        
        self.capsProblem.checkStatus(msg = "while creating HTML tree")
        
    ## Get bounding boxes for all geometric bodies.
    # 
    # \return Returns a dictionary of the bounding boxes, {"BodyName", [ x_min, y_min, z_min, x_max, y_max, z_max]}.
    def getBoundingBox(self):
        cdef: 
            int numBody = 0
            cEGADS.ego body
            char *units
            
            double box[6]
            
        boundingBox = {}
        nameList = []
        
        # Force regeneration of the geometry
        self.buildGeometry()
        
        self.capsProblem.status = cCAPS.caps_size(self.capsProblem.problemObj, cCAPS.BODIES, cCAPS.NONE, &numBody)
        self.capsProblem.checkStatus(msg = "while geting bounding box"
                                     + " (during a call to caps_size)")

        if numBody < 1:
            self.capsProblem.status = cCAPS.CAPS_SOURCEERR
            self.capsProblem.checkStatus(msg = "No bodies available to retrieve the bounding box for!")
                

        # Make copies of the bodies to hand to the model
        for i in range(numBody):
            self.capsProblem.status = cCAPS.caps_bodyByIndex(self.capsProblem.problemObj, i+1, &body, &units)
            self.capsProblem.checkStatus(msg = "while geting bounding box"
                                         + " (during a call to caps_bodyByIndex)")

            self.capsProblem.status = cEGADS.EG_getBoundingBox(body, box)
            self.capsProblem.checkStatus(msg = "while geting bounding box"
                                         + " (during a call to EG_getBoundingBox)")
        
            
            name, nameList = self.capsProblem.__getBodyName(i, nameList, body)
            
            boundingBox[name] = [box[0], box[1], box[2], 
                                 box[3], box[4], box[5]]
        
        return boundingBox
                            
## Defines a CAPS bound object.
# May be created directly by providing a problem object or indirectly through \ref capsProblem.createDataBound or 
# \ref capsProblem.createDataTransfer .
cdef class capsBound:

    cdef: 
        capsProblem capsProblem
        cCAPS.capsObj boundObj
        
        readonly object boundName
        readonly object dataSetSrc
        readonly object dataSetDest
        
        readonly object vertexSet
        readonly object variables 
     
#     @property 
#     def dateSet(self):
#         return self.dataSetSrc
    
    ## Initializes a bound object. See \ref capsProblem.createDataBound for additional information and caveats. 
    #
    # \param problemObject capsProblem object in which the analysis object will be added to. 
    #
    # \param **kwargs See below.
   #
    # Valid keywords: 
    # \param capsBound Name of capsBound to use for the data bound. 
    # \param variableName Single or list of variables names to add. 
    # \param aimSrc Single or list of AIM names that will be the data sources for the bound.
    # \param aimDest Single or list of AIM names that will be the data destinations during the transfer. 
    # \param transferMethod Single or list of transfer methods to use during the transfer.
    # \param initValueDest Single or list of initial values for the destication data.
    def __cinit__(self, problemObject, **kwargs):
        
        # Check problem object
        if not isinstance(problemObject, capsProblem) and not isinstance(problemObject, CapsProblem):
            raise TypeError("Provided problemObject is not a capsProblem or CapsProblem class")
        
        ## Reference to the problem object that loaded the AIM during a call to \ref capsProblem.loadAIM .
        self.capsProblem = problemObject
        
        # Look for wild card fieldname and compare variable
        def checkFieldWildCard(variable, fieldName):
            #variable = _byteify(variable)
            for i in fieldName:
                try:
                    i_ind = i.index("*")
                except ValueError:
                    i_ind = -1

                try:
                    if variable[:i_ind] == i[:i_ind]:  
                        return True
                except:
                    pass
                
            return False
        
        def checkAnalysisName(aimList):
            for i in aimList:
                if i not in self.capsProblem.analysis.keys():
                    self.capsProblem.status = cCAPS.CAPS_BADVALUE
                    self.capsProblem.checkStatus(msg="AIM name = " + i + " not found in analysis dictionary")
            
        if "capsTransfer" in kwargs.keys():
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg="The keyword, capsTransfer, should be changed to capsBound")
        
        if "variableRank" in kwargs.keys():
            print("WARNING: Keyword variableRank is not necessary. Please remove it from your script. May raise an error in future releases.")
   
        
        capsBound = kwargs.pop("capsBound", None)
        varName = kwargs.pop("variableName", None)
        aimSrc = kwargs.pop("aimSrc", None)
        aimDest = kwargs.pop("aimDest", [])
        transferMethod = kwargs.pop("transferMethod", [])
        
        initValueDest = kwargs.pop("initValueDest", None)
         
        if not capsBound:
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg="No transfer name (capsBound) provided")
        
        if not varName:
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg="No variable name (variableName) provided")
        
        if not aimSrc:
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg="No source AIM(s) (aimSrc) provided")
            
        if not isinstance(varName, list):
            varName = [varName] # Convert to list
            
        numDataSet = len(varName)

        if not isinstance(aimSrc, list):
            aimSrc = [aimSrc]
            
        temp = []
        for i in aimSrc:
            if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                temp.append(i.aimName)
            elif isinstance(i, str):
                temp.append(i)
            else:    
                raise TypeError("Provided aimSrc is not a capsAnalysis or CapsAnalysis class or a string name")   
        aimSrc = temp[:]
            
        checkAnalysisName(aimSrc)
             
        if len(aimSrc) != numDataSet:
            self.capsProblem.status = cCAPS.CAPS_MISMATCH
            self.capsProblem.checkStatus(msg="aimSrc dimensional mismatch between inputs!")
    
        # Check the field names in source analysis
        for var, src in zip(varName, aimSrc):
            
            fieldName = self.capsProblem.analysis[src].getAnalysisInfo(printInfo = False, infoDict = True)["fieldName"]

            if var not in fieldName:
                if not checkFieldWildCard(var, fieldName):
                    self.capsProblem.status = cCAPS.CAPS_BADVALUE
                    self.capsProblem.checkStatus(msg="Variable (field name) " + var + " not found in aimSrc "+ src)
        
        if aimDest:
            if not isinstance(aimDest, list):
                aimDest = [aimDest]
                
            temp = []
            for i in aimDest:
                if isinstance(i, capsAnalysis) or isinstance(i,CapsAnalysis):
                    temp.append(i.aimName)
                elif isinstance(i, str):
                    temp.append(i)
                else:    
                    raise TypeError("Provided aimDest is not a capsAnalysis or CapsAnalysis class or a string name")   
                
            aimDest = temp[:]
                
            checkAnalysisName(aimDest)
            
            if len(aimDest) != numDataSet:
                self.capsProblem.status = cCAPS.CAPS_MISMATCH
                self.capsProblem.checkStatus(msg="aimDest dimensional mismatch between inputs!")
                
            if transferMethod:
                
                if not isinstance(transferMethod, list):
                    transferMethod = [transferMethod]
                
                temp = []
                for i in transferMethod:
                    if "conserve" == i.lower():
                        temp.append(cCAPS.Conserve)
                        
                    elif "interpolate" == i.lower():
                        temp.append(cCAPS.Interpolate)
                    else:
                        print("Unreconized data transfer method - defaulting to \"Interpolate\"")
                        temp.append(cCAPS.Interpolate)
                
                transferMethod = temp[:]
                      
            else:
                print("No data transfer method (transferMethod) provided - defaulting to \"Interpolate\"")
                transferMethod =  [cCAPS.Interpolate]*numDataSet
                
            if len(transferMethod) != numDataSet:
                    self.capsProblem.status = cCAPS.CAPS_MISMATCH
                    self.capsProblem.checkStatus(msg="transferMethod dimensional mismatch between inputs!")
     
        if initValueDest:
            if not isinstance(initValueDest, list):
                initValueDest = [initValueDest]
                
            if len(initValueDest) != numDataSet:
                self.capsProblem.status = cCAPS.CAPS_MISMATCH
                self.capsProblem.checkStatus(msg="initValueDest dimensional mismatch between inputs!")
        else:
            if aimDest:
                initValueDest = [None]*numDataSet
            else:
                initValueDest = []

        # Create the object
        self.__createObj(capsBound, varName, aimSrc, aimDest, transferMethod, initValueDest)
        
        # Set the bound in the problem
        self.capsProblem.dataBound[capsBound] = self
       
    def __createObj(self, boundName, varName, aimSrc, aimDest, transferMethod, initValueDest):
        
        ## Bound/transfer name used to set up the data bound.
        self.boundName = boundName
        
        ## List of variables in the bound. 
        self.variables = varName
        
        ## Dictionary of vertex set object (\ref pyCAPS._capsVertexSet) in the bound. 
        self.vertexSet = {}
        
        ## Dictionary of "source" data set objects (\ref pyCAPS._capsDataSet) in the bound.
        self.dataSetSrc = {}
        
        ## Dictionary of "destination" data set objects (\ref pyCAPS._capsDataSet) in the bound.
        self.dataSetDest = {}
        
        if aimDest == None:
            aimDest = []
        if transferMethod == None:
            transferMethod = []
        if initValueDest == None:
            initValueDest = [None]*len(aimDest)

        # Open bound 
        boundName = _byteify(boundName)
        
        self.capsProblem.status = cCAPS.caps_makeBound(self.capsProblem.problemObj,
                                                       2, #int dim,
                                                       <char *> boundName,
                                                       &self.boundObj)
        self.capsProblem.checkStatus(msg = "while creating bound - " + str(self.boundName) 
                                     + " (during a call to caps_makeBound)")

        # Make all of our vertex sets 
        for i in set(aimSrc+aimDest):
            #print"Vertex Set - ", i
            self.vertexSet[i] = _capsVertexSet(i, 
                                               self.capsProblem, 
                                               self, 
                                               self.capsProblem.analysis[i])
        

        if len(varName) != len(aimSrc):
            print( "ERROR: varName and aimSrc must all be the same length")
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus()

        # Make all of our source data sets
        for var, src in zip(varName, aimSrc):
            #print"DataSetSrc - ", var, src
            
            rank = self.__getFieldRank(src, var)
                
            self.dataSetSrc[var] = _capsDataSet(var,
                                                self.capsProblem,
                                                self,
                                                self.vertexSet[src],
                                                cCAPS.Analysis, 
                                                rank)
        # Make destination data set
        if len(aimDest) > 0:

            if len(varName) != len(aimDest) or \
               len(varName) != len(transferMethod) or \
               len(varName) != len(aimSrc) or \
               len(varName) != len(initValueDest):
                print( "ERROR: varName, aimDest, transferMethod, aimSrc, and initValueDest must all be the same length")
                self.capsProblem.status = cCAPS.CAPS_BADVALUE
                self.capsProblem.checkStatus()

            for var, dest, method, src, initVal in zip(varName, aimDest, transferMethod, aimSrc, initValueDest):
                #print"DataSetDest - ", var, dest, method
            
                rank = self.__getFieldRank(src, var)
                
                self.dataSetDest[var] = _capsDataSet(var,
                                                     self.capsProblem,
                                                     self,
                                                     self.vertexSet[dest],
                                                     method, 
                                                     rank,
                                                     initVal)
            
        # Close bound 
        self.capsProblem.status = cCAPS.caps_closeBound(self.boundObj)
        self.capsProblem.checkStatus(msg = "while creating data bound - " + str(self.boundName) 
                                     + " (during a call to caps_closeBound)")
        
    cdef cCAPS.capsObj __boundObj(self):
        return self.boundObj
    
    # Get the field rank based on source and field name 
    def __getFieldRank(self, src, var):
        
        # Look for wild card fieldname and compare variable
        def __checkFieldWildCard(variable, fieldName):
            #variable = _byteify(variable)
            for i in fieldName:
                try:
                    i_ind = i.index("*")
                except ValueError:
                    i_ind = -1

                try:
                    if variable[:i_ind] == i[:i_ind]:  
                        return fieldName.index(i)
                except:
                    pass
                    
            return -1
        
        fieldName = self.capsProblem.analysis[src].getAnalysisInfo(printInfo=False, infoDict=True)["fieldName"]
        
        try: 
           
           index = fieldName.index(var)
        except:
            index  = __checkFieldWildCard(var, fieldName)
            if  index < 0:
                print("ERROR: Variable (field name)", var, "not found in aimSrc", src )
                self.capsProblem.status = cCAPS.CAPS_BADVALUE
                self.capsProblem.checkStatus()
                
        return self.capsProblem.analysis[src].getAnalysisInfo(printInfo=False, infoDict=True)["fieldRank"][index]
         
    ## Gets bound information for the bound object.
    #
    # \param printInfo Print information to sceen if True.
    # \param **kwargs See below.
    #
    # \return State of bound object.
    #
    # Valid keywords:
    # \param infoDict Return a dictionary containing bound information instead 
    # of just the state (default - False) 
    def getBoundInfo(self, printInfo=True, **kwargs):
        
        cdef:
            cCAPS.capsState state
            int dimension
            double limits[4]
            
            object boundInfo = {-2 : "multiple BRep entities – Error in reparameterization!y",
                                -1 : "Open",
                                 0 : "Empty and Closed", 
                                 1 : "Single BRep entity", 
                                 2 : "Multiple BRep entities"}

            
        self.capsProblem.status = cCAPS.caps_boundInfo(self.boundObj, 
                                                       &state, 
                                                       &dimension, 
                                                       limits)
        self.capsProblem.checkStatus(msg = "while getting data transfer information - " + str(self.boundName)
                                         + " (during a call to caps_boundInfo)")
        if printInfo:
            print("Bound Info:")
            print("\tName           = ", self.boundName)
            print("\tDimension      = ", dimension)
            # Add something for limits 
            
            if state in boundInfo:
                print("\tState          = ", state, " -> ", boundInfo[state] )
            else: 
                print("\tState          = ", state, " ",)
        
        if kwargs.pop("infoDict", False):
            return {"name"     : self.boundName, 
                    "dimension": dimension, 
                    "status"   : state
                    }
            
        return state

    ## Populates VertexSets for the bound. Must be called to finalize
    # the bound after all mesh generation aim's have been executed
    def fillVertexSets(self):
        
        cdef:
            int nErr
            cCAPS.capsErrs *errors

        self.capsProblem.status = cCAPS.caps_fillVertexSets(self.boundObj, 
                                                            &nErr, 
                                                            &errors)
        self.capsProblem.checkStatus(msg = "while fillin vertex sets - " + str(self.boundName)
                                         + " (during a call to caps_fillVertexSets)")

    ## Execute data transfer for the bound. 
    # \param variableName Name of variable to implement the data transfer for. If no name is provided 
    # the first variable in bound is used.
    def executeTransfer(self, variableName = None):
        
        if not self.dataSetDest:
            print("No destination data exists in bound!")
            return
        
        self.getDestData(variableName)


    def getDestData(self, variableName = None):
        
        if not self.dataSetDest:
            print("No destination data exists in bound!")
            return
        
        if variableName is None:
            print("No variable name provided - defaulting to variable - " + self.variables[0])
            variableName = self.variables[0]
       
        elif variableName not in self.variables:
            print("Unreconized varible name " + str(variableName))
            return

        self.dataSetDest[variableName].getData()
    
    ## Visualize data in the bound. 
    # The function currently relies on matplotlib or the capViewer class (webviewer) to plot the data. 
    # \param variableName Name of variable to visualize. If no name is provided 
    # the first variable in the bound is used.
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param viewerType What viewer should be used (default - capsViewer). Options: capsViewer or matplotlib (options are
    # case insensitive). Important: if $filename isn't set to None, the default to changed to matplotlib.
    # \param portNumber Port number to start the server listening on (default - 7681).
    # \param filename Save image(s) to file specified (default - None). Not available when using the webviewer
    # \param colorMap  Valid string for a, matplotlib::cm,  colormap (default - 'Blues'). Not as options are 
    # available when using the webviewer (see \ref capsViewer for additional details). 
    # \param showImage Show image(s) (default - True).
    def viewData(self, variableName=None, **kwargs):
        
        filename = kwargs.pop("filename", None)       
        if (filename is not None and "." not in filename):
            filename += ".png"
        
        viewer = kwargs.pop("viewerType", "capsViewer").lower()
        
        if variableName is None:
            print("No variable name provided - defaulting to variable - " + self.variables[0])
            variableName = self.variables[0]
           
        elif variableName not in self.variables:
         
            self.capsProblem.status = cCAPS.CAPS_NOTFOUND
            self.capsProblem.checkStatus(msg = "Unrecognized varible name " + str(variableName))
                
        # Are we going to use the viewer?       
        if filename is None and viewer == "capsviewer":
            
            # Initialize the viewer class
            viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))
            
            # Add src data set to viewer
            viewer.addDataSet(self.dataSetSrc[variableName])
            
            # Add dest data set to viewer
            if self.dataSetDest:
                viewer.addDataSet(self.dataSetDest[variableName])
            
            # Start up server 
            viewer.startServer()
            
            del viewer
        
        else:
        
            try:
                import matplotlib.pyplot as plt
            except:
                print("Error: Unable to import matplotlib - viewing the data is not possible")
                raise ImportError  
             

            
            showImage = kwargs.pop("showImage", True)
            kwargs["showImage"] =  False
            
            kwargs["filename"] = None
              
            fig = plt.figure() #figsize=(12,15)
            
            
            if self.dataSetDest:
            
                self.dataSetSrc[variableName].viewData(fig,  2, 0, **kwargs)
                self.dataSetDest[variableName].viewData(fig, 2, 1, **kwargs)
            
            else:
                self.dataSetSrc[variableName].viewData(fig,  1, 0, **kwargs)
            
            if filename is not None:
                print("Saving figure - ", filename)
                fig.savefig(filename)
                
            if showImage:
                plt.show()
            
            # Cleanup 
            fig.clf()
            del fig 

    ## Write a Tecplot compatiable file for the data in the bound. 
    # \param filename Name of file to save data to.
    # \param variableName Single or list of variables to write data for. If no name is provided 
    # all variables in the bound are used.
    def writeTecplot(self, filename, variableName=None):

        if filename is None:
            print("No filename provided")
            return
        
        if ("." not in filename):
            filename += ".dat"
        
        try:        
            file = open(filename, "w")
        except:
            self.capsProblem.status = cCAPS.CAPS_IOERR
            self.capsProblem.checkStatus(msg = "while open file for writing Tecplot file for data set - " 
                                         + str(self.dataSetName))

        if variableName is None:
            print("No variable name provided - all variables will be writen to file")
            variableName = self.variables
       
        if not isinstance(variableName,list):
            variableName = [variableName]
             
        for var in variableName:
            if var not in self.variables:
                print("Unreconized varible name " + str(var))
                continue
           
            self.dataSetSrc[var].writeTecplot(file)
            
            if self.dataSetDest:
            
                self.dataSetDest[var].writeTecplot(file)
            
        # Close the file 
        file.close()
        
    ## Create a HTML dendrogram/tree of the current state of the bound. The HTML file relies on 
    # the open-source JavaScript library, D3, to visualize the data. 
    # This library is freely available from https://d3js.org/ and is dynamically loaded within the HTML file. 
    # If running on a machine without internet access a (miniaturized) copy of the library may be written to 
    # a file alongside the generated HTML file by setting the internetAccess keyword to False. If set to True, internet 
    # access will be necessary to view the tree.  
    #
    # \param filename Filename to use when saving the tree (default - "boundName"). Note an ".html" is 
    # automatically appended to the name (same with ".json" if embedJSON = False).
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param embedJSON Embed the JSON tree data in the HTML file itself (default - True). If set to False 
    #  a seperate file is generated for the JSON tree data.
    #
    # \param internetAccess Is internet access available (default True)? If set to True internet access 
    # will be necessary to view the tree.  
    def createTree(self, filename = "boundName", **kwargs):
        
        if filename == "boundName":
            filename = self.boundName
            
        embedJSON = kwargs.pop("embedJSON", True)
        
        if embedJSON is False:
            try:
                # Write JSON data file 
                fp = open(filename + ".json",'w')
                fp.write(createTree(self, False, False, False))
                fp.close()
                
                # Write the HTML file
                self.capsProblem.status = writeTreeHTML(filename +".html", 
                                                        filename + ".json", 
                                                        None, 
                                                        kwargs.pop("internetAccess", True))
            except:
                self.capsProblem.status = cCAPS.CAPS_IOERR
        else: 
            self.capsProblem.status = writeTreeHTML(filename+".html", 
                                                    None, 
                                                    createTree(self, False, False, False),
                                                    kwargs.pop("internetAccess", True))
        
        self.capsProblem.checkStatus(msg = "while creating HTML tree")

## Functions to interact with a CAPS vertexSet object.
# Should be initiated within \ref pyCAPS.capsBound (not a standalone class)
cdef class _capsVertexSet:
    
    cdef:
        capsProblem capsProblem
        capsBound capsBound
        capsAnalysis capsAnalysis
        
        cCAPS.capsObj vertexSetObj
        readonly object vertexSetName
         
    def __cinit__(self, name, problemParent, boundParent, analysisParent):
        
        cdef:
            cCAPS.capsState state
            int dimension
            double limits[4]
        
        ## Vertex set name (analysis name the vertex belongs to).  
        self.vertexSetName = name  
        
        # Problem object the vertex set belongs to.  
        self.capsProblem  = problemParent
        
        ## Bound object (\ref pyCAPS.capsBound)the vertex set belongs to. 
        self.capsBound    = boundParent
        
        ## Analysis object (\ref pyCAPS.capsAnalysis) the vertex set belongs to. 
        self.capsAnalysis = analysisParent
        
        self.capsProblem.status = cCAPS.caps_makeVertexSet(self.capsBound.boundObj,
                                                           self.capsAnalysis.analysisObj,
                                                           <char *> NULL,
                                                           &self.vertexSetObj)
        self.capsProblem.checkStatus(msg = "while making vertex set - " + str(self.vertexSetName) 
                                     + " (during a call to caps_makeVertexSet)")

    cdef cCAPS.capsObj __vertexSetObj(self):
        return self.vertexSetObj

## Functions to interact with a CAPS dataSet object.
# Should be initiated within \ref pyCAPS.capsBound (not a standalone class)
cdef class _capsDataSet:
    
    cdef:
        capsProblem capsProblem
        
        readonly capsBound capsBound
        readonly _capsVertexSet capsVertexSet
        
        cCAPS.capsObj dataSetObj
        
        readonly object dataSetName
        readonly object dataSetMethod
         
        # Data values
        double *dataValue
        readonly int dataRank
        int dataNumPoint
        char *dataUnits
        cCAPS.capsdMethod dataMethod
        
        # Data xyz 
        int numNode 
        double *xyz
        char *xyzUnits
        
    def __cinit__(self, name, problemParent, boundParent, vertexSetParent, dataMethod, dataRank, initValue = None):
        
        ## Data set name (variable name). 
        self.dataSetName = name  
        
        # Reference to the problem object that data set pertains to.
        self.capsProblem = problemParent
        
        ## Reference to the bound object (\ref pyCAPS.capsBound) that data set pertains to.
        self.capsBound = boundParent
        
        ## Reference to the vertex set object (\ref pyCAPS._capsVertexSet) that data set pertains to.
        self.capsVertexSet = vertexSetParent
         
        ## Data method: Analysis, Interpolate, Conserve. 
        self.dataSetMethod = None
        
        self.dataMethod = dataMethod
        
        if self.dataMethod == cCAPS.Conserve:
            self.dataSetMethod = "Conserve"
        elif self.dataMethod == cCAPS.Interpolate:
            self.dataSetMethod = "Interpolate"
        elif self.dataMethod == cCAPS.Analysis:
            self.dataSetMethod = "Analysis"
        
        ## Rank of data set.
        self.dataRank = dataRank
        
        dataSetName = _byteify(self.dataSetName)
        self.capsProblem.status = cCAPS.caps_makeDataSet(self.capsVertexSet.vertexSetObj,
                                                         <char *> dataSetName,
                                                         self.dataMethod,
                                                         self.dataRank,
                                                         &self.dataSetObj)
        self.capsProblem.checkStatus(msg = "while making data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_makeDataSet)")

        # Set an initial constant value for the data set
        if initValue is not None:
            self.initDataSet(initValue)

    cdef cCAPS.capsObj __dataSetObj(self):
        return self.dataSetObj

    ## Executes caps_initDataSet on data set object to set an inital constant value needed for cyclic data transfers. 
    def initDataSet(self, initValue):
        
        cdef:
            int rank
            double *startup = NULL
            
        if (isinstance(initValue, list) or 
            isinstance(initValue, tuple)):
            rank = len(initValue)
        else:
            rank = 1
            initValue = [initValue]
            
        startup = <double *> malloc(rank*sizeof(double))
        for i in range(rank):
            startup[i] = initValue[i]

        self.capsProblem.status = cCAPS.caps_initDataSet(self.dataSetObj, 
                                                         rank, 
                                                         startup)
        free(startup)
        self.capsProblem.checkStatus(msg = "while setting initial values for data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_iniDataSet)")

    ## Executes caps_getData on data set object to retrieve data set variable, \ref dataSetName. 
    # \return Optionally returns a list of data values. Data with a rank greater than 1 returns a list of lists (e.g. 
    # data representing a displacement would return [ [Node1_xDisplacement, Node1_yDisplacement, Node1_zDisplacement], 
    # [Node2_xDisplacement, Node2_yDisplacement, Node2_zDisplacement], etc. ]   
    def getData(self):
        self.capsProblem.status = cCAPS.caps_getData(self.dataSetObj, 
                                                     &self.dataNumPoint, 
                                                     &self.dataRank, 
                                                     &self.dataValue,
                                                     &self.dataUnits)
        self.capsProblem.checkStatus(msg = "while getting data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_getData)")
        
        dataOut = []
        for i in range(self.dataNumPoint):
            
            if self.dataRank > 1:
                temp = []
                for j in range(self.dataRank):
                    temp.append(self.dataValue[self.dataRank*i + j])
            else:
                temp = self.dataValue[i]
            
            dataOut.append(temp)
                   
        return dataOut
    
    ## Executes caps_getData on data set object to retrieve XYZ coordinates of the data set. 
    # \return Optionally returns a list of lists of x,y, z values (e.g. [ [x2, y2, z2], [x2, y2, z2],  
    # [x3, y3, z3], etc. ] )      
    def getDataXYZ(self):
        
        cdef:
            int tempRank
            cCAPS.capsObj tempDataObj
            
        self.capsProblem.status = cCAPS.caps_childByName(self.capsVertexSet.vertexSetObj, 
                                                         cCAPS.DATASET, 
                                                         cCAPS.NONE,
                                                         <char *> "xyz", 
                                                         &tempDataObj)
        self.capsProblem.checkStatus(msg = "while getting XYZ associated with data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_childByName)")
        
        self.capsProblem.status = cCAPS.caps_getData(tempDataObj, 
                                                     &self.numNode, 
                                                     &tempRank, 
                                                     &self.xyz,
                                                     &self.xyzUnits)
        
        self.capsProblem.checkStatus(msg = "while getting XYZ associated with data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_getData)")
      
        return [[self.xyz[tempRank*i + j] for j in range(tempRank)] for i in range(self.numNode)]
    
    ## Executes caps_triangulate on data set's vertex set to retrieve the connectivity (triangles only) information 
    # for the data set. 
    # \return Optionally returns a list of lists of connectivity values 
    # (e.g. [ [node1, node2, node3], [node2, node3, node7], etc. ] ) and a list of lists of data connectivity (not this is    
    # an empty list if the data is node-based) (eg. [ [node1, node2, node3], [node2, node3, node7], etc. ]
    def getDataConnect(self):
        cdef:
            int numTri, numData
            int *triConn = NULL 
            int *dataConn = NULL
     
        self.capsProblem.status = cCAPS.caps_triangulate(self.capsVertexSet.vertexSetObj, 
                                                         &numTri, 
                                                         &triConn, 
                                                         &numData, 
                                                         &dataConn)
        self.capsProblem.checkStatus(msg = "while getting data connectivity associated with data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_triangulate)")
        
        connectivity = [[triConn[3*i + j] for j in range(3)] for i in range(numTri)]
        
        if numData == 0:
            dconnectivity = []
        else:
            dconnectivity = [[dataConn[3*i + j] for j in range(3)] for i in range(numData)]
            
        if triConn:
            cEGADS.EG_free(triConn)
        
        if dataConn:
            cEGADS.EG_free(dataConn)
            
        return (connectivity, dconnectivity)
        
    ## Visualize data set.
    # The function currently relies on matplotlib to plot the data. 
    # \param fig Figure object (matplotlib::figure) to append image to. 
    # \param numDataSet Number of data sets in \ref _capsDataSet.viewData$fig.  
    # \param dataSetIndex Index of data set being added to \ref _capsDataSet.viewData$fig.
    # \param **kwargs See below.
    #
    # Valid keywords: 
    # \param filename Save image(s) to file specified (default - None).  
    # \param colorMap  Valid string for a, matplotlib::cm,  colormap (default - 'Blues').
    # \param showImage Show image(s) (default - True).
    # \param title Set a custom title on the plot (default - VertexSet= 'name', DataSet = 'name', (Var. '#') ).
    def viewData(self, fig=None, numDataSet=1, dataSetIndex=0, **kwargs):
        
        #viewer = kwargs.pop("viewerType", "capsViewer").lower()
        
        # Are we going to use the viewer?       
        #if viewer == "capsviewer":
            
            # Initialize the viewer class
            #viewer = capsViewer(portNumber = kwargs.pop("portNumber", 7681))
            
            # Add src data set to viewer
            #viewer.addDataSet(self)
            
            # Start up server 
            #viewer.startServer()
            
            #del viewer
        
        #else:
            #self._viewDataMatplotLib(self,fig,numDataSet,dataSetIndex,**kwargs)

    #def _viewDataMatplotLib(self, fig=None, numDataSet=1, dataSetIndex=0, **kwargs):

        cdef:
            int i, j
            int numTri, numData
            int *triConn 
            int *dataConn

        try:
            import matplotlib.pyplot as plt
            from mpl_toolkits.mplot3d import Axes3D
            from mpl_toolkits.mplot3d.art3d import Poly3DCollection
           
            from matplotlib import colorbar
            from matplotlib import colors
            from matplotlib import cm as colorMap
        except:
            print("Error: Unable to import matplotlib - viewing the data is not possible")
            raise ImportError 
        
        filename = kwargs.pop("filename", None)       
        if (filename is not None and "." not in filename):
            filename += ".png"
            
        showImage =  kwargs.pop("showImage", True)
        title = kwargs.pop("title", None)
        
        if (filename is None and fig is None and showImage is False):
            print("No 'filename' or figure instance provided and 'showImage' is false. Nothing to do here!")
            return 
       
        cMap = kwargs.pop("colorMap", "Blues")
        
        try:
            colorMap.get_cmap(cMap)
        except ValueError:
            print("Colormap ",  cMap, "is not recognized. Defaulting to 'Blues'")
            cMap = "Blues"
            
        # Initialize values
        colorArray = []
        triConn = NULL
        dataConn = NULL
        
        self.getData()
        self.getDataXYZ()

        self.capsProblem.status = cCAPS.caps_triangulate(self.capsVertexSet.vertexSetObj, 
                                                         &numTri, 
                                                         &triConn, 
                                                         &numData, 
                                                         &dataConn)
        self.capsProblem.checkStatus(msg = "while viewing data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_triangulate)")
       
        if not ((self.dataNumPoint == self.numNode) or (self.dataNumPoint == numTri)):
            self.capsProblem.status = cCAPS.CAPS_MISMATCH
            self.capsProblem.checkStatus(msg = "viewData only supports node-centered or cell-centered data!")

        if numData != 0:
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "viewData does NOT currently support cell-centered based data sets!")
        
        if fig is None:
            fig = plt.figure() #figsize=(9,9)
        
        cellCenter = (self.dataNumPoint == numTri)
        
        ax = []
        for j in range(self.dataRank):
            
            if self.dataRank > 1 or numDataSet > 1:
                ax.append(fig.add_subplot(self.dataRank, numDataSet, numDataSet*j + 1 + dataSetIndex, projection='3d'))
            else:
                ax.append(fig.gca(projection='3d'))
            
            colorArray = []
            for i in range(self.dataNumPoint):
                colorArray.append(self.dataValue[self.dataRank*i + j])
            
            xLimits = [9E6, -9E6] # [minX, maxX] - initial values are arbitary
            yLimits = [9E6, -9E6]
            zLimits = [9E6, -9E6]
            
            minColor = min(colorArray)
            maxColor = max(colorArray)
            dColor   = maxColor - minColor
            
            p = 0.0
            if (maxColor == 0.0 and minColor == 0.0):
                zeroP = True
            else:
                zeroP = False
                
            for i in range (numTri):
                tri = []
                
                index = triConn[3*i+0]-1
                tri.append([self.xyz[3*index+0], self.xyz[3*index+1], self.xyz[3*index+2]])
                
                if zeroP == False and not cellCenter:
                    p = (colorArray[index] - minColor)/dColor
                
                index = triConn[3*i+1]-1
                tri.append([self.xyz[3*index+0], self.xyz[3*index+1], self.xyz[3*index+2]])
                
                if zeroP == False and not cellCenter:
                    p += (colorArray[index] - minColor)/dColor
                
                index = triConn[3*i+2]-1
                tri.append([self.xyz[3*index+0], self.xyz[3*index+1], self.xyz[3*index+2]])
                
                if zeroP == False and not cellCenter:
                    p += (colorArray[index] - minColor)/dColor
                
                if cellCenter:
                    p = 255*(colorArray[i] - minColor)/dColor
                else:
                    p = p/3.0*255
                    
                xLimits[0] = min(xLimits[0], tri[-1][0], tri[-2][0], tri[-3][0])
                xLimits[1] = max(xLimits[1], tri[-1][0], tri[-2][0], tri[-3][0])
                
                yLimits[0] = min(yLimits[0], tri[-1][1], tri[-2][1], tri[-3][1])
                yLimits[1] = max(yLimits[1], tri[-1][1], tri[-2][1], tri[-3][1])
                
                zLimits[0] = min(zLimits[0], tri[-1][2], tri[-2][2], tri[-3][2])
                zLimits[1] = max(zLimits[1], tri[-1][2], tri[-2][2], tri[-3][2])
                
                tri = Poly3DCollection([tri])
                 
                facecolor = colorMap.get_cmap(cMap)(int(p))
               
                tri.set_facecolor(facecolor)
                tri.set_edgecolor(facecolor)
                ax[-1].add_collection3d(tri)
            
            #ax[-1].autoscale_view()
            if not title:
                title = "VertexSet = " + str(self.capsVertexSet.vertexSetName) + \
                        "\nDataSet = " + str(self.dataSetName)
                if self.dataRank > 1:
                    title += " (Var. " + str(j+1) + ")"
                     
            ax[-1].set_title(title, size="small")
            ax[-1].set_axis_off() # Turn the axis off
            #ax[-1].invert_xaxis() # Invert x- axis
            ax[-1].view_init(azim=120, elev=30) # Isometric
            
            center = []
            center.append((xLimits[1] + xLimits[0])/2.0)
            center.append((yLimits[1] + yLimits[0])/2.0)
            center.append((zLimits[1] + zLimits[0])/2.0)
         
            scale = max(abs(xLimits[1] - xLimits[0]),
                        abs(yLimits[1] - yLimits[0]),
                        abs(zLimits[1] - zLimits[0]))      
         
            scale = 0.6*(scale/2.0)
             
            ax[-1].set_xlim(-scale + center[0], scale + center[0])
            ax[-1].set_ylim(-scale + center[1], scale + center[1])
            ax[-1].set_zlim(-scale + center[2], scale + center[2])
            
            cax, kw = colorbar.make_axes(ax[-1], location="right")
            norm = colors.Normalize(vmin=min(colorArray), vmax=max(colorArray))
            colorbar.ColorbarBase(cax, cmap=cMap, norm=norm)
            
            del colorArray[:]
            
        if triConn:
            cEGADS.EG_free(triConn)
        
        if dataConn:
            cEGADS.EG_free(dataConn)
        
        #fig.subplots_adjust(hspace=0.1) #left=0.0, bottom=0.0, right=1.0, top=.93, wspace=0.0 )
        
        if showImage:
            plt.show()
            
        if filename is not None:
            print("Saving figure - ", filename)
            fig.savefig(filename)
    
           
    ## Write data set to a Tecplot compatible data file. A triagulation of the data set will be used 
    # for the connectivity. 
    # \param file Optional open file object to append data to. If not provided a 
    # filename must be given via the keyword arguement \ref _capsDataSet.writeTecplot$filename.
    # \param filename Write Tecplot file with the specified name.  
    def writeTecplot(self, file=None,  filename=None):
        
        cdef:
            object xyz, data, connectivity, dconnectivity
   
        # If no file 
        if file == None:
            
            if filename == None:
                print("No file name or open file object provided")
                self.capsProblem.status = cCAPS.CAPS_BADVALUE
                self.capsProblem.checkStatus(msg = "while writing Tecplot file for data set - " + str(self.dataSetName) 
                                     + " (during a call to caps_triangulate)")
            
            if "." not in filename:
                filename += ".dat"
            
            try:        
                fp = open(filename, "w")
            except:
                self.capsProblem.status = cCAPS.CAPS_IOERR
                self.capsProblem.checkStatus(msg = "while open file for writing Tecplot file for data set - " 
                                             + str(self.dataSetName))
        else:
            fp = file
            
        data = self.getData()
        if not data: # Do we have empty data
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "the data set is empty for - " 
                                             + str(self.dataSetName))
        else:
            if not isinstance(data[0], list): # Only returns a list if rank 1 
                temp = []
                for d in data:
                    temp.append([d])
                data = temp
            
        xyz = self.getDataXYZ()
        
        connectivity, dconnectivity = self.getDataConnect()
     
        #TODO: dconnectivity does not appear to indicate a cell centered data set?
        if dconnectivity:
            self.capsProblem.status = cCAPS.CAPS_BADVALUE
            self.capsProblem.checkStatus(msg = "writeTecplot does NOT currently support dconnectivity based data sets!")
       
        if not ((len(data) == len(xyz)) or (len(data) == len(connectivity))):
            self.capsProblem.status = cCAPS.CAPS_MISMATCH
            self.capsProblem.checkStatus(msg = "writeTecplot only supports node-centered or cell-centered data!")
        
        title = 'TITLE = "VertexSet: ' + str(self.capsVertexSet.vertexSetName) + \
                ', DataSet = ' + str(self.dataSetName) + '"\n'
        fp.write(title)
        
        title = 'VARIABLES = "x", "y", "z"'
        for i in range(self.dataRank):
            title += ', "Variable_' + str(i) + '"'
        title += "\n"
        fp.write(title)

        numNode = self.numNode if len(data) == len(xyz) else 3*len(connectivity)
        
        title = "ZONE T = " + \
                '"VertexSet: ' + str(self.capsVertexSet.vertexSetName) + \
                ', DataSet = ' + str(self.dataSetName) + '"' +\
                 ", N = " + str(numNode) + \
                 ", E = " + str(len(connectivity)) + \
                 ", DATAPACKING = POINT, ZONETYPE= FETRIANGLE\n"
        fp.write(title)
       
        if len(data) == len(xyz): # Node centered data set
            # Write coordinates + data values 
            for i, j in zip(xyz, data):
            
                dataString = str(i[0]) +" "+ str(i[1]) +" "+ str(i[2])
                for k in j:
                    dataString += " " + str(k)
            
                fp.write(dataString + "\n")
        
            # Write connectivity
            for i in connectivity:
            
                dataString = str(i[0]) +" "+ str(i[1]) +" "+ str(i[2])
            
                fp.write(dataString + "\n")
                
        elif len(data) == len(connectivity): # Cell centered data set
            # Write coordinates + data values 
            for i, j in zip(connectivity, data):
           
                for t in i: 
                    dataString = str(xyz[t-1][0]) +" "+ str(xyz[t-1][1]) +" "+ str(xyz[t-1][2])
                    for k in j:
                        dataString += " " + str(k)
            
                    fp.write(dataString + "\n")
        
            # Write connectivity
            for i in xrange(len(connectivity)):
            
                dataString = str(3*i+1) +" "+ str(3*i+2) +" "+ str(3*i+3)
            
                fp.write(dataString + "\n")
        
        if file == None:
            fp.close()
            
## Functions to interact with a CAPS value object.
# Should be created with \ref capsProblem.createValue (not a standalone class)
cdef class _capsValue:
    
    cdef:
        readonly capsProblem capsProblem
        
        cCAPS.capsObj valueObj
    
        cCAPS.capssType subType # PARAMETER or USER - if creating the value from scratch
        readonly object name
        readonly object units 
    
        object _value
        object _limits
    
    def __cinit__(self, problemParent, 
                  name, value, # Use if creating a value from scratch, else set to None and call __setupObj
                  subType="parameter",
                  units=None, limits=None, 
                  fixedLength=True, fixedShape=True): 
                  
        
        cdef:
            cCAPS.capsvType valueType
            int numRow = 0
            int numCol = 0
            void *data = NULL
            
            int dim
            cCAPS.capsFixed lfixed, sfixed
           
        ## Reference to the problem object that loaded the value.
        self.capsProblem = problemParent
        
        # If the CAPS value object is provided we don't need to make it!
        if name is not None or value is not None :
            
            ## Variable name
            self.name = name
        
            ## Units of the variable
            self.units = units
        
            # Save subtype for quick access     
            if subType.lower() == "parameter":
                self.subType = cCAPS.PARAMETER 
            elif subType.lower() == "user":
                self.subType = cCAPS.USER
            else:
                print("Unreconized subType for value object!")
                self.capsProblem.status = cCAPS.CAPS_BADTYPE
                self.capsProblem.checkStatus(msg = "while creating value object - " + str(self.name))
        
            valueType = sortPythonObjType(value) 
            
            numRow, numCol = getPythonObjShape(value)
            
            byteValue  = []
            data = castValue2VoidP(value, byteValue, valueType)
            
            unitSys = _byteify(self.units)
            name = _byteify(self.name)
            
            try:
                if unitSys:
                    self.capsProblem.status = cCAPS.caps_makeValue(self.capsProblem.problemObj, 
                                                                   <char *> name,
                                                                   self.subType,
                                                                   valueType,
                                                                   numRow, 
                                                                   numCol,
                                                                   data,
                                                                   <char *> unitSys,
                                                                   &self.valueObj)
                                                              
                else:
                    self.capsProblem.status = cCAPS.caps_makeValue(self.capsProblem.problemObj, 
                                                                   <char *> name,
                                                                   self.subType,
                                                                   valueType, 
                                                                   numRow, 
                                                                   numCol,
                                                                   data,
                                                                   NULL,
                                                                   &self.valueObj)
                                                             
            finally:
                if isinstance(byteValue,list):
                    del byteValue[:]
                
                if data and valueType != cCAPS.String:
                    free(data)
             
            self.capsProblem.checkStatus(msg = "while creating value object - " + str(self.name) 
                                                 + " (during a call to caps_makeValue)")
            
            ## Acceptable limits for the value. Limits may be set directly. 
            #See \ref value3.py for a representative use case.
            self.limits = limits
            
            ## Value of the variable. Value may be set directly. 
            #See \ref value2.py for a representative use case.
            self.value = value # redundant - here just for the documentation purposes 
            
            # Get shape parameters
            if numRow == 1 and numCol == 1:
                dim = 0
            elif numRow == 1 or numCol ==1:
                dim = 1
            else:
                dim = 2
            
            if fixedLength:
                lfixed = cCAPS.Fixed
            else:
                lfixed = cCAPS.Change
            
            if fixedShape:
                sfixed = cCAPS.Fixed
            else:
                sfixed = cCAPS.Change
            
            # Set the shape of the value object  
            self.capsProblem.status = cCAPS.caps_setValueShape(self.valueObj, 
                                                               dim, 
                                                               lfixed, 
                                                               sfixed, 
                                                               cCAPS.NotNull)
            self.capsProblem.checkStatus(msg = "while setting the shape of value object - " + str(self.name) 
                                                + " (during a call to caps_setValueShape)")
        
    cdef __setupObj(self, cCAPS.capsObj valueObj):
        cdef:
            #  Get value infomation varianles 
            char *valueName = NULL
            cCAPS.capsoType objType
            cCAPS.capsObj link, parent, 
            cCAPS.capsOwn owner 
        
        self.valueObj = valueObj
            
        # Get the name of object
        self.capsProblem.status = cCAPS.caps_info(self.valueObj, 
                                                  &valueName, 
                                                  &objType, 
                                                  &self.subType, 
                                                  &link, &parent, &owner)
        self.capsProblem.checkStatus(msg = "while getting info on value object - " 
                                     + " (during a call to caps_info)")
            
        if valueName:
            self.name = <object> _strify(<char *> valueName)
        else:
            self.name = None
            
        # Get rest of value
        self.value
        
        return self
        
    def __del__(self):
        #print"Deleting value", self.name
        
        if self.subType == cCAPS.USER:
            self.capsProblem.status = cCAPS.caps_delete(self.valueObj)
            
            self.capsProblem.checkStatus(msg = "while deleting value object - " + str(self.name) 
                                             + " (during a call to caps_delete)")
        
    def __dealloc__(self):
        #print"Deallocating value", self.name
        
        if self.subType == cCAPS.USER:
            self.capsProblem.status = cCAPS.caps_delete(self.valueObj)
            
            self.capsProblem.checkStatus(msg = "while deleting value object - " + str(self.name) 
                                             + " (during a call to caps_delete)")
   
    cdef cCAPS.capsObj __valueObj(self):
        return self.valueObj
 
    @property 
    def value(self):
        cdef:
            
            cCAPS.capsvType vtype
            int vlen
            const void *data
        
            const char *valunits

            int nErr
            cCAPS.capsErrs *errors
            
            int vdim, vnrow, vncol
            cCAPS.capsFixed vlfixed, vsfixed
            cCAPS.capsNull  vntype
            
        self.capsProblem.status = cCAPS.caps_getValue(self.valueObj, 
                                                      &vtype, 
                                                      &vlen, 
                                                      &data, 
                                                      &valunits, 
                                                      &nErr, 
                                                      &errors)

        # Check errors
        report_Errors(nErr, errors)
        
        self.capsProblem.checkStatus(msg = "while getting value of value object - " + str(self.name) 
                                     + " (during a call to caps_getValue)")
        
        if valunits:
            self.units = <object> _strify(<char *> valunits)
        else:
            self.units = None
                
        self.capsProblem.status = cCAPS.caps_getValueShape(self.valueObj,
                                                           &vdim,
                                                           &vlfixed,
                                                           &vsfixed,
                                                           &vntype,
                                                           &vnrow,
                                                           &vncol)
        self.capsProblem.checkStatus(msg = "while getting shape of value object - " + str(self.name) 
                                             + " (during a call to caps_getValueShape)")
        
        if vntype == cCAPS.IsNull:
            self._value = None
        else:
#             if vnrow == 0: # Possibly remove after changes are made to getValueShape for attributes
#                 if vtype == cCAPS.String:
#                     vnrow = 1
#                 else:
#                     vnrow = vlen
#             if vncol == 0:
#                 vncol = 1
            
            self._value = castValue2PythonObj(data, vtype, vnrow, vncol)
        
        return self._value
    
    @value.setter
    def value(self, dataValue):
        
        cdef: 
            cCAPS.capsvType valueType
            
            int numRow, numCol
            void *data 
        
        numRow, numCol = getPythonObjShape(dataValue)
        
        valueType = sortPythonObjType(self.value)
       
        byteValue  = []
        data = castValue2VoidP(dataValue,
                               byteValue,
                               valueType)
        
        try:
            self.capsProblem.status = cCAPS.caps_setValue(self.valueObj, 
                                                          <int> numRow, 
                                                          <int> numCol, 
                                                          <const void *> data)
        finally:
            if isinstance(byteValue,list):
                del byteValue[:]
            
            if data and valueType != cCAPS.String:
                free(data)
                
        self.capsProblem.checkStatus(msg = "while setting analysis variable - " + str(self.name) 
                                         + " (during a call to caps_setValue)")
        
        self._value = dataValue
    
    @property
    def limits(self):
        
        cdef: 
            const void *limits = NULL
            cCAPS.capsvType valueType
        
        self.capsProblem.status = cCAPS.caps_getLimits(self.valueObj, &limits)
        self.capsProblem.checkStatus(msg = "while getting limits of value object - " + str(self.name) 
                                     + " (during a call to caps_getLimits)")
        
        valueType = sortPythonObjType(self.value)
        
        self._limits = castValue2PythonObj(limits, valueType, <int> 2, <int> 1)
        
        return self._limits
    
    @limits.setter
    def limits(self, limitsValue):
        
        cdef:
            void *limits = NULL
            
            cCAPS.capsvType valueType
        
        self._limits = None
        
        if limitsValue:
        
            self.capsProblem.status = cCAPS.CAPS_SUCCESS
            
            # Make sure limitValue is the right shape and dimension
            if not isinstance(limitsValue, list):
                print("Limits value should be 2 element list - [min value, max value]!")
                self.capsProblem.status = cCAPS.CAPS_BADVALUE
            else:
                if len(limitsValue) != 2:
                    print("Limits value should be 2 element list - [min value, max value]!")
                    self.capsProblem.status = cCAPS.CAPS_BADVALUE
                
                for i in limitsValue:
                    if not isinstance(i, (int, float)):
                        print("Invalid element type for limits value, only integer or float values are valid!")
                        self.capsProblem.status = cCAPS.CAPS_BADVALUE
                    
            self.capsProblem.checkStatus(msg = "while checking the limit of value object - " + str(self.name))
            
            valueType = sortPythonObjType(self.value) 
            
            byteValue  = []
            limits = castValue2VoidP(limitsValue, byteValue, valueType)
           
            try:
                self.capsProblem.status = cCAPS.caps_setLimits(self.valueObj, limits)
                
            finally:
                if isinstance(byteValue,list):
                    del byteValue[:]
                
                if limits and valueType != cCAPS.String:
                    free(limits)
                    
            self.capsProblem.checkStatus(msg = "while setting limits of value object - " + str(self.name) 
                                         + " (during a call to caps_setLimits)")
    ## \example value2.py
    # Example use cases for interacting the a pyCAPS._capsValue object for setting the value.
    
    ## \example value3.py
    # Example use cases for interacting the a pyCAPS._capsValue object for setting the limits.
    
    
    ## Change the value of the object. See \ref value2.py for a representative use case.
    # \param data Data value(s) for the variable. Note that data will be type casted to match the type 
    # used to original create the capsValue object.
    def setVal(self, data):
        self.value = data
    
    ## Get the current value of object. See \ref value2.py for a representative use case.
    # \return Current value of set for object. 
    def getVal(self):
        return self.value
    
    ## Set new limits. See \ref value3.py for a representative use case.
    # \param newLimits New values to set for the limits. Should be 2 element list - [min value, max value].
    def setLimits(self, newLimits):
        self.limits = newLimits
    
    ## Get the current value for the limits. See \ref value3.py for a representative use case.
    # \return Current value for the limits.
    def getLimits(self):
        return self.limits
    
    ## Return the current value of the object in the desired, specified units. Note that this neither
    # changes the value or units of the object, only returns a converted value. 
    # See \ref value4.py for a representative use case.
    # \return Current value of the object in the specified units.  
    def convertUnits(self, toUnits):
        ## \example value4.py
        # Example use cases for interacting the with pyCAPS._capsValue.convertUnits() function.
        
        if self.units == None:
            print("Value doesn't have defined units")
            self.capsProblem.status = cCAPS.CAPS_UNITERR
            self.capsProblem.checkStatus(msg = "while converting units for value object - " + str(self.name))
            
        
        try:
            data = capsConvert(self.value, self.units, toUnits)
        except:
            self.capsProblem.status = cCAPS.CAPS_UNITERR
            self.capsProblem.checkStatus(msg = "while converting units for value object - " + str(self.name) 
                                         + " (during a call to pyCAPS.capsConvert)")
        
        if data == None:
            self.capsProblem.status = cCAPS.CAPS_UNITERR
            self.capsProblem.checkStatus(msg = "while converting units for value object - " + str(self.name) 
                                         + " (during a call to pyCAPS.capsConvert)")
            
        return data

## Pythonic alias to capsProblem
cdef class CapsProblem(capsProblem):
    pass

## Pythonic alias to capsAnalysis
cdef class CapsAnalysis(capsAnalysis):
    pass

## Pythonic alias to capsGeometry
cdef class CapsGeometry(capsGeometry):
    pass

## Pythonic alias to capsBound
cdef class CapsBound(capsBound):
    pass
