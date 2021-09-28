#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

# C-Imports 
cimport cUDUNITS
    
## Convert a value from one unit to another using the UDUNITS-2 library.
# See \ref units.py for example use cases.
# Note that UDUNITS-2 is packaged with CAPS, so no additional dependencies 
# are necessary. Please refer to the UDUNITS-2 documentation for 
# specifics regarding the syntax for valid unit strings,  
# [UDUNITS-2 Manual](http://www.unidata.ucar.edu/software/udunits/udunits-current/doc/udunits/udunits2.html)
# or [NIST Units](http://physics.nist.gov/cuu/Units/index.html). The following table 
# was taken from the UDUNITS-2 manual to assist in some of the unit specifications. 
# 
# String Type | Using Names | Using Symbols | Comment|
# ------------|-------------|---------------|---------|
# Simple  |     meter       |      m   | |
# Raised  |     meter^2     |      m2  | higher precedence than multiplying or dividing|
# Product |     newton meter|      N.m | |
# Quotient|     meter per second | m/s | |
# Scaled  |     60 second   |      60 s| |
# Prefixed|     kilometer   |      km  | |
# Offset  |     kelvin from 273.15   |  K @ 273.15  | lower precedence than multiplying or dividing|
# Logarithmic |    lg(re milliwatt)  |  lg(re mW)   | "lg" is base 10, "ln" is base e, and "lb" is base 2|
# Grouped |    (5 meter)/(30 second) |  (5 m)/(30 s)| |
#
# \param value Input value to convert. Value may be an integer, float/double, or list of 
# integers and floats/doubles. Note that integers are automatically cast to floats/doubles
# \param fromUnits Current units of the input value (see [UDUNITS-2 Manual] for valid string). 
# \param toUnits Desired units to convert the input value to (see [UDUNITS-2 Manual] for valid string).  
# \param ignoreWarning Ignore UDUNITS verbose warnings (default - True). Errors during 
# unit conversions are still reported.
#
# \return Return the input value(s) in the specified units.    
#
# \exception TypeError Wrong type [float(s)/double(s)/integer(s)] for value or non-string value for from/to-Units.
# \exception ValueError Error during unit conversion. See raised message for additional details.
cpdef capsConvert(value, fromUnits, toUnits, ignoreWarning = True):
    ## \example units.py
    # Example use cases for pyCAPS.capsConvert() function.
    cdef: 
        cUDUNITS.ut_system *unitSys=NULL
        
        cUDUNITS.ut_unit *unitIn=NULL
        cUDUNITS.ut_unit *unitOut=NULL
        
        cUDUNITS.cv_converter *converter=NULL
       
        double tempValue
        
        int status
    
    if not isinstance(value, (list, float, int)):
        #print "Value should either be an integer, float/double or a list of integers/floats/doubles"
        #print "Note that integers will be converted to floats/doubles"
        raise TypeError("Value should either be an integer, float/double or a list of integers/floats/doubles")
    
    if fromUnits == None or not isinstance(fromUnits, str):
        #print "No valid units specified for fromUnits!"
        raise TypeError("No valid units specified for fromUnits!")
        #return None
    
    if toUnits == None or not isinstance(toUnits, str):
        #print "No valid units specified for toUnits!"
        raise TypeError("No valid units specified for toUnits!")
        #return None
    
    if ignoreWarning:
        <void> cUDUNITS.ut_set_error_message_handler(cUDUNITS.ut_ignore)
        
    # Set default value
    valueOut = None
        
    unitSys = cUDUNITS.ut_read_xml(NULL)
    if unitSys == NULL:
        
        status = <int> cUDUNITS.ut_get_status()
        
        if status == <int> cUDUNITS.UT_OPEN_ENV:
            msg = "Environment variable UDUNITS2_XML_PATH is set but file couldn't be opened!"
        
        elif status == <int> cUDUNITS.UT_OPEN_DEFAULT:
            msg = "Environment variable UDUNITS2_XML_PATH is unset, and the " \
                  + "installed, default, unit database couldn't be opened!"
        
        elif status == <int> cUDUNITS.UT_PARSE:
            msg = "Couldn't parse unit database!"
        
        elif status == <int> cUDUNITS.UT_OS:
            msg = "Operating-system error!"
   
        print("Unable to create new unit sytem!")
        
        raise ValueError(msg)
        #return None
    
    # Byteify strings 
    fromUnits = _byteify(fromUnits)
    toUnits = _byteify(toUnits)
    
    unitIn  = cUDUNITS.ut_parse(unitSys, <const char *> fromUnits, cUDUNITS.UT_UTF8)
    unitOut = cUDUNITS.ut_parse(unitSys, <const char *> toUnits, cUDUNITS.UT_UTF8)
    
    status  = cUDUNITS.ut_are_convertible(unitOut, unitIn)
    if status == <int> 0:
        
        status = <int> cUDUNITS.ut_get_status()
      
        if status == <int> cUDUNITS.UT_BAD_ARG:
            msg = "fromUnits or toUnits is NULL!"
            
        elif status == <int> cUDUNITS.UT_NOT_SAME_SYSTEM: 
             msg = "fromUnits and toUnits belong to different unit-sytems!"
        
        elif status == <int> cUDUNITS.UT_SUCCESS:
            msg = "Conversion between the units (" +str(fromUnits) + " and " +str(toUnits) \
                  + ") is meaningless/not possible (e.g. 'fromUnits' is 'meter' " \
                  + "and 'toUnits' is 'kilogram')!"
    
        print("Unit conversion failure!") # Add better information here
        
        if unitIn:
            cUDUNITS.ut_free(unitIn)
        
        if unitOut:
            cUDUNITS.ut_free(unitOut)
        
        if unitSys:
            cUDUNITS.ut_free_system(unitSys)
            
        raise ValueError(msg)
        #return None
    
    converter = cUDUNITS.ut_get_converter( unitIn, unitOut)
    if converter == NULL:
        
        status = <int> cUDUNITS.ut_get_status()
        
        if status == <int> cUDUNITS.UT_BAD_ARG:
            msg = "fromUnits or toUnits is NULL!"
            
        elif status == <int> cUDUNITS.UT_NOT_SAME_SYSTEM: 
             msg = "fromUnits and toUnits belong to different unit-sytems!"
        
        elif status == <int> cUDUNITS.UT_MEANINGLESS:
             msg = "Conversion between the units (" +str(fromUnits) + " and " +str(toUnits) \
                 + ") is meaningless/not possible (e.g. 'fromUnits' is 'meter' " \
                  + "and 'toUnits' is 'kilogram')!"
    
        print("Issue getting converter!")
        
        if unitIn:
            cUDUNITS.ut_free(unitIn)
        
        if unitOut:
            cUDUNITS.ut_free(unitOut)
        
        if unitSys:
            cUDUNITS.ut_free_system(unitSys)
        
        raise ValueError(msg)
        #return None
    
    try:
        if isinstance(value, list):
            valueOut = []
            
            for i in value:
                
                if isinstance(i, list):
                    valueOut2 = []
                    for j in i:
                        j = float(j)
                        
                        tempDouble = cUDUNITS.cv_convert_double(converter, <double> j)
                        valueOut2.append(<object> tempDouble)
                    
                    valueOut.append(valueOut2)
                else:
                    i = float(i)
                    
                    tempDouble = cUDUNITS.cv_convert_double(converter, <double> i)
                    valueOut.append(<object> tempDouble)
        else:
            value = float(value)
            tempDouble = cUDUNITS.cv_convert_double(converter, <double> value)
            valueOut = <object> tempDouble
            
    finally:
        if unitIn:
            cUDUNITS.ut_free(unitIn)
        
        if unitOut:
            cUDUNITS.ut_free(unitOut)
        
        if converter:
            cUDUNITS.cv_free(converter)
        
        if unitSys:
            cUDUNITS.ut_free_system(unitSys)
        
    #Convert memory cleanups to a closure once support in cpdef
    
    # Return value
    return valueOut

# Convert values based on units from capsValue object. Returns a capsStatus and the data value
cdef object convertUnitPythonObj(cCAPS.capsObj valueObj,
                                 cCAPS.capsvType valueType,
                                 object dataUnit,
                                 object data):
    
    cdef:
        double valueOut
        int status 
        
    if not dataUnit: # No units - thats ok return 
        return cCAPS.CAPS_SUCCESS, data
    
    dataUnit = _byteify(dataUnit)
    
    if valueType != cCAPS.Double and valueType != cCAPS.Integer: # We only want to convert doubles and integers
        return cCAPS.CAPS_SUCCESS, data
    
    numRow, numCol = getPythonObjShape(data)
    
    for i in range(numCol): # Loop data through matrix 
        
        for j in range(numRow):
        
            if numCol == 1 and numRow == 1:
                
                status = cCAPS.caps_convert(valueObj, <char *> dataUnit, <double> data, NULL, &valueOut)
                if status != cCAPS.CAPS_SUCCESS: 
                    return status, data
                
                data = <object> valueOut
                
            elif numCol == 1:
                
                status = cCAPS.caps_convert(valueObj, <char *> dataUnit, <double> data[j], NULL, &valueOut)
                if status != cCAPS.CAPS_SUCCESS: 
                    return  status, data
                
                data[j] = <object> valueOut
                
            else:
                
                status = cCAPS.caps_convert(valueObj, <char *> dataUnit, <double> data[j][i], NULL, &valueOut)
                if status != cCAPS.CAPS_SUCCESS: 
                    return status, data
                
                data[j][i] = <object> valueOut
    
    return cCAPS.CAPS_SUCCESS, data