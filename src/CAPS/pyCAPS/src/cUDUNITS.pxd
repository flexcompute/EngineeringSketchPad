#
# Written by Dr. Ryan Durscher AFRL/RQVC
# 
# This software has been cleared for public release on 27 Oct. 2020, case number 88ABW-2020-3328.

cdef extern from "stdarg.h":
    ctypedef struct va_list

cdef extern from "udunits.h":
    
    ctypedef struct ut_system
    ctypedef union ut_unit
    ctypedef union cv_converter
    
    cdef enum utEncoding:
        UT_ASCII, UT_ISO_8859_1, UT_LATIN1 = UT_ISO_8859_1, UT_UTF8 = 2
        
    ctypedef utEncoding ut_encoding

    ctypedef int (*ut_error_message_handler)(const char* fmt, va_list args)
    
    cdef enum utStatus:
        UT_SUCCESS = 0,    # Success 
        UT_BAD_ARG,        # An argument violates the function's contract 
        UT_EXISTS,         # Unit, prefix, or identifier already exists 
        UT_NO_UNIT,        # No such unit exists 
        UT_OS,             # Operating-system error.  See "errno". 
        UT_NOT_SAME_SYSTEM,# The units belong to different unit-systems 
        UT_MEANINGLESS,    # The operation on the unit(s) is meaningless 
        UT_NO_SECOND,      # The unit-system doesn't have a unit named "second" 
        UT_VISIT_ERROR,    # An error occurred while visiting a unit 
        UT_CANT_FORMAT,    # A unit can't be formatted in the desired manner 
        UT_SYNTAX,         # string unit representation contains syntax error 
        UT_UNKNOWN,        # string unit representation contains unknown word 
        UT_OPEN_ARG,       # Can't open argument-specified unit database 
        UT_OPEN_ENV,       # Can't open environment-specified unit database 
        UT_OPEN_DEFAULT,   # Can't open installed, default, unit database 
        UT_PARSE           # Error parsing unit specification 

    ctypedef utStatus ut_status


    # Returns the binary representation of a unit corresponding to a string
    # representation.
    #
    # Arguments:
    #    system        Pointer to the unit-system in which the parsing will
    #            occur.
    #    string        The string to be parsed (e.g., "millimeters").  There
    #            should be no leading or trailing whitespace in the
    #            string.  See ut_trim().
    #    encoding    The encoding of "string".
    # Returns:
    #    NULL        Failure.  "ut_get_status()" will be one of
    #                UT_BAD_ARG        "system" or "string" is NULL.
    #                UT_SYNTAX        "string" contained a syntax
    #                        error.
    #                UT_UNKNOWN        "string" contained an unknown
    #                        identifier.
    #                UT_OS        Operating-system failure.  See
    #                        "errno".
    #    else        Pointer to the unit corresponding to "string".
    ut_unit* ut_parse( const ut_system* const    system,
                       const char* const         string,
                       const ut_encoding         encoding)
    
    # Returns a converter of numeric values in one unit to numeric values in
    # another unit.  The returned converter should be passed to cv_free() when it is
    # no longer needed by the client.
    #
    # NOTE:  Leap seconds are not taken into account when converting between
    # timestamp units.
    #
    # Arguments:
    #    from        Pointer to the unit from which to convert values.
    #    to        Pointer to the unit to which to convert values.
    # Returns:
    #    NULL        Failure.  "ut_get_status()" will be:
    #                UT_BAD_ARG        "from" or "to" is NULL.
    #                UT_NOT_SAME_SYSTEM    "from" and "to" belong to
    #                        different unit-systems.
    #                UT_MEANINGLESS    Conversion between the units is
    #                        not possible.  See
    #                        "ut_are_convertible()".
    #    else        Pointer to the appropriate converter.  The pointer
    #            should be passed to cv_free() when no longer needed by
    #            the client.
    cv_converter* ut_get_converter( ut_unit* const    fromVal,
                                    ut_unit* const    toVal)
    
    # Converts a double.
    # ARGUMENTS:
    #    converter    The converter.
    #    value        The value to be converted.
    # RETURNS:
    #    The converted value.
    double cv_convert_double( const cv_converter*    converter,
                              const double    value)
    
    # Indicates if numeric values in one unit are convertible to numeric values in
    # another unit via "ut_get_converter()".  In making this determination, 
    # dimensionless units are ignored.
    #
    # Arguments:
    #    unit1        Pointer to a unit.
    #    unit2        Pointer to another unit.
    # Returns:
    #    0        Failure.  "ut_get_status()" will be
    #                    UT_BAD_ARG        "unit1" or "unit2" is NULL.
    #                UT_NOT_SAME_SYSTEM    "unit1" and "unit2" belong to
    #                        different unit-sytems.
    #                UT_SUCCESS        Conversion between the units is
    #                        not possible (e.g., "unit1" is
    #                        "meter" and "unit2" is
    #                        "kilogram").
    #    else    Numeric values can be converted between the units.
    int ut_are_convertible( const ut_unit* const    unit1,
                            const ut_unit* const    unit2)
    
    # Frees resources associated with a unit.  This function should be invoked on
    # all units that are no longer needed by the client.  Use of the unit upon
    # return from this function will result in undefined behavior.
    #
    # Arguments:
    #    unit    Pointer to the unit to have its resources freed or NULL.
    void ut_free( ut_unit* const unit) # "udunits2.h":
    
    # Frees resources associated with a converter.
    # ARGUMENTS:
    #    conv    The converter to have its resources freed or NULL.
    void cv_free( cv_converter* const conv) # "converter.h"
    
    # Returns the unit-system corresponding to an XML file.  This is the usual way
    # that a client will obtain a unit-system.
    #
    # Arguments:
    #    path    The pathname of the XML file or NULL.  If NULL, then the
    #        pathname specified by the environment variable UDUNITS2_XML_PATH
    #        is used if set; otherwise, the compile-time pathname of the
    #        installed, default, unit database is used.
    # Returns:
    #    NULL    Failure.  "ut_get_status()" will be
    #            UT_OPEN_ARG        "path" is non-NULL but file couldn't be
    #                    opened.  See "errno" for reason.
    #            UT_OPEN_ENV        "path" is NULL and environment variable
    #                    UDUNITS2_XML_PATH is set but file
    #                    couldn't be opened.  See "errno" for
    #                    reason.
    #            UT_OPEN_DEFAULT    "path" is NULL, environment variable
    #                    UDUNITS2_XML_PATH is unset, and the
    #                    installed, default, unit database
    #                    couldn't be opened.  See "errno" for
    #                    reason.
    #            UT_PARSE        Couldn't parse unit database.
    #            UT_OS        Operating-system error.  See "errno".
    #    else    Pointer to the unit-system defined by "path".
    ut_system* ut_read_xml( const char*    path)

    # Frees a unit-system.  All unit-to-identifier and identifier-to-unit mappings
    # will be removed.
    #
    # Arguments:
    #    system        Pointer to the unit-system to be freed.  Use of "system"
    #            upon return results in undefined behavior.
    void ut_free_system( ut_system*    system)
    
    # Returns the previously-installed error-message handler and optionally
    # installs a new handler.  The initial handler is "ut_write_to_stderr()".
    #
    # Arguments:
    #      handler        NULL or pointer to the error-message handler.  If NULL,
    #            then the handler is not changed.  The 
    #            currently-installed handler can be obtained this way.
    # Returns:
    #    Pointer to the previously-installed error-message handler.
    ut_error_message_handler ut_set_error_message_handler( ut_error_message_handler handler)
    
    # Does nothing with an error-message.
    #
    # Arguments:
    #    fmt    The format for the error-message.
    #    args    The arguments of "fmt".
    # Returns:
    #    0    Always.
    int ut_ignore( const char* const    fmt, va_list  args)
    
    # Returns the status of the last operation by the units module.  This function
    # will not change the status.
    ut_status ut_get_status()
