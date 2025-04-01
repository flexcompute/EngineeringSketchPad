#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TKExpress" for configuration "Release"
set_property(TARGET TKExpress APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TKExpress PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libTKExpress.so.7.8.1"
  IMPORTED_SONAME_RELEASE "libTKExpress.so.7.8"
  )

list(APPEND _IMPORT_CHECK_TARGETS TKExpress )
list(APPEND _IMPORT_CHECK_FILES_FOR_TKExpress "${_IMPORT_PREFIX}/lib/libTKExpress.so.7.8.1" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
