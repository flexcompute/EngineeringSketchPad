#----------------------------------------------------------------
# Generated CMake target import file for configuration "Release".
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Import target "TKService" for configuration "Release"
set_property(TARGET TKService APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TKService PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libTKService.so.7.7.0"
  IMPORTED_SONAME_RELEASE "libTKService.so.7.7"
  )

list(APPEND _IMPORT_CHECK_TARGETS TKService )
list(APPEND _IMPORT_CHECK_FILES_FOR_TKService "${_IMPORT_PREFIX}/lib/libTKService.so.7.7.0" )

# Import target "TKV3d" for configuration "Release"
set_property(TARGET TKV3d APPEND PROPERTY IMPORTED_CONFIGURATIONS RELEASE)
set_target_properties(TKV3d PROPERTIES
  IMPORTED_LOCATION_RELEASE "${_IMPORT_PREFIX}/lib/libTKV3d.so.7.7.0"
  IMPORTED_SONAME_RELEASE "libTKV3d.so.7.7"
  )

list(APPEND _IMPORT_CHECK_TARGETS TKV3d )
list(APPEND _IMPORT_CHECK_FILES_FOR_TKV3d "${_IMPORT_PREFIX}/lib/libTKV3d.so.7.7.0" )

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
