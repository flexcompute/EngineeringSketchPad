# Generated by CMake

if("${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION}" LESS 2.5)
   message(FATAL_ERROR "CMake >= 2.6.0 required")
endif()
cmake_policy(PUSH)
cmake_policy(VERSION 2.6)
#----------------------------------------------------------------
# Generated CMake target import file.
#----------------------------------------------------------------

# Commands may need to know the format version.
set(CMAKE_IMPORT_FILE_VERSION 1)

# Protect against multiple inclusion, which would fail when already imported targets are added once more.
set(_targetsDefined)
set(_targetsNotDefined)
set(_expectedTargets)
foreach(_expectedTarget TKCDF TKLCAF TKCAF TKBinL TKXmlL TKBin TKXml TKStdL TKStd TKTObj TKBinTObj TKXmlTObj TKVCAF)
  list(APPEND _expectedTargets ${_expectedTarget})
  if(NOT TARGET ${_expectedTarget})
    list(APPEND _targetsNotDefined ${_expectedTarget})
  endif()
  if(TARGET ${_expectedTarget})
    list(APPEND _targetsDefined ${_expectedTarget})
  endif()
endforeach()
if("${_targetsDefined}" STREQUAL "${_expectedTargets}")
  unset(_targetsDefined)
  unset(_targetsNotDefined)
  unset(_expectedTargets)
  set(CMAKE_IMPORT_FILE_VERSION)
  cmake_policy(POP)
  return()
endif()
if(NOT "${_targetsDefined}" STREQUAL "")
  message(FATAL_ERROR "Some (but not all) targets in this export set were already defined.\nTargets Defined: ${_targetsDefined}\nTargets not yet defined: ${_targetsNotDefined}\n")
endif()
unset(_targetsDefined)
unset(_targetsNotDefined)
unset(_expectedTargets)


# Compute the installation prefix relative to this file.
get_filename_component(_IMPORT_PREFIX "${CMAKE_CURRENT_LIST_FILE}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
get_filename_component(_IMPORT_PREFIX "${_IMPORT_PREFIX}" PATH)
if(_IMPORT_PREFIX STREQUAL "/")
  set(_IMPORT_PREFIX "")
endif()

# Create imported target TKCDF
add_library(TKCDF SHARED IMPORTED)

set_target_properties(TKCDF PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKernel"
)

# Create imported target TKLCAF
add_library(TKLCAF SHARED IMPORTED)

set_target_properties(TKLCAF PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel"
)

# Create imported target TKCAF
add_library(TKCAF SHARED IMPORTED)

set_target_properties(TKCAF PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKernel;TKGeomBase;TKBRep;TKTopAlgo;TKMath;TKG2d;TKG3d;TKCDF;TKLCAF;TKBO"
)

# Create imported target TKBinL
add_library(TKBinL SHARED IMPORTED)

set_target_properties(TKBinL PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKLCAF"
)

# Create imported target TKXmlL
add_library(TKXmlL SHARED IMPORTED)

set_target_properties(TKXmlL PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKMath;TKLCAF"
)

# Create imported target TKBin
add_library(TKBin SHARED IMPORTED)

set_target_properties(TKBin PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKBRep;TKMath;TKernel;TKG2d;TKG3d;TKCAF;TKCDF;TKLCAF;TKBinL"
)

# Create imported target TKXml
add_library(TKXml SHARED IMPORTED)

set_target_properties(TKXml PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKMath;TKBRep;TKG2d;TKGeomBase;TKG3d;TKLCAF;TKCAF;TKXmlL"
)

# Create imported target TKStdL
add_library(TKStdL SHARED IMPORTED)

set_target_properties(TKStdL PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKernel;TKCDF;TKLCAF"
)

# Create imported target TKStd
add_library(TKStd SHARED IMPORTED)

set_target_properties(TKStd PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKernel;TKCDF;TKCAF;TKLCAF;TKBRep;TKMath;TKG2d;TKG3d;TKStdL"
)

# Create imported target TKTObj
add_library(TKTObj SHARED IMPORTED)

set_target_properties(TKTObj PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKMath;TKLCAF"
)

# Create imported target TKBinTObj
add_library(TKBinTObj SHARED IMPORTED)

set_target_properties(TKBinTObj PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKTObj;TKMath;TKLCAF;TKBinL"
)

# Create imported target TKXmlTObj
add_library(TKXmlTObj SHARED IMPORTED)

set_target_properties(TKXmlTObj PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKCDF;TKernel;TKTObj;TKMath;TKLCAF;TKXmlL"
)

# Create imported target TKVCAF
add_library(TKVCAF SHARED IMPORTED)

set_target_properties(TKVCAF PROPERTIES
  INTERFACE_INCLUDE_DIRECTORIES "${_IMPORT_PREFIX}/include/opencascade"
  INTERFACE_LINK_LIBRARIES "TKernel;TKGeomBase;TKBRep;TKTopAlgo;TKMath;TKService;TKG2d;TKG3d;TKCDF;TKLCAF;TKBO;TKCAF;TKV3d"
)

if(CMAKE_VERSION VERSION_LESS 2.8.12)
  message(FATAL_ERROR "This file relies on consumers using CMake 2.8.12 or greater.")
endif()

# Load information for each installed configuration.
get_filename_component(_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(GLOB CONFIG_FILES "${_DIR}/OpenCASCADEApplicationFrameworkTargets-*.cmake")
foreach(f ${CONFIG_FILES})
  include(${f})
endforeach()

# Cleanup temporary variables.
set(_IMPORT_PREFIX)

# Loop over all imported files and verify that they actually exist
foreach(target ${_IMPORT_CHECK_TARGETS} )
  foreach(file ${_IMPORT_CHECK_FILES_FOR_${target}} )
    if(NOT EXISTS "${file}" )
      message(FATAL_ERROR "The imported target \"${target}\" references the file
   \"${file}\"
but this file does not exist.  Possible reasons include:
* The file was deleted, renamed, or moved to another location.
* An install or uninstall procedure did not complete successfully.
* The installation package was faulty and contained
   \"${CMAKE_CURRENT_LIST_FILE}\"
but not all the files it references.
")
    endif()
  endforeach()
  unset(_IMPORT_CHECK_FILES_FOR_${target})
endforeach()
unset(_IMPORT_CHECK_TARGETS)

# Make sure the targets which have been exported in some other 
# export set exist.
unset(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets)
foreach(_target "TKernel" "TKGeomBase" "TKBRep" "TKTopAlgo" "TKMath" "TKG2d" "TKG3d" "TKBO" "TKService" "TKV3d" )
  if(NOT TARGET "${_target}" )
    set(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets "${${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets} ${_target}")
  endif()
endforeach()

if(DEFINED ${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets)
  if(CMAKE_FIND_PACKAGE_NAME)
    set( ${CMAKE_FIND_PACKAGE_NAME}_FOUND FALSE)
    set( ${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE "The following imported targets are referenced, but are missing: ${${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets}")
  else()
    message(FATAL_ERROR "The following imported targets are referenced, but are missing: ${${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets}")
  endif()
endif()
unset(${CMAKE_FIND_PACKAGE_NAME}_NOT_FOUND_MESSAGE_targets)

# Commands beyond this point should not need to know the version.
set(CMAKE_IMPORT_FILE_VERSION)
cmake_policy(POP)
