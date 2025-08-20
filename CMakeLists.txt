
cmake_minimum_required(VERSION 3.11)

project (egads LANGUAGES CXX C)

add_library(egads INTERFACE)

target_include_directories(egads SYSTEM
  INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/include"
)

target_link_directories(egads
  INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/lib"
)
target_link_libraries(egads INTERFACE libegads.so)
set_property(TARGET egads PROPERTY IMPORTED_LOCATION "${CMAKE_CURRENT_SOURCE_DIR}/lib/libegads.so")
set_property(TARGET egads PROPERTY IMPORTED_NO_SONAME TRUE) # strip path from IMPORTED_LOCATION
set_property(TARGET egads PROPERTY CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
set_property(TARGET egads PROPERTY INSTALL_RPATH "$ORIGIN:$ORIGIN/../lib:${CMAKE_CURRENT_SOURCE_DIR}/lib")

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/lib/ DESTINATION lib)

add_subdirectory(test)
