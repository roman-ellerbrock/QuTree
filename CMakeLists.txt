###########################################################################
#
# --- QuantumDynamics library (QDlib) ---
#
# Description:
# QDlib is a linear algebra library suited for quantum dynamics simulations.
# The main intended use of the library are tensor tree applications.
#
# Main author: Roman Ellerbrock
# Contributing authors: Thomas Weike, Tim Lenzen, Stefan Seritan, K. Grace Johnson
###########################################################################

cmake_minimum_required(VERSION 2.8)
project(QDlib C CXX)
enable_language(C)
enable_language(CXX)
set(QDlib_VERSION_MAJOR 0)
set(QDlib_VERSION_MINOR 1)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
cmake_policy(SET "CMP0042" NEW)

#####################################################################
# Find Required packages
#####################################################################
find_package(Eigen3)
if (NOT EIGEN3_FOUND)
  set( EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
  if (NOT EIGEN3_INCLUDE_DIR)
    message( FATAL_ERROR "No Eigen3 found, either point CMAKE_PREFIX_PATH or EIGEN3_INCLUDE_DIR to your Eigen installation." )
  else ()
    message ( STATUS "Using EIGEN3_INCLUDE_DIR environment variable for Eigen 3 location instead: ${EIGEN3_INCLUDE_DIR}" )
  endif ()
endif ()
include_directories("${EIGEN3_INCLUDE_DIR}")

#####################################################################
# Set Compiler Flags
#####################################################################
# Ensure C++ 11 standard
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
else ()
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif ()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math ")
set(FLAGS_VECTORIZEE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect' ")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${FLAGS_VECTORIZE}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} -mcmodel=large")
set(FFLAGS "-heap-arrays -mcmodel=large")
if (APPLE)
    message("On OSX")
    set(CMAKE_CXX_FLAGS_RELEASE "-mmacosx-version-min=10.14 ${CMAKE_CXX_FLAGS_RELEASE}")
    set(CMAKE_MACOSX_RPATH 1)
elseif (UNIX)
    message("On UNIX")
    set(CMAKE_CXX_FLAGS_RELEASE "-ftree-parallelize-loops=8  ${CMAKE_CXX_FLAGS_RELEASE}")
else ()
    message("Other OS")
    set(CMAKE_CXX_FLAGS_RELEASE " ${CMAKE_CXX_FLAGS_RELEASE}")
endif ()


#####################################################################
# Set Build Options
#####################################################################
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

#####################################################################
# Recurse into subdirectories
#####################################################################
include_directories("include")
include_directories("include/Tree")
include_directories("include/TensorTreeBasis/LeafTypes")
include_directories("include/Operators")
add_subdirectory("src")
add_subdirectory("include")
add_subdirectory("tests")

#####################################################################
# Install QDlib
#####################################################################
install (TARGETS QDlib
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
    RUNTIME DESTINATION bin)
