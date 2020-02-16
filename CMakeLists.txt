###########################################################################
#
# --- QuantumDynamics library (QDlib) ---
#
# Description:
# QDlib is a linear algebra library suited for quantum dynamics simulations.
# The main intended use of the library are tensor tree applications.
#
# Main author: Roman Ellerbrock
# Contributing authors: Stefan Seritan, K. Grace Johnson,
#                       Thomas Weike, Tim Lenzen
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
  message ( STATUS "Eigen not found, switching to internal submodule" )
  execute_process(COMMAND git submodule update --init -- external/eigen
      WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

  # Set env and install headers
  set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
      CACHE PATH "Eigen include directory")
  install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen DESTINATION include)

  # Convenience target for exports (TODO: Figure out what this actually does)
  add_library(Eigen INTERFACE)
  target_include_directories(Eigen INTERFACE
      $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
      $<INSTALL_INTERFACE:include>)
  install(TARGETS Eigen EXPORT eigen_export DESTINATION include)
endif ()
include_directories("${EIGEN3_INCLUDE_DIR}")

set(openmp OFF)
if (openmp)
    if (APPLE)
        # Assumes you've installed llvm openmp using homebrew (brew install llvm)
        message("On MAC OS - Setting Clan-compiler & OpenMP lib")
        set(CMAKE_C_COMPILER "/usr/local/opt/llvm/bin/clang")
        set(CMAKE_CXX_COMPILER "/usr/local/opt/llvm/bin/clang++")
        set(OPENMP_LIBRARIES "/usr/local/opt/llvm/lib")
        set(OPENMP_INCLUDES "/usr/local/opt/llvm/include")

        OPTION(USE_OpenMP "Use OpenMP to enable <omp.h>" ON)

        # Find OpenMP
        if (USE_OpenMP)
            if (CMAKE_C_COMPILER_ID MATCHES "Clang")
                set(OpenMP_C "${CMAKE_C_COMPILER}")
                set(OpenMP_C_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
                set(OpenMP_C_LIB_NAMES "libomp" "libgomp" "libiomp5")
                set(OpenMP_libomp_LIBRARY ${OpenMP_C_LIB_NAMES})
                set(OpenMP_libgomp_LIBRARY ${OpenMP_C_LIB_NAMES})
                set(OpenMP_libiomp5_LIBRARY ${OpenMP_C_LIB_NAMES})
            endif ()
            if (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
                set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
                set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
                set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
                set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
                set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
                set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})
            endif ()
        endif ()
    endif ()
    find_package(OPENMP REQUIRED)
    if (OPENMP_FOUND)
        include_directories("${OPENMP_INCLUDES}")
        link_directories("${OPENMP_LIBRARIES}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        # set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif (OPENMP_FOUND)
endif ()


#####################################################################
# Set Compiler Flags
#####################################################################
# Ensure C++ 11 standard
if (CMAKE_VERSION VERSION_LESS "3.1")
    if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
        set(CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")
    endif ()
else ()
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED True)
endif ()

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -g")
if (openmp)
    set(CMAKE_CXX_FLAGS_RELEASE "-fopenmp -O3 -ffast-math ")
    set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect' ")
endif()
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
# @TODO: Remove include_directories, instead use target_include_directories

include_directories("include")
include_directories("include/TreeClasses")
include_directories("include/TreeShape/LeafTypes")
include_directories("include/TreeOperators")
include_directories("tests/benchmarks")
add_subdirectory("src")
add_subdirectory("include")
add_subdirectory("tests")
add_subdirectory("examples")

