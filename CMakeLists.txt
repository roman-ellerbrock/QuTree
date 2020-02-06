###########################################################################
#
# --- QuantumDynamics library (QDlib) ---
#
# Description:
# QDlib is a linear algebra library suited for quantum dynamics simulations.
# The main intended use of the library are tensor tree applications.
#
# Main author: Roman Ellerbrock
# Contributing authors: Thomas Weike, Tim Lenzen, Stefan Seritan,
#                       K. Grace Johnson
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
# Target setup
#####################################################################
set(QDlib_INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/QDlib)
set(QDlib_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(QDlib_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/QDlib)

include(src/CMakeLists.txt)
include(include/CMakeLists.txt)

add_library(QDlib ${QDlib_SOURCE_FILES})

# Set up headers in include/ for install & build and src/ for just building
target_include_directories(QDlib
    PUBLIC
        $<INSTALL_INTERFACE:${QDlib_INSTALL_INC_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

# Install rules
install(TARGETS QDlib
    EXPORT QDlib-export
    LIBRARY DESTINATION ${QDlib_INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${QDlib_INSTALL_LIB_DIR}
    )

install(EXPORT QDlib-export
    FILE QDlibTargets.cmake
    NAMESPACE QDlib::
    DESTINATION ${QDlib_INSTALL_CMAKE_DIR}
    )

#####################################################################
# Set Compiler Flags
#####################################################################
# Ensure C++14 standard
target_compile_features(QDlib PRIVATE cxx_std_14)

set(QDlib_DEBUG_FLAGS " ")
set(QDlib_RELEASE_FLAGS "-O3 -ffast-math")
set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect'")
set(QDlib_RELEASE_FLAGS "${QDlib_RELEASE_FLAGS} ${FLAGS_VECTORIZE}")
if (APPLE)
    message("On OSX")
    set(QDlib_DEBUG_FLAGS "-mmacosx-version-min=10.14 ${QDlib_DEBUG_FLAGS}")
    set(QDlib_RELEASE_FLAGS "-mmacosx-version-min=10.14 ${QDlib_RELEASE_FLAGS}")
    set(CMAKE_MACOSX_RPATH 1)
elseif (UNIX)
    message("On UNIX")
    set(QDlib_RELEASE_FLAGS "-ftree-parallelize-loops=8  ${QDlib_RELEASE_FLAGS}")
else ()
    message("Other OS")
endif ()

# Set debug/release flags: https://stackoverflow.com/a/23995391/3052876
target_compile_options(QDlib PUBLIC "$<$<CONFIG:DEBUG>:${QDlib_DEBUG_FLAGS}>")
target_compile_options(QDlib PUBLIC "$<$<CONFIG:RELEASE>:${QDlib_RELEASE_FLAGS}>")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()

#####################################################################
# Find Required packages
#####################################################################
find_package(Eigen3 QUIET)
if (NOT EIGEN3_FOUND)
    message ( STATUS "Eigen not found, switching to internal submodule" )
    execute_process(COMMAND git submodule update --init -- external/eigen
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # Set env and install headers
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
        CACHE PATH "Eigen include directory")
    install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen DESTINATION include)

    # Convenience target for exports
    add_library(Eigen INTERFACE)
    target_include_directories(Eigen INTERFACE
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:include>)
    install(TARGETS Eigen EXPORT eigen_export DESTINATION include)
endif ()
include_directories("${EIGEN3_INCLUDE_DIR}")

#####################################################################
# Recurse into subdirectories
#####################################################################
#include_directories("tests/benchmarks")
#add_subdirectory("tests")
#add_subdirectory("examples")

