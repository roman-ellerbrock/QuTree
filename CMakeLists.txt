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
# Target setup
#####################################################################
set(QDlib_INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/QDlib)
set(QDlib_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(QDlib_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/QDlib)

include(src/CMakeLists.txt)
include(include/CMakeLists.txt)

add_library(QDlib SHARED ${QDlib_SOURCE_FILES})

# Set up headers in include/ for install & build and src/ for just building
target_include_directories(QDlib
    PUBLIC
        $<INSTALL_INTERFACE:${QDlib_INSTALL_INC_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )
set_target_properties(QDlib PROPERTIES PUBLIC_HEADER "${QDlib_INCLUDE_FILES}")

# Install rules
install(TARGETS QDlib
    EXPORT QDlib-export
    PUBLIC_HEADER DESTINATION ${QDlib_INSTALL_INC_DIR}
    LIBRARY DESTINATION ${QDlib_INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${QDlib_INSTALL_LIB_DIR}
    )

install(EXPORT QDlib-export
    FILE QDlibTargets.cmake
    NAMESPACE QDlib::
    DESTINATION ${QDlib_INSTALL_CMAKE_DIR}
    )

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
# For debugging
#set(CMAKE_VERBOSE_MAKEFILE ON)

# Ensure C++14 standard
target_compile_features(QDlib PRIVATE cxx_std_14)

#set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect'")
set(QDlib_DEBUG_FLAGS "-g")
set(QDlib_RELEASE_FLAGS "-O3 -ffast-math")

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

# Get command-line ready options for target_compile: https://stackoverflow.com/a/27651464/3052876
separate_arguments(QDlib_DEBUG_FLAGS UNIX_COMMAND "${QDlib_DEBUG_FLAGS}")
separate_arguments(QDlib_RELEASE_FLAGS UNIX_COMMAND "${QDlib_RELEASE_FLAGS}")

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

    # Set env and install headers as subdir of QDlib to avoid clashes
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
        CACHE PATH "Eigen include directory")
    install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen
        DESTINATION ${QDlib_INSTALL_INC_DIR})

    # Convenience target for exports
    add_library(Eigen INTERFACE)
    target_include_directories(Eigen INTERFACE
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:${QDlib_INSTALL_INC_DIR}>)
    install(TARGETS Eigen
        EXPORT eigen_export
        DESTINATION ${QDlib_INSTALL_INC_DIR})
endif ()
include_directories("${EIGEN3_INCLUDE_DIR}")

#####################################################################
# Recurse into subdirectories
#####################################################################
add_subdirectory("tests")
add_subdirectory("examples")

