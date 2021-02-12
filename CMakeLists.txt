###########################################################################
#
# --- QuTree library ---
#
# Description:
# QuTree is a linear algebra library suited for tensor tree applications.
# The library was originally designed for quantum dynamics simulations.
#
# Main author: Roman Ellerbrock
# Contributing authors: Stefan Seritan, K. Grace Johnson,
#                       Thomas Weike, Tim Lenzen
###########################################################################

cmake_minimum_required(VERSION 3.0)
cmake_policy(SET "CMP0048" NEW)
project(QuTree LANGUAGES C CXX VERSION 0.1.0)
enable_language(C)
enable_language(CXX)
enable_language(Fortran)
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
cmake_policy(SET "CMP0042" NEW)

#####################################################################
# Target setup
#####################################################################
set(QuTree_INSTALL_INC_DIR ${CMAKE_INSTALL_PREFIX}/include/QuTree)
set(QuTree_INSTALL_LIB_DIR ${CMAKE_INSTALL_PREFIX}/lib)
set(QuTree_INSTALL_CMAKE_DIR ${CMAKE_INSTALL_PREFIX}/lib/cmake/QuTree)

include(src/CMakeLists.txt)
include(include/CMakeLists.txt)

add_library(QuTree SHARED ${QuTree_SOURCE_FILES})
add_library(QuTree::QuTree ALIAS QuTree)

#find_package(Boost REQUIRED)
#include_directories(${Boost_INCLUDE_DIRS})
#####################################################################
# Find Required packages
#####################################################################
option(USE_INTERNAL_EIGEN "Use internal Eigen module" OFF)
if (NOT USE_INTERNAL_EIGEN)
    find_package(Eigen3 NO_MODULE)
endif()
if (USE_INTERNAL_EIGEN OR NOT TARGET Eigen3::Eigen)
    message ( STATUS "Eigen not found, switching to internal submodule" )
    set(USE_INTERNAL_EIGEN ON)

    execute_process(COMMAND git submodule update --init -- external/eigen
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

    # Set env and install headers as subdir of QuTree to avoid clashes
    set(EIGEN3_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/external/eigen
        CACHE PATH "Eigen include directory")
    install(DIRECTORY ${EIGEN3_INCLUDE_DIR}/Eigen
        DESTINATION ${QuTree_INSTALL_INC_DIR})

    # Convenience target for exports
    add_library(Eigen INTERFACE)
    add_library(Eigen3::Eigen ALIAS Eigen)
    target_include_directories(Eigen INTERFACE
        $<BUILD_INTERFACE:${EIGEN3_INCLUDE_DIR}>
        $<INSTALL_INTERFACE:${QuTree_INSTALL_INC_DIR}>)
    install(TARGETS Eigen
        EXPORT QuTree-export
        DESTINATION ${QuTree_INSTALL_INC_DIR})
endif()
target_link_libraries(QuTree Eigen3::Eigen)

#####################################################################
# Target Installation setup
#####################################################################
# Set up headers in include/ for install & build and src/ for just building
target_include_directories(QuTree
    PUBLIC
        $<INSTALL_INTERFACE:${QuTree_INSTALL_INC_DIR}>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
    )

# Install rules
install(TARGETS QuTree
    EXPORT QuTree-export
    LIBRARY DESTINATION ${QuTree_INSTALL_LIB_DIR}
    ARCHIVE DESTINATION ${QuTree_INSTALL_LIB_DIR}
    )
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/
    DESTINATION ${QuTree_INSTALL_INC_DIR}
    FILES_MATCHING PATTERN "*.h*")

install(EXPORT QuTree-export
    FILE QuTreeTargets.cmake
    NAMESPACE QuTree::
    DESTINATION ${QuTree_INSTALL_CMAKE_DIR})

# Exports for external find_package
include(CMakePackageConfigHelpers)
configure_file(QuTreeConfig.cmake.in QuTreeConfig.cmake @ONLY)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfigVersion.cmake
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY SameMajorVersion
    )

install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/QuTreeConfigVersion.cmake
    DESTINATION ${QuTree_INSTALL_CMAKE_DIR})

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
    find_package(OpenMP REQUIRED)
#    find_package(OPENMP REQUIRED)
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
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Ensure C++14 standard (force downstream as well)
target_compile_features(QuTree PUBLIC cxx_std_14)

#set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect'")
set(FLAGS_VECTORIZE "-fopenmp-simd -march=native -loop-vectorize")
#set(FLAGS_VECTORIZE "-ftree-parallelize-loops=8 ${FLAGS_VECTORIZE}")
set(QuTree_RELEASE_FLAGS "${FLAGS_VECTORIZE} -O3 -ffast-math")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -openmp -Ofast -ffast-math")

#include(vectorization.cmake)

set(QuTree_DEBUG_FLAGS "-g")
#set(QuTree_RELEASE_FLAGS " -O3 -ffast-math")

if (APPLE)
    message("On OSX")
    set(QuTree_DEBUG_FLAGS "-mmacosx-version-min=10.14 ${QuTree_DEBUG_FLAGS}")
    set(QuTree_RELEASE_FLAGS "-mmacosx-version-min=10.14 ${QuTree_RELEASE_FLAGS}")
    set(CMAKE_MACOSX_RPATH 1)
elseif (UNIX)
    message("On UNIX")
else ()
    message("Other OS")
endif ()

# Get command-line ready options for target_compile: https://stackoverflow.com/a/27651464/3052876
separate_arguments(QuTree_DEBUG_FLAGS UNIX_COMMAND "${QuTree_DEBUG_FLAGS}")
separate_arguments(QuTree_RELEASE_FLAGS UNIX_COMMAND "${QuTree_RELEASE_FLAGS}")

# Set debug/release flags: https://stackoverflow.com/a/23995391/3052876
target_compile_options(QuTree PUBLIC "$<$<CONFIG:DEBUG>:${QuTree_DEBUG_FLAGS}>")
target_compile_options(QuTree PRIVATE "$<$<CONFIG:RELEASE>:${QuTree_RELEASE_FLAGS}>")

if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif ()


#####################################################################
# Recurse into subdirectories
#####################################################################
# Excluded from default builds
add_subdirectory(tests EXCLUDE_FROM_ALL)
add_subdirectory(examples EXCLUDE_FROM_ALL)

