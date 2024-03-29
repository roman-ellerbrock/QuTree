###########################################################################
#
# --- QuTree library ---
#
# Description:
# QuTree is a linear algebra library suited for tensor tree applications.
# The library was originally designed for quantum dynamics simulations.
#
# Main author: Roman Ellerbrock
# Contributing authors: K. Grace Johnson, Stefan Seritan,
#                       H. J. Zhang, H. Hoppe
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
set(CMAKE_VERBOSE_MAKEFILE ON) # For debugging

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

#####################################################################
# Find Required packages
#####################################################################
option(USE_INTERNAL_EIGEN "Use internal Eigen module" OFF)
if (NOT USE_INTERNAL_EIGEN)
    find_package(Eigen3 NO_MODULE)
endif ()
if (USE_INTERNAL_EIGEN OR NOT TARGET Eigen3::Eigen)
    message(STATUS "Eigen not found, switching to internal submodule")
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
endif ()

#####################################################################
# Target Installation setup
#####################################################################
# Set up BLAS / LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
if (APPLE)
   if (CMAKE_HOST_SYSTEM_PROCESSOR STREQUAL arm64)
      message(STATUS "USING MAC OS with M1!")
      set(BLAS_LIBRARIES "-L/opt/homebrew/Cellar/openblas/0.3.21/lib -lblas")
      target_include_directories(QuTree PUBLIC "/opt/homebrew/Cellar/openblas/0.3.21/include")
      include_directories("/opt/homebrew/Cellar/openblas/0.3.21/include")
      set(LAPACK_LIBRARIES "-L/opt/homebrew/Cellar/lapack/3.11/lib -llapacke -llapack")
      target_include_directories(QuTree PUBLIC "/opt/homebrew/Cellar/lapack/3.11/include")
      include_directories("/opt/homebrew/Cellar/lapack/3.11/include")
   else ()
      message(STATUS "USING MAC OS!")
      set(BLAS_LIBRARIES "/usr/local/opt/openblas/lib/libblas.dylib")
      target_include_directories(QuTree PUBLIC "/usr/local/opt/openblas/include")
      include_directories("/usr/local/opt/openblas/include")
   endif ()	
elseif (UNIX)
    #    set(BLA_VENDOR OpenBLAS)
    message(STATUS "USING UNIX!")

    #   empirically found blas locations
    target_include_directories(QuTree PUBLIC "/usr/include/x86_64-linux-gnu")
    include_directories("/usr/include/x86_64-linux-gnu")
endif ()

target_link_libraries(QuTree PUBLIC lapacke blas lapack)

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

if (NOT openmp)
    set(openmp OFF)
    message("Setting openMP off (${openmp})")
else ()
    message("openmp=${openmp}")
endif ()

if (openmp)
    find_package(OpenMP REQUIRED)
    target_link_libraries(QuTree PRIVATE OpenMP::OpenMP_CXX OpenMP::OpenMP_Fortran)
    set(QuTree_RELEASE_FLAGS "${OpenMP_CXX_FLAGS} ${QuTree_RELEASE_FLAGS}")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${OpenMP_Fortran_FLAGS}")
endif ()

target_link_libraries(QuTree PUBLIC Eigen3::Eigen ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})

#####################################################################
# Set Compiler Flags
#####################################################################

# Ensure C++14 standard (force downstream as well)
target_compile_features(QuTree PUBLIC cxx_std_17)

set(FLAGS_VECTORIZE "-fopenmp-simd")
set(QuTree_RELEASE_FLAGS "${FLAGS_VECTORIZE} -O3 -ffast-math")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast -ffast-math")

#include(vectorization.cmake)

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
add_subdirectory(contrib EXCLUDE_FROM_ALL)
