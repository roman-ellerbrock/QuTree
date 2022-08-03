###########################################################################
# 
#  ________   ___  ___   _________   ________   _______    _______      
# |\   __  \ |\  \|\  \ |\___   ___\|\   __  \ |\  ___ \  |\  ___ \     
# \ \  \|\  \\ \  \\\  \\|___ \  \_|\ \  \|\  \\ \   __/| \ \   __/|    
#  \ \  \\\  \\ \  \\\  \    \ \  \  \ \   _  _\\ \  \_|/__\ \  \_|/__  
#   \ \  \\\  \\ \  \\\  \    \ \  \  \ \  \\  \|\ \  \_|\ \\ \  \_|\ \ 
#    \ \_____  \\ \_______\    \ \__\  \ \__\\ _\ \ \_______\\ \_______\
#     \|___| \__\\|_______|     \|__|   \|__|\|__| \|_______| \|_______|
#           \|__|                                                       
#                                                                       
# Description:
# QuTree is a linear algebra library suited for tensor tree applications.
# The library was originally designed for quantum dynamics simulations.
#
# Main author: Roman Ellerbrock
# Contributing authors: Stefan Seritan, K. Grace Johnson,
#                       Thomas Weike, Tim Lenzen
###########################################################################

#set(CMAKE_VERBOSE_MAKEFILE ON) # For debugging

cmake_minimum_required(VERSION 3.0)
project(QuTree LANGUAGES C CXX VERSION 0.1.0)
cmake_policy(SET "CMP0048" NEW)
cmake_policy(SET "CMP0042" NEW)
enable_language(C)
enable_language(CXX)
enable_language(Fortran)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})

include(src/CMakeLists.txt)

add_library(QuTree SHARED ${QuTree_SOURCE_FILES} ${QuTree_SOURCE_FILES_CUDA})
add_library(QuTree::QuTree ALIAS QuTree)

# Build-dependent Compiler flags
include(cmake/CompilerFlags.cmake)

# Build polymorphic memory
include(cmake/PolymorphicMemory.cmake)

# optional cuda support
include(cmake/cuda.cmake)

# find BLAS++/LAPACK++ backends
#include(cmake/FindBLASPP.cmake)

include(cmake/VectorizationSimple.cmake)

include(cmake/OpenMP.cmake)

# perform the linking
target_link_libraries(QuTree PUBLIC ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
#target_link_libraries(QuTree PUBLIC Eigen3::Eigen ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
target_include_directories(QuTree 
	PUBLIC 
	${BLAS_DIRECTORIES}
	${PolymorphicMemory_DIR}
)

# Install instructions
include(cmake/Install.cmake)

# Unit tests
add_subdirectory(tests EXCLUDE_FROM_ALL)

# Benchmarking suite
add_subdirectory(benchmark EXCLUDE_FROM_ALL)

# Examples of how to use QuTree
add_subdirectory(examples EXCLUDE_FROM_ALL)

