
cmake_minimum_required(VERSION 3.13)
project(QDlib)

set(CMAKE_CXX_STANDARD 14)

# Required packages
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR})
find_package(Eigen3)
if (NOT EIGEN3_FOUND)
    set( EIGEN3_INCLUDE_DIR $ENV{EIGEN3_INCLUDE_DIR})
    if (NOT EIGEN3_INCLUDE_DIR)
        message( FATAL_ERROR "No Eigen3 found, either point CMAKE_PREFIX_PATH or EIGEN3_INCLUDE_DIR to your Eigen installation." )
    else ()
        message ( STATUS "Using EIGEN3_INCLUDE_DIR environment variable for Eigen 3 location instead" )
    endif ()
endif ()
include_directories("${EIGEN3_INCLUDE_DIR}")

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math ")
#set(FLAGS_VECTORIZEE "-ftree-vectorize -mavx -mavx2 -msse2 -mfma -mbmi")
set(FLAGS_VECTORIZEE "-fopenmp-simd -march=native -Rpass='loop|vect' -Rpass-missed='loop|vect' -Rpass-analysis='loop|vect' ")

# Link the libraries "Operators" and "PESs" into a single library "Hamiltonians"
add_library(QDlib
	TensorDim.cpp
	TensorDim_Extension.cpp
	TemplateInstantiations_1.cpp
	TemplateInstantiations_2.cpp	
	FFT.cpp
	JacobiRotationFramework.cpp 
	SimultaneousDiagonalization.cpp
	stdafx.cpp
	CWrapper.f
	long_integer.cpp long_integer.h
	TreeNode.h Tree.h
	QMConstants.h)

target_include_directories(QDlib PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

