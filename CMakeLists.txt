
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

