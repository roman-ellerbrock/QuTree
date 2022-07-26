



#####################################################################
# CUDA support
#####################################################################

option (CUDA_FEATS "Set to On to use CUDA features" ON)

if (CUDA_FEATS)
	message("CUDA compilation enabled.")

	enable_language(CUDA)
	find_package(CUDA REQUIRED)

	# Set CUDA_ARCHITECTURES
#	cmake_policy(SET "CMP0104" NEW)
	if(NOT DEFINED ${CMAKE_CUDA_ARCHITECTURES})
		set(CMAKE_CUDA_ARCHITECTURES 52 61 75)
	endif()

	set_target_properties(QuTree PROPERTIES 
		CUDA_SEPARABLE_COMPILATION ON
	)

	message("CUDA EXEC: ${CUDA_NVCC_EXECUTABLE}")
	message("CUDA COMP: ${CMAKE_CUDA_COMPILER}")
	set_target_properties(QuTree PROPERTIES 
		IMPORTED_LOCATION "/usr/local/cuda/lib64/"
		INTERFACE_INCLUDE_DIRECTORIES "/usr/local/cuda/include"
		CUDA_ARCHITECTURES ${CMAKE_CUDA_ARCHITECTURES}
	)

	target_link_options(QuTree PUBLIC ${CUDA_FLAGS} -lcublas)

#	find_library(CUBLAS
#	NAMES cublas
#	HINTS /usr/local/cuda/lib64
#	OPTIONAL)
#
#	find_library(CURAND
#	NAMES curand
#	HINTS /usr/local/cuda/lib64
#	OPTIONAL)
#	set(CUDA_FLAGS "${CUBLAS} ${CURAND}")

	#	separate_arguments(CUDA_FLAGS UNIX_COMMAND "${CUDA_FLAGS}")

	message("CUDA Flags: ${CUDA_FLAGS}")
	message("CUDA Include directories: ${CUDA_INCLUDE_DIRS}")

endif()




