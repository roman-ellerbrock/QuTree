
# GoogleTest requires at least C++14
#set(CMAKE_CXX_STANDARD 14)

    message("PROJECT NAME:         ${CMAKE_PROJECT_DIR}")
include(FetchContent)
FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG release-1.12.1
)

# For Windows: Prevent overriding the parent project's compiler/linker settings
set(gtest_force_shared_crt ON CACHE BOOL "" FORCE)
FetchContent_MakeAvailable(googletest)

enable_testing()

add_executable(
  TestQuTree
  tests.cc

  Tensor/cuTensor_test.cpp
  Tensor/cuTensormxp_test.cpp
  Tensor/cuTensorBLAS1_test.cpp
  Tensor/Tensor_test.cpp
  Tensor/Tensormxp_test.cpp
  Tensor/TensorBLAS1_test.cpp
  Tensor/TensorBLAS2_test.cpp
  Tensor/TensorLapack_test.cpp
  Tensor/ManybodyTensor_test.cpp
#  Tensor/TensorTemplating_test.cpp

  Tree/PrimitiveBasis_test.cpp
  Tree/Node_test.cpp
  Tree/Edge_test.cpp
  Tree/Tree_test.cpp
  Tree/sparse_vector_test.cpp
  Tree/SubTree_test.cpp

  TensorNetwork/TensorTree_test.cpp
  TensorNetwork/ApplySCF_test.cpp

  Operators/LeafOperator_test.cpp
  Operators/ProductOperator_test.cpp

  Util/blaspp_device_extension_test.cpp
)

target_include_directories(TestQuTree 
	PUBLIC
	$<BUILD_INTERFACE:${CMAKE_PROJECT_DIR}/src>
	${BLAS_DIRECTORIES}
	${PolymorphicMemory_DIR}
        ${BLAS_INCLUDE}
	${CUDA_INCLUDE_DIRS}
	)

target_compile_features(TestQuTree INTERFACE cxx_std_17)

if (CUDA_FEATS)
	set_target_properties(TestQuTree PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	target_link_options(TestQuTree INTERFACE ${CUDA_FLAGS})
endif()

target_link_directories(TestQuTree
	PUBLIC	
	${PolymorphicMemory_LINK}
)

target_link_libraries(
        TestQuTree 
	PUBLIC 
        gtest_main
	QuTree 
	${BLAS_LIBRARIES} -lblaspp
	${LAPACK_LIBRARIES} -llapackpp
	polymorphicMemory
)

include(GoogleTest)
gtest_discover_tests(TestQuTree)


add_executable(benchmark_qutree
  benchmarks/benchmark_helper.cpp
  benchmarks/benchmark_tensor.cpp
  benchmarks/main_benchmarks.cpp
)

target_compile_features(benchmark_qutree PRIVATE cxx_std_17)
target_link_libraries(
        benchmark_qutree 
	PUBLIC 
	QuTree 
	${BLAS_LIBRARIES} -lblaspp
	${LAPACK_LIBRARIES} -llapackpp
  OpenMP::OpenMP_CXX 
	polymorphicMemory
)
