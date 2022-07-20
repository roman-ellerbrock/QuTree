
# GoogleTest requires at least C++14
#set(CMAKE_CXX_STANDARD 14)

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
  Tensor/Tensor_test.cpp
  Tensor/TensorBLAS1_test.cpp
  Tensor/TensorBLAS2_test.cpp
  Tensor/TensorLapack_test.cpp

  Tree/PrimitiveBasis_test.cpp
  Tree/Node_test.cpp
  Tree/Edge_test.cpp
  Tree/Tree_test.cpp
  Tree/sparse_vector_test.cpp
  Tree/SubTree_test.cpp
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


set(QuTree_tests2
        # Tensor tests
        Tensor/test_Tensor.cpp
        Tensor/test_TensorBLAS1.cpp
        Tensor/test_TensorBLAS2.cpp
        Tensor/test_TensorLapack.cpp

        Tensor/test_cuTensor.cpp

        # Tree tests
        Tree/test_Edge.cpp
        Tree/test_Node.cpp
        Tree/test_PrimitiveBasis.cpp
        Tree/test_sparse_vector.cpp
        Tree/test_SubTree.cpp
        Tree/test_Tree.cpp
        Tree/test_TreeFactory.cpp

        # Tensor Network classes
        TensorNetwork/test_ApplySCF.cpp
        TensorNetwork/test_TensorTree.cpp
        TensorNetwork/test_contraction.cpp

        # Tensor Network classes
        Operators/test_LeafOperator.cpp
        Operators/test_ProductOperator.cpp

        # all tests
        tests.cpp

#        test_GradientDescent.cpp
#        test_Integrator.cpp
#        test_MatrixTensorTree.cpp
#        test_TensorTree.cpp
#        test_TensorBLAS.cpp
#        test_Operators.cpp
#        test_SparseTree.cpp
#        test_SymTensorTree.cpp
#        test_TensorTreeBasis.cpp
#        test_TensorTreeOperator.cpp
#        test_MatrixTree.cpp
#        test_SparseMatrixTree.cpp
#        test_RandomMatrices.cpp
        )

#target_link_libraries(TestQuTree UnitTest++)


# run unit tests as post build step
#add_custom_command(TARGET TestQuTree
#    POST_BUILD COMMAND TestQuTree
#    COMMENT "Running unit tests")


#set(QuTree_benchmarks
#        benchmarks/benchmark_tensor.cpp
#        benchmarks/benchmark_tensor.h
#        benchmarks/benchmark_tree.cpp benchmarks/benchmark_tree.h
#        benchmarks/benchmark_helper.cpp benchmarks/benchmark_helper.h
#        benchmarks/optimize_matrixtensor.cpp benchmarks/optimize_matrixtensor.h
#        )
#
#add_executable(benchmarks benchmarks/main_benchmarks.cpp ${QuTree_benchmarks})
#target_link_libraries(benchmarks QuTree /usr/local/opt/openblas/lib/libblas.dylib)
#target_compile_features(benchmarks PRIVATE cxx_std_14)
#
#
