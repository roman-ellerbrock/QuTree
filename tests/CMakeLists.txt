


find_package(UnitTest++ QUIET)

if (NOT UnitTest++_FOUND)
  message ( STATUS "UnitTest++ not found, falling back to PkgConfig" )
  include(FindPkgConfig)
  pkg_check_modules(UnitTest++ QUIET UnitTest++)
endif()

if (NOT UnitTest++_FOUND)
  message ( STATUS "UnitTest++ not found, switching to internal submodule" )
  execute_process(COMMAND git submodule update --init -- tests/external/unittest-cpp
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../)
  add_subdirectory("external/unittest-cpp")
  include_directories("external/unittest-cpp")
endif()
include_directories(${UTPP_INCLUDE_DIRS})

set(QuTree_tests
        # Tensor tests
        Tensor/test_Tensor.cpp
        Tensor/test_PolymorphicMemory.cpp
	Tensor/test_PolymorphicMemory.cu
	#	Tensor/test_PolymorphicMemory_CUDA.cpp
        Tensor/test_TensorBLAS1.cpp
        Tensor/test_TensorBLAS2.cpp
        Tensor/test_TensorLapack.cpp

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

add_executable(TestQuTree ${QuTree_tests})

#target_link_libraries(TestQuTree QuTree)
#target_link_libraries(TestQuTree UnitTest++)

target_include_directories(TestQuTree 
	PUBLIC
	$<BUILD_INTERFACE:${CMAKE_PROJECT_DIR}/src>
        ${BLAS_INCLUDE}
	${CUDA_INCLUDE_DIRS}
	)

target_compile_features(TestQuTree INTERFACE cxx_std_17)

if (CUDA_FEATS)
	set_target_properties(TestQuTree PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
	target_link_options(TestQuTree INTERFACE ${CUDA_FLAGS})
endif()

target_link_libraries(TestQuTree 
	PUBLIC 
	QuTree 
	UnitTest++ 
	${BLAS_LIBRARIES} -lblaspp
	${LAPACK_LIBRARIES} -llapackpp
)
#target_link_libraries(TestQuTree UnitTest++)

target_include_directories(TestQuTree
	PUBLIC 
	${BLAS_DIRECTORIES}
)

#target_link_libraries(TestQuTree PUBLIC QuTree UnitTest++)
##target_link_libraries(TestQuTree PUBLIC QuTree UnitTest++ ${CUDA_FLAGS})

# run unit tests as post build step
add_custom_command(TARGET TestQuTree
    POST_BUILD COMMAND TestQuTree
    COMMENT "Running unit tests")


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
