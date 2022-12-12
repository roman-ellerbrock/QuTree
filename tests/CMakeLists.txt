include(FetchContent)
FetchContent_Declare(
        googletest
        URL https://github.com/google/googletest/archive/refs/heads/main.zip
)
FetchContent_MakeAvailable(googletest)

enable_testing()

set(QuTree_tests
        test_GradientDescent.cpp
        test_Integrator.cpp
        test_Matrix.cpp
        test_MatrixBLAS.cpp
        test_MatrixTensorTree.cpp
        test_Tensor.cpp
        tests.cpp
        test_TensorTree.cpp
        test_TensorBLAS.cpp
        test_Operators.cpp
        test_SparseTree.cpp
        test_SymTensorTree.cpp
        test_TensorTreeBasis.cpp
#        test_TensorTreeOperator.cpp
        test_MatrixTree.cpp
        test_SparseMatrixTree.cpp
        test_RandomMatrices.cpp
#        test_Node.cpp
        )
add_executable(TestQuTree ${QuTree_tests})
target_link_libraries(TestQuTree QuTree)
target_link_libraries(TestQuTree gtest_main)
target_compile_features(TestQuTree PRIVATE cxx_std_14)

set(QuTree_benchmarks
        benchmarks/benchmark_tensor.cpp
        benchmarks/benchmark_tensor.h
        benchmarks/benchmark_tree.cpp benchmarks/benchmark_tree.h
        benchmarks/benchmark_helper.cpp benchmarks/benchmark_helper.h
        benchmarks/optimize_matrixtensor.cpp benchmarks/optimize_matrixtensor.h
        )

add_executable(benchmarks benchmarks/main_benchmarks.cpp ${QuTree_benchmarks})
target_link_libraries(benchmarks QuTree /usr/local/opt/openblas/lib/libblas.dylib)
target_compile_features(benchmarks PRIVATE cxx_std_14)

