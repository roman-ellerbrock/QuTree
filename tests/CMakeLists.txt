find_package(UnitTest++ QUIET)
if (NOT UnitTest++_FOUND)
  message ( STATUS "UnitTest++ not found, switching to internal submodule" )
  execute_process(COMMAND git submodule update --init -- external/unittest-cpp
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory("external/unittest-cpp")
  include_directories("external/unittest-cpp")
endif()
include_directories(${UTPP_INCLUDE_DIRS})

set(QuTree_tests
        test_GradientDescent.cpp
        test_Integrator.cpp
        test_Matrix.cpp
        test_Tensor.cpp
        tests.cpp
        test_TensorTree.cpp
        test_Operators.cpp
        test_TensorTreeBasis.cpp
        test_MatrixTree.cpp
        test_SparseMatrixTree.cpp)

add_executable(TestQuTree ${QuTree_tests})
target_link_libraries(TestQuTree QuTree)
target_link_libraries(TestQuTree UnitTest++)
target_compile_features(TestQuTree PRIVATE cxx_std_14)

# run unit tests as post build step
add_custom_command(TARGET TestQuTree
    POST_BUILD COMMAND TestQuTree
    COMMENT "Running unit tests")


set(QuTree_benchmarks
        benchmarks/benchmark_tensor.cpp
        benchmarks/benchmark_tensor.h
        benchmarks/benchmark_tree.cpp benchmarks/benchmark_tree.h
        benchmarks/benchmark_helper.cpp benchmarks/benchmark_helper.h)

add_executable(benchmarks benchmarks/main_benchmarks.cpp ${QuTree_benchmarks})
target_link_libraries(benchmarks QuTree)
target_compile_features(benchmarks PRIVATE cxx_std_14)
