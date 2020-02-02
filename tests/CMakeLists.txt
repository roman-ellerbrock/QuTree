
include_directories("include")
find_package(UnitTest++)
if (NOT UnitTest++_FOUND)
  message ( STATUS "UnitTest++ not found, switching to internal submodule" )
  execute_process(COMMAND git submodule update --init -- external/unittest-cpp
	  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
  add_subdirectory("external/unittest-cpp")
  include_directories("external/unittest-cpp")
endif()
include_directories(${UTPP_INCLUDE_DIRS})

set(QDlib_tests
        test_Matrix.cpp
        test_Tensor.cpp
        tests.cpp
        test_Tree.cpp
        test_Operators.cpp
        test_.FactorMatrixTree.cpp
        test_SparseFactorMatrixTree.cpp
        test_Core.cpp)

add_executable(TestQDlib ${QDlib_tests})
target_link_libraries(TestQDlib QDlib)
target_link_libraries(TestQDlib UnitTest++)

# run unit tests as post build step
add_custom_command(TARGET TestQDlib
    POST_BUILD COMMAND TestQDlib
    COMMENT "Running unit tests")
