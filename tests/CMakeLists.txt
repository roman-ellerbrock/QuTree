
include_directories("include")
find_package(UnitTest++)
include_directories(${UTPP_INCLUDE_DIRS})

set(QDlib_tests
        test_Vector.cpp
        test_Matrix.cpp
        test_Tensor.cpp
        tests.cpp
        test_Tree.cpp
        test_Operators.cpp)

add_executable(tests ${QDlib_tests})
target_link_libraries(tests QDlib)
target_link_libraries(tests UnitTest++)

