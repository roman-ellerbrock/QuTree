
include_directories("include")
find_package(UnitTest++)
include_directories(${UTPP_INCLUDE_DIRS})

set(QDlib_tests
        Tensor_test_1.cpp
        Tensor_test_1.h
        tests.cpp
        )

add_executable(tests ${QDlib_tests})
target_link_libraries(tests QDlib)
target_link_libraries(tests UnitTest++)

