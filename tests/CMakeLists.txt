include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        release-1.11.0
)

FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

add_executable(
    tests
    tests.cpp
    test_Graph.cpp
    test_GraphFactory.cpp
    test_TensorNetwork.cpp
    graph_util.cpp
)

target_include_directories(tests PUBLIC
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/tests
)

target_link_libraries(tests PRIVATE
    qutreelib
    GTest::GTest
)

set_target_properties(tests PROPERTIES 
    CXX_STANDARD 17
)

gtest_discover_tests(tests)

add_test(NAME AllTests
    COMMAND tests
)
