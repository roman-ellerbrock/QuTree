add_executable(example_1 example_1.cpp)
target_link_libraries(example_1 QuTree)

add_executable(example_2 example_2.cpp)
target_link_libraries(example_2 QuTree::QuTree)

add_executable(example_3 example_3.cpp)
target_link_libraries(example_3 QuTree)

# run unit tests as post build step
add_custom_target(TestExamples
    COMMAND "./example_1"
    COMMAND "./example_2"
    COMMAND "./example_3"
    DEPENDS example_1 example_2 example_3
    COMMENT "Run all examples")
target_link_libraries(example_3 QuTree::QuTree)
