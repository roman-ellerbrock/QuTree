include_directories("include")

add_executable(example_1 example_1.cpp)
target_link_libraries(example_1 QDlib)

add_executable(example_2 example_2.cpp)
target_link_libraries(example_2 QDlib)

add_executable(example_3 example_3.cpp)
target_link_libraries(example_3 QDlib)

add_executable(supervised_learning supervised_learning.cpp)
target_link_libraries(supervised_learning QDlib)