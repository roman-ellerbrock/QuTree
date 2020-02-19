add_executable(example_1 example_1.cpp)
target_link_libraries(example_1 QDlib)
target_compile_features(example_1 PRIVATE cxx_std_14)

add_executable(example_2 example_2.cpp)
target_link_libraries(example_2 QDlib)
target_compile_features(example_2 PRIVATE cxx_std_14)

add_executable(example_3 example_3.cpp)
target_link_libraries(example_3 QDlib)
target_compile_features(example_3 PRIVATE cxx_std_14)
