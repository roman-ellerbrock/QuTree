
add_library(PES_CH5P PES_CH5P.cpp PES_CH5P.f)
target_include_directories(PES_CH5P PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_CH5P MCTDHlib)

