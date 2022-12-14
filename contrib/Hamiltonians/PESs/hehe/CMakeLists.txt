
add_library(hehe_PES hehepotential.cpp hehepotential.f hehepotentials.f90)
target_include_directories(hehe_PES PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(hehe_PES MCTDHlib)

