
add_library(NOCl_PES NOClPotential.cpp)
target_include_directories(NOCl_PES PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(NOCl_PES MCTDHlib)

