
add_library(PES_CH3_1 CH3Potential.cpp CH3Potential.f)
target_include_directories(PES_CH3_1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_CH3_1 PESs)

