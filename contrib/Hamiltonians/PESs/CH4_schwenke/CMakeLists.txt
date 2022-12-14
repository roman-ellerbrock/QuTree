
add_library(PES_CH4 PES_CH4.cpp PES_CH4.f)
target_include_directories(PES_CH4 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_CH4 MCTDHlib)

