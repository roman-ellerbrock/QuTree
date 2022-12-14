
add_library(PES_CH4SM PES_CH4SM.cpp PES_CH4SM.f90)
target_include_directories(PES_CH4SM PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(PES_CH4SM MCTDHlib)

