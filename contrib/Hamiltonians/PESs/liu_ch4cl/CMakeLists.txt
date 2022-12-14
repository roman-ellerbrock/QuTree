add_library(liu_ch4cl pipNN/potential.f pipNN/ch4cl-pipnn.f liuch4cl.cpp)
target_include_directories(liu_ch4cl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(liu_ch4cl PESs)