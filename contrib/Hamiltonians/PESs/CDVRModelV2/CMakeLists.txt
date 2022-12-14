
add_library(CDVRModelV2 CDVRModelV2.cpp)
target_include_directories(CDVRModelV2 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(CDVRModelV2 MCTDHlib)

