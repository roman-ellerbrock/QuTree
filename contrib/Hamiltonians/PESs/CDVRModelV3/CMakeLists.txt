
add_library(CDVRModelV3 CDVRModelV.cpp)
target_include_directories(CDVRModelV3 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(CDVRModelV3 MCTDHlib)

