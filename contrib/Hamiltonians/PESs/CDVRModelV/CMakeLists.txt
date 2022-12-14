
add_library(CDVRModelV CDVRModelV.cpp)
target_include_directories(CDVRModelV PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/src)
target_compile_features(CDVRModelV PRIVATE cxx_std_14)
target_link_libraries(CDVRModelV MCTDHlib)

