
add_library(Model1 Model1.cpp)
target_include_directories(Model1 PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(Model1 MCTDHlib)

