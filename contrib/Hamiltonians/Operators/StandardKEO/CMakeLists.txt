
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(StandardKEO StandardKEO.cpp)
target_include_directories (StandardKEO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (StandardKEO MCTDHlib)

