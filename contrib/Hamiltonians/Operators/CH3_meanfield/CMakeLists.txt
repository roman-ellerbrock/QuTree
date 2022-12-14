
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH3_meanfield CH3_meanfield.cpp)
target_include_directories (CH3_meanfield PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (CH3_meanfield Operators)

