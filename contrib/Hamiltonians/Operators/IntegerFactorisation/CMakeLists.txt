#####################################################################
# Compile operator
#####################################################################
# Gather QD-library source files and compile a library
add_library(IntegerFactorisation IntegerFactorisation.cpp IntegerFactorisation.h)
target_include_directories (IntegerFactorisation PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (IntegerFactorisation MCTDHlib)

