
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(NOCl NOCl.cpp)
target_include_directories (NOCl PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(NOCl PRIVATE cxx_std_14)
target_link_libraries (NOCl MCTDHlib)

