
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(scattering scattering.cpp)
target_include_directories (scattering PUBLIC .)
target_link_libraries (scattering MCTDHlib)

