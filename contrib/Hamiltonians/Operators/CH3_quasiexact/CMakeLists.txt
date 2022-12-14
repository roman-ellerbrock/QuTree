
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(CH3_quasiexact CH3_quasiexact.cpp CH3_quasiexact_wrong.f TrafoCH3Quasiexact.cpp TrafoCH3Quasiexact.h)
target_include_directories (CH3_quasiexact PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (CH3_quasiexact Operators)

