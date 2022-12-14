
#####################################################################
# Compile operator
#####################################################################
# Gather QD-library source files and compile a library
add_library(CoupledHO CoupledHO.cpp CoupledHO.h)
target_include_directories (CoupledHO PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(CoupledHO PRIVATE cxx_std_14)
target_link_libraries (CoupledHO MCTDHlib)
target_link_libraries(CoupledHO QuTree::QuTree)

