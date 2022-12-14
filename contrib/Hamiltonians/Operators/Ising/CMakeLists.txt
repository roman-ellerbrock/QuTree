
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(Ising Ising.cpp Ising.h Pauli.cpp Vertex2D.cpp)
target_include_directories (Ising PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (Ising MCTDHlib)

