#####################################################################
# Compile operator
#####################################################################
# Gather QD-library source files and compile a library
add_library(NumberPartitioning NumberPartitioning.cpp NumberPartitioning.h NumberPartitioning_MLO.cpp NumberPartitioning_MLO.h)
target_include_directories (NumberPartitioning PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (NumberPartitioning MCTDHlib)

