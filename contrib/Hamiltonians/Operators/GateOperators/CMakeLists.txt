
#####################################################################
# Compile PES
#####################################################################
# Gather QD-library source files and compile a library
add_library(GateOperators GateOperators.cpp
        QFT.cpp QFT.h
        Arithmetic.cpp Arithmetic.h Register.cpp Register.h)
target_include_directories (GateOperators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (GateOperators MCTDHlib)

