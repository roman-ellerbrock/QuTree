
#####################################################################
# Compile
#####################################################################
# Gather QD-library source files and compile a library
add_library(circuits
        Arithmetic.cpp
        Arithmetic.h
        GateOperators.cpp
        GateOperators.h
        QFT.h
        QFT.cpp
        Register.cpp
        Register.h
        Random.cpp
        Random.h
        QuantumCircuits.cpp
        QuantumCircuits.h
        NumberTheory.h
        NumberTheory.cpp
        Shor.h
        Shor.cpp
        ShorSparsetensor.cpp
        ShorSparsetensor.h
#       TensorOperatorTreeCircuits.cpp
        )

target_include_directories(circuits
        PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}
        ${CMAKE_CURRENT_SOURCE_DIR}/QVM/Circuits
        ${CMAKE_CURRENT_SOURCE_DIR}/QVM/Util
        ${CMAKE_CURRENT_SOURCE_DIR}/QVM
        )
target_link_libraries (circuits LINK_PUBLIC QuTree::QuTree yaml-cpp ${Boost_LIBRARIES})

