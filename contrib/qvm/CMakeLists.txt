
find_package(yaml-cpp REQUIRED)

include_directories(${YAML_CPP_INCLUDE_DIR})
add_library(Yaml INTERFACE IMPORTED)

find_package(Boost REQUIRED)

add_executable(qvm
        qvm_main.cpp
        QVM.cpp
        ConditionalQuantumInstruction.cpp
        IfElse.cpp
        InTimeQuantumInstruction.cpp
        MeasurementInstruction.cpp
        Measurements.cpp
        NumberTheory.cpp
        OutputInstruction.cpp
        QuantumCircuit.cpp
        QuantumInstruction.cpp
        FullRank.cpp
        RandomProjectorInstruction.cpp

        Circuits/Arithmetic.cpp
        Circuits/GateOperators.cpp
        Circuits/QFT.cpp
        Circuits/Register.cpp
        Circuits/Random.cpp
        Circuits/QuantumCircuits.cpp
        Circuits/Shor.cpp
        Circuits/ShorSparsetensor.cpp

        Util/filenames.cpp
        Util/long_integer.cpp
        Util/OverlapUtilities.cpp
        Util/QuadraticSieve.cpp
        Util/statisticalWavefunction.cpp
        )

target_compile_features(qvm PRIVATE cxx_std_14)
target_compile_options(qvm PUBLIC "$<$<CONFIG:DEBUG>:${qvm_DEBUG_FLAGS}>")
target_compile_options(qvm PRIVATE "$<$<CONFIG:RELEASE>:${qvm_RELEASE_FLAGS}>")
target_include_directories(qvm PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_options(qvm PUBLIC "-ffast-math" "-O3")
target_link_libraries(qvm QuTree::QuTree yaml-cpp Hamiltonians)

