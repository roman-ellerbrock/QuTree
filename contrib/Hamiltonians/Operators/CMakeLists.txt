
#add_library(Operators CWrapper.f Operators.cpp GenericKEO/GenericKEO.cpp GenericKEO/GenericKEO.h Operator_lib.cpp Operator_lib.h NumberPartitioning/NumberPartitioning.cpp NumberPartitioning/NumberPartitioning.h IntegerFactorisation/IntegerFactorisation.h IntegerFactorisation/IntegerFactorisation.cpp UnsupervisedLearning/PictureInterface.cpp UnsupervisedLearning/PictureInterface.h)
#target_include_directories(Operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

#############################################
# List all available Operators here
# Note: Remember to link them to the Operators-library
#############################################
add_library(Operators Operators.h
        CoupledHO/CoupledHO.cpp
        CoupledHO/CoupledHO.h

        CH3_meanfield/CH3_meanfield.cpp
        CH3_meanfield/CH3_meanfield.h

        CH3_quasiexact/CH3_quasiexact.f
#        CH3_quasiexact/CH3_quasiexact.wrong.f
        CH3_quasiexact/CH3_quasiexact.cpp
        CH3_quasiexact/CH3_quasiexact.h
        CH3_quasiexact/TrafoCH3Quasiexact.cpp
        CH3_quasiexact/TrafoCH3Quasiexact.h

        ElectronicStructure/ElectronicStructure.cpp
        ElectronicStructure/ElectronicStructure.h
        ElectronicStructure/FermiGas.cpp
        ElectronicStructure/FermiGas.h
        ElectronicStructure/JordanWigner.cpp
        ElectronicStructure/JordanWigner.h

        NDI/NDIHamiltonian.cpp
        NDI/NDIHamiltonian.h

        Molecules/NormalModes.h
        Molecules/NOCl.cpp
        Molecules/NOCl.h

        PhysicalSystems/StandardOperatorLibrary.cpp
        PhysicalSystems/StandardOperatorLibrary.h
        PhysicalSystems/SpinBoson.cpp
        PhysicalSystems/SpinBoson.h

        PortfolioOptimization/PortfolioOptimization.cpp
        PortfolioOptimization/PortfolioOptimization.h

        schaepers/fortran/schaepers.f
        schaepers/schaepers.cpp
        schaepers/schaepers.h)

target_include_directories(Operators PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/../../src)

target_compile_features(Operators PRIVATE cxx_std_14)
target_compile_options(Operators PRIVATE "-Ofast")
target_link_libraries(Operators QuTree::QuTree)

#target_link_libraries(Operators LINK_PUBLIC
#        CoupledHO
#        )

#add_subdirectory("CH3_meanfield")
#add_subdirectory("CH3_quasiexact")
#add_subdirectory("CoupledHO")
#add_subdirectory("GateOperators")
#add_subdirectory("Ising")
#add_subdirectory("StandardKEO")
#add_subdirectory("NumberPartitioning")
#add_subdirectory("IntegerFactorisation")
#add_subdirectory("scattering")
#add_subdirectory("UnsupervisedLearning")
#
#target_link_libraries(Operators LINK_PUBLIC
#        CH3_quasiexact
#        CH3_meanfield
#        CoupledHO
#        GateOperators
#        Ising
#        NOCl
#        StandardKEO
#        NumberPartitioning
#        IntegerFactorisation
#        scattering
#        UnsupervisedLearning
#        )
