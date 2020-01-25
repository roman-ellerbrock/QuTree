
set(QDlib_SOURCE_FILES
        TensorDim.cpp
        TensorDim_Extension.cpp
        FFT.cpp
        JacobiRotationFramework.cpp
        SimultaneousDiagonalization.cpp
        CWrapper.f
        long_integer.cpp
        Vector_Instantiations.cpp
        Matrix_Instantiations.cpp
        Tensor_Instantiations.cpp
        Tensor_Extension_Instantiations.cpp
        FactorMatrix_Instantiations.cpp
        TensorTree/Leaf.cpp
        TensorTree/LinearizedLeaves.cpp
        TensorTree/Node.cpp
        TensorTree/NodePosition.cpp
        TensorTree/TensorTreeBasis.cpp
        TensorTree/TensorTree_Instantiation.cpp
        TensorTree/DenseOverlap.cpp
        TensorTree/HoleOverlap.cpp
        Operators/SingleParticleOperatorFunction.cpp
        Operators/MultiParticleOperator.cpp
        stdafx.cpp)

set(LEAF_TYPE_CPP_FILES
        TensorTree/LeafTypes/BosonNumberBasis.cpp
        TensorTree/LeafTypes/DVRBasis.cpp
        TensorTree/LeafTypes/FermionNumberBasis.cpp
        TensorTree/LeafTypes/FFTGrid.cpp
        TensorTree/LeafTypes/HO_Basis.cpp
        TensorTree/LeafTypes/LegendrePolynomials.cpp
        TensorTree/LeafTypes/LogicalBasis.cpp
        TensorTree/LeafTypes/NumberBasis.cpp
        TensorTree/LeafTypes/SpinGroup.cpp
        )

add_library(QDlib ${QDlib_SOURCE_FILES} ${LEAF_TYPE_CPP_FILES})
install(TARGETS QDlib DESTINATION lib)

