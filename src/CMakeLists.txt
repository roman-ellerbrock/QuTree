
set(QDlib_SOURCE_FILES
        Core/TensorDim.cpp
        Core/TensorDim_Extension.cpp
        Core/FFT.cpp
        Core/JacobiRotationFramework.cpp
        Core/SimultaneousDiagonalization.cpp
        Core/long_integer.cpp
        Core/Vector_Instantiations.cpp
        Core/Matrix_Instantiations.cpp
        Core/Tensor_Instantiations.cpp
        Core/Tensor_Extension_Instantiations.cpp
        Core/FactorMatrix_Instantiations.cpp
        TensorTreeBasis/Leaf.cpp
        TensorTreeBasis/LinearizedLeaves.cpp
        TensorTreeBasis/Node.cpp
        TensorTreeBasis/NodePosition.cpp
        TensorTreeBasis/TensorTreeBasis.cpp
        Tree/TensorTree_Instantiation.cpp
        Tree/DenseOverlap.cpp
        Tree/HoleOverlap.cpp
        Tree/FactorMatrixTree.cpp
        Tree/HoleMatrixTree.cpp
        Operators/MultiParticleOperator.cpp
        Operators/SingleParticleOperatorFunction.cpp
        Operators/SingleParticleOperatorMatrix.cpp
        Operators/SumOfProductsOperator.cpp
        Core/stdafx.cpp
        )

set(LEAF_TYPE_CPP_FILES
        TensorTreeBasis/LeafTypes/BosonNumberBasis.cpp
        TensorTreeBasis/LeafTypes/DVRBasis.cpp
        TensorTreeBasis/LeafTypes/FermionNumberBasis.cpp
        TensorTreeBasis/LeafTypes/FFTGrid.cpp
        TensorTreeBasis/LeafTypes/HO_Basis.cpp
        TensorTreeBasis/LeafTypes/LegendrePolynomials.cpp
        TensorTreeBasis/LeafTypes/LogicalBasis.cpp
        TensorTreeBasis/LeafTypes/NumberBasis.cpp
        TensorTreeBasis/LeafTypes/SpinGroup.cpp
        )

# Compile objects only once for both static and shared libraries
# From https://stackoverflow.com/a/29824424/3052876
add_library(objlib OBJECT ${QDlib_SOURCE_FILES} ${LEAF_TYPE_CPP_FILES})
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)
add_library(QDlib_static STATIC $<TARGET_OBJECTS:objlib>)
add_library(QDlib SHARED $<TARGET_OBJECTS:objlib>)

# CMake <=2.8.7 compatible version to get both static and shared
#add_library(QDlib_static STATIC  ${QDlib_SOURCE_FILES} ${LEAF_TYPE_CPP_FILES})
#add_library(QDlib SHARED  ${QDlib_SOURCE_FILES} ${LEAF_TYPE_CPP_FILES})

SET_TARGET_PROPERTIES(QDlib_static PROPERTIES OUTPUT_NAME QDlib CLEAN_DIRECT_OUTPUT 1)

install(TARGETS QDlib DESTINATION lib)

