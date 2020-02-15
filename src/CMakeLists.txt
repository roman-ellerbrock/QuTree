
set(QDlib_SOURCE_FILES
        Core/TensorDim.cpp
        Core/TensorDim_Extension.cpp
        Core/FactorMatrix_Instantiations.cpp
        Core/JacobiRotationFramework.cpp
        Core/Vector_Instantiations.cpp
        Core/Matrix_Instantiations.cpp
        Core/RandomMatrices.cpp
        Core/Tensor_Instantiations.cpp
        Core/Tensor_Extension_Instantiations.cpp
        Core/stdafx.cpp
        Util/FFT.cpp
        Util/SimultaneousDiagonalization.cpp
        Util/long_integer.cpp
        Util/WeightedSimultaneousDiagonalization.cpp
        TreeOperators/MultiLeafOperator.cpp
        TreeOperators/LeafFunction.cpp
        TreeOperators/LeafMatrix.cpp
        TreeOperators/SumOfProductsOperator.cpp
        TreeShape/Leaf.cpp
        TreeShape/LinearizedLeaves.cpp
        TreeShape/Node.cpp
        TreeShape/NodePosition.cpp
        TreeShape/Tree.cpp
        TreeShape/TTBasisFactory.cpp
        TreeClasses/IOTree.cpp
        TreeClasses/MatrixTree.cpp
        TreeClasses/MatrixTreeFunctions.cpp
        TreeClasses/SparseMatrixTree.cpp
        TreeClasses/SparseMatrixTreeFunctions.cpp
        TreeClasses/SpectralDecompositionTree.cpp
        TreeClasses/TensorTree_Instantiation.cpp
        TreeClasses/SparseTree.cpp
        )

set(LEAF_TYPE_CPP_FILES
        TreeShape/LeafTypes/BosonNumberBasis.cpp
        TreeShape/LeafTypes/DVRBasis.cpp
        TreeShape/LeafTypes/FermionNumberBasis.cpp
        TreeShape/LeafTypes/FFTGrid.cpp
        TreeShape/LeafTypes/HO_Basis.cpp
        TreeShape/LeafTypes/LegendrePolynomials.cpp
        TreeShape/LeafTypes/LogicalBasis.cpp
        TreeShape/LeafTypes/NumberBasis.cpp
        TreeShape/LeafTypes/SpinGroup.cpp
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

