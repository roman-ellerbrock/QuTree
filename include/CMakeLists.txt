
set(QDlib_INCLUDE_FILES
        Util/BS_integrator.h
        Util/DoubleParticleOperator.h
        Util/FFT.h
        Util/FFTCooleyTukey.h
        Util/JacobiRotationFramework.h
        Util/Lanzcos.h
        Util/long_integer.h
        Core/Matrix.h
        Core/Matrix_Implementation.h
        Util/MultiIndex.h
        QDlib.h
        Util/QMConstants.h
        Util/SimultaneousDiagonalization.h
        Core/FactorMatrix.h
        Util/RandomMatrices.h
        Core/Tensor.h
        Core/TensorDim.h
        Core/TensorDim_Extension.h
        Core/Tensor_Extension.h
        Core/Tensor_Extension_Implementation.h
        Core/Tensor_Implementation.h
        Util/Tree.h
        Util/TreeNode.h
        Core/Vector.h
        Core/Vector_Implementation.h
        Util/WeightedSimultaneousDiagonalization.h
        Core/stdafx.h
        )

set(TENSOR_TREE_FILES
        TensorTreeBasis/AbstractNode.h
        TensorTreeBasis/Leaf.h
        TensorTreeBasis/LinearizedLeaves.h
        TensorTreeBasis/Node.h
        TensorTreeBasis/NodePosition.h
        TensorTreeBasis/TensorTreeBasis.h
        TensorTreeBasis/TTBasisFactory.h
        Tree/TreeStructuredObject.h
        Tree/TensorTree.h
        Tree/TensorTree_Implementation.h
        Tree/TreeStructuredObject.h
        Tree/TreeMarker.h
        Tree/SparseFactorMatrixTree.h
        Tree/SparseHoleMatrixTree.h
        Tree/SpectralDecompositionTree.h
        Tree/IOTree.h
        Tree/MatrixTree.h
        Tree/MatrixTreeFunctions.h
        Tree/MatrixTreeFunctions_Implementation.h
        Tree/SparseMatrixTree.h
        )

set(OPERATOR_FILES
        Operators/MultiLeafOperator.h
        Operators/LeafOperator.h
        Operators/LeafMatrix.h
        Operators/LeafFunction.h
        Operators/SumOfProductsOperator.h
        )

set(LEAF_TYPE_FILES
        TensorTreeBasis/LeafTypes/LeafInterface.h
        TensorTreeBasis/LeafTypes/BosonNumberBasis.h
        TensorTreeBasis/LeafTypes/DVRBasis.h
        TensorTreeBasis/LeafTypes/FermionNumberBasis.h
        TensorTreeBasis/LeafTypes/FFTGrid.h
        TensorTreeBasis/LeafTypes/HO_Basis.h
        TensorTreeBasis/LeafTypes/LegendrePolynomials.h
        TensorTreeBasis/LeafTypes/LogicalBasis.h
        TensorTreeBasis/LeafTypes/NumberBasis.h
        TensorTreeBasis/LeafTypes/SpinGroup.h
        )

install(FILES ${QDlib_INCLUDE_FILES} DESTINATION include/QDlib)
install(FILES ${TENSOR_TREE_FILES} DESTINATION include/QDlib/TensorTree)
install(FILES ${OPERATOR_FILES} DESTINATION include/QDlib/Operators)

