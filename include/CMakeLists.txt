
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
        TreeHandling/AbstractNode.h
        TreeHandling/Leaf.h
        TreeHandling/LinearizedLeaves.h
        TreeHandling/Node.h
        TreeHandling/NodePosition.h
        TreeHandling/Tree.h
        TreeHandling/TTBasisFactory.h
        Tree/NodeAttribute.h
        Tree/TensorTree.h
        Tree/TensorTree_Implementation.h
        Tree/NodeAttribute.h
        Tree/SparseTree.h
        Tree/SpectralDecompositionTree.h
        Tree/IOTree.h
        Tree/MatrixTree.h
        Tree/MatrixTreeFunctions.h
        Tree/MatrixTreeFunctions_Implementation.h
        Tree/SparseMatrixTree.h
        Tree/SparseMatrixTreeFunctions.h
        Tree/SparseMatrixTreeFunctions_Implementation.h
        )

set(OPERATOR_FILES
        Operators/MultiLeafOperator.h
        Operators/LeafOperator.h
        Operators/LeafMatrix.h
        Operators/LeafFunction.h
        Operators/SumOfProductsOperator.h
        )

set(LEAF_TYPE_FILES
        TreeHandling/LeafTypes/LeafInterface.h
        TreeHandling/LeafTypes/BosonNumberBasis.h
        TreeHandling/LeafTypes/DVRBasis.h
        TreeHandling/LeafTypes/FermionNumberBasis.h
        TreeHandling/LeafTypes/FFTGrid.h
        TreeHandling/LeafTypes/HO_Basis.h
        TreeHandling/LeafTypes/LegendrePolynomials.h
        TreeHandling/LeafTypes/LogicalBasis.h
        TreeHandling/LeafTypes/NumberBasis.h
        TreeHandling/LeafTypes/SpinGroup.h
        )

install(FILES ${QDlib_INCLUDE_FILES} DESTINATION include/QDlib)
install(FILES ${TENSOR_TREE_FILES} DESTINATION include/QDlib/TensorTree)
install(FILES ${OPERATOR_FILES} DESTINATION include/QDlib/Operators)

