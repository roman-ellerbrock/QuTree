
set(QDlib_INCLUDE_FILES
        Core/BS_integrator.h
        Core/DoubleParticleOperator.h
        Core/FFT.h
        Core/FFTCooleyTukey.h
        Core/JacobiRotationFramework.h
        Core/Lanzcos.h
        Core/long_integer.h
        Core/Matrix.h
        Core/Matrix_Implementation.h
        Core/MultiIndex.h
        QDlib.h
        Core/QMConstants.h
        Core/SimultaneousDiagonalization.h
        Core/FactorMatrix.h
        Core/Tensor.h
        Core/TensorDim.h
        Core/TensorDim_Extension.h
        Core/Tensor_Extension.h
        Core/Tensor_Extension_Implementation.h
        Core/Tensor_Implementation.h
        Core/Tree.h
        Core/TreeNode.h
        Core/Vector.h
        Core/Vector_Implementation.h
        Core/WeightedSimultaneousDiagonalization.h
        Core/stdafx.h
        )

set(TENSOR_TREE_FILES
        TensorTreeBasis/AbstractNode.h
        TensorTreeBasis/Leaf.h
        TensorTreeBasis/LinearizedLeaves.h
        TensorTreeBasis/Node.h
        TensorTreeBasis/NodePosition.h
        TensorTreeBasis/TensorTreeBasis.h
        Tree/TreeStructuredObject.h
        Tree/TensorTree.h
        Tree/TensorTree_Implementation.h
        Tree/FactorMatrixTree.h
        Tree/TreeStructuredObject.h
        Tree/TreeMarker.h
        Tree/SparseFactorMatrixTree.h
        Tree/SparseHoleMatrixTree.h
        Tree/SpectralDecompositionTree.h
        )

set(OPERATOR_FILES
        Operators/MultiParticleOperator.h
        Operators/SingleParticleOperator.h
        Operators/SingleParticleOperatorMatrix.h
        Operators/SingleParticleOperatorFunction.h
        Operators/SumOfProductsOperator.h
        )

set(LEAF_TYPE_FILES
        TensorTreeBasis/LeafTypes/PrimitiveBasis.h
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

