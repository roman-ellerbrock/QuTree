
set(QDlib_INCLUDE_FILES
        Util/BS_integrator.h
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
        QDlib.h
        )

set(TENSOR_TREE_FILES
        TreeShape/AbstractNode.h
        TreeShape/Leaf.h
        TreeShape/LinearizedLeaves.h
        TreeShape/Node.h
        TreeShape/NodePosition.h
        TreeShape/Tree.h
        TreeShape/TreeFactory.h
        TreeClasses/NodeAttribute.h
        TreeClasses/TensorTree.h
        TreeClasses/TensorTree_Implementation.h
        TreeClasses/NodeAttribute.h
        TreeClasses/SparseTree.h
        TreeClasses/SpectralDecompositionTree.h
        TreeClasses/IOTree.h
        TreeClasses/MatrixTree.h
        TreeClasses/MatrixTreeFunctions.h
        TreeClasses/MatrixTreeFunctions_Implementation.h
        TreeClasses/SparseMatrixTree.h
        TreeClasses/SparseMatrixTreeFunctions.h
        TreeClasses/SparseMatrixTreeFunctions_Implementation.h
        )

set(OPERATOR_FILES
        TreeOperators/MultiLeafOperator.h
        TreeOperators/LeafOperator.h
        TreeOperators/LeafMatrix.h
        TreeOperators/LeafFunction.h
        TreeOperators/SumOfProductsOperator.h
        )

set(LEAF_TYPE_FILES
        TreeShape/LeafTypes/LeafInterface.h
        TreeShape/LeafTypes/BosonNumberBasis.h
        TreeShape/LeafTypes/DVRBasis.h
        TreeShape/LeafTypes/FermionNumberBasis.h
        TreeShape/LeafTypes/FFTGrid.h
        TreeShape/LeafTypes/HO_Basis.h
        TreeShape/LeafTypes/LegendrePolynomials.h
        TreeShape/LeafTypes/LogicalBasis.h
        TreeShape/LeafTypes/NumberBasis.h
        TreeShape/LeafTypes/SpinGroup.h
        )

install(FILES ${QDlib_INCLUDE_FILES} DESTINATION include/QDlib)
install(FILES ${TENSOR_TREE_FILES} DESTINATION include/QDlib/TensorTree)
install(FILES ${OPERATOR_FILES} DESTINATION include/QDlib/Operators)

