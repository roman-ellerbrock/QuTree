
set(QDlib_INCLUDE_FILES
        Core/BS_integrator.h
        Core/DoubleParticleOperator.h
        Core/FFT.h
        Core/FFTCooleyTukey.h
        Core/JacobiRotationFramework.h
        Core/Lanzcos.h
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
        Core/long_integer.h
        Core/stdafx.h
        )

set(TENSOR_TREE_FILES
        TensorTree/AbstractNode.h
        TensorTree/Leaf.h
        TensorTree/LinearizedLeaves.h
        TensorTree/Node.h
        TensorTree/NodePosition.h
        TensorTree/TensorTreeBasis.h
        TensorTree/TreeStructuredObject.h
        TensorTree/TensorTree.h
        TensorTree/TensorTree_Implementation.h
        TensorTree/DenseOverlap.h
        TensorTree/TreeStructuredObject.h
        TensorTree/TreeMarker.h
        TensorTree/FactorMatrixTree.h
        TensorTree/HoleMatrixTree.h
        )

set(OPERATOR_FILES
        Operators/MultiParticleOperator.h
        Operators/SingleParticleOperator.h
        Operators/SingleParticleOperatorMatrix.h
        Operators/SingleParticleOperatorFunction.h
        Operators/SumOfProductsOperator.h
        )

set(LEAF_TYPE_FILES
        TensorTree/LeafTypes/PrimitiveBasis.h
        TensorTree/LeafTypes/BosonNumberBasis.h
        TensorTree/LeafTypes/DVRBasis.h
        TensorTree/LeafTypes/FermionNumberBasis.h
        TensorTree/LeafTypes/FFTGrid.h
        TensorTree/LeafTypes/HO_Basis.h
        TensorTree/LeafTypes/LegendrePolynomials.h
        TensorTree/LeafTypes/LogicalBasis.h
        TensorTree/LeafTypes/NumberBasis.h
        TensorTree/LeafTypes/SpinGroup.h
        )

install(FILES ${QDlib_INCLUDE_FILES} DESTINATION include/QDlib)
install(FILES ${TENSOR_TREE_FILES} DESTINATION include/QDlib/TensorTree)
install(FILES ${OPERATOR_FILES} DESTINATION include/QDlib/Operators)

