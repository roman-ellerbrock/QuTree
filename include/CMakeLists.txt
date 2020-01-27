
set(QDlib_INCLUDE_FILES
        BS_integrator.h
        DoubleParticleOperator.h
        FFT.h
        FFTCooleyTukey.h
        JacobiRotationFramework.h
        Lanzcos.h
        Matrix.h
        Matrix_Implementation.h
        MultiIndex.h
        QDlib.h
        QMConstants.h
        SimultaneousDiagonalization.h
        FactorMatrix.h
        Tensor.h
        TensorDim.h
        TensorDim_Extension.h
        Tensor_Extension.h
        Tensor_Extension_Implementation.h
        Tensor_Implementation.h
        Tree.h
        TreeNode.h
        Vector.h
        Vector_Implementation.h
        long_integer.h
        stdafx.h
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
        TensorTree/HMatrices.h
        TensorTree/HHoleMatrices.h
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

