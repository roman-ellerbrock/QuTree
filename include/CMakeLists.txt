# Easily regenerated with find include -name '*.h' | sort
set(QuTree_INCLUDE_FILES
        QuTree.h

        include/Check.h

        # Tensor
        include/Tensor/Tensor.h
        include/Tensor/TensorBLAS1.h
        include/Tensor/TensorBLAS2.h
        include/Tensor/TensorShape.h
        include/Tensor/Tensor_Implementation.h
        include/Tensor/stdafx.h

        # Tree
        include/Tree/PrimitiveBasis/BasisAPI.h
        include/Tree/PrimitiveBasis/BasisParameters.h
        include/Tree/PrimitiveBasis/DVR.h
        include/Tree/PrimitiveBasis/FFTGrid.h
        include/Tree/PrimitiveBasis/LegendrePolynomials.h
        include/Tree/PrimitiveBasis/PrimitiveBasis.h
        include/Tree/PrimitiveBasis/SpinGroup.h
        include/Tree/Leaf.h
        include/Tree/NodePosition.h
        include/Tree/Node.h
        include/Tree/NodeArray.h
        include/Tree/Edge.h
        include/Tree/EdgeArray.h
        include/Tree/Tree.h
        include/Tree/Factory.h

        # TreeNetwork
        include/TensorNetwork/TensorTree.h
        include/TensorNetwork/MatrixTree.h

        # Operators
        include/Operator/LeafMatrix.h
        include/Operator/LeafFunction.h
        include/Operator/LeafOperator.h
        include/Operator/ProductOperator.h


        #    include/Core/Matrix.h
        #    include/Core/Matrix_Implementation.h
        #    include/Core/TensorMath.h
        #    include/Core/Tensor_Functions.h
        #    include/Core/Tensor_Functions_Implementation.h
        #    include/Core/Vector.h
        #    include/Core/Vector_Implementation.h
        #    include/QuTree.h

        #    include/TreeClasses/EdgeAttribute.h
        #    include/TreeClasses/MatrixTensorTree.h
        #    include/TreeClasses/MatrixTensorTreeFunctions.h
        #    include/TreeClasses/MatrixTree.h
        #    include/TreeClasses/MatrixTreeFunctions.h
        #    include/TreeClasses/MatrixTreeFunctions_Implementation.h
        #    include/TreeClasses/MatrixTreeTransformations.h
        #    include/TreeClasses/MatrixTreeTransformations_Implementation.h
        #    include/TreeClasses/NodeAttribute.h
        #    include/TreeClasses/SparseMatrixTree.h
        #    include/TreeClasses/SparseMatrixTreeFunctions.h
        #    include/TreeClasses/SparseMatrixTreeFunctions_Implementation.h
        #    include/TreeClasses/SparseNodeAttribute.h
        #    include/TreeClasses/SparseTensorTree.h
        #    include/TreeClasses/SparseTensorTree_Implementation.h
        #    include/TreeClasses/SparseTree.h
        #    include/TreeClasses/SOPMatrixTrees.h
        #    include/TreeClasses/SpectralDecompositionTree.h
        #    include/TreeClasses/SymTensorTree.h
        #    include/TreeClasses/TensorTree.h
        #    include/TreeClasses/TensorTree_Implementation.h
        #    include/TreeClasses/TreeIO.h
        #
        #    include/TreeOperators/LeafFunction.h
        #    include/TreeOperators/LeafMatrix.h
        #    include/TreeOperators/LeafOperator.h
        #    include/TreeOperators/MultiLeafOperator.h
        #    include/TreeOperators/PotentialOperator.h
        #    include/TreeOperators/SOPVector.h
        #    include/TreeOperators/SumOfProductsOperator.h
        #    include/TreeOperators/SumOfProductsOperator_Implementation.h
        #    include/TreeOperators/TensorOperators/TensorTreeOperator.h
        #    include/TreeOperators/TensorOperators/TensorOperatorRepresentation.h
        #    include/TreeOperators/TensorOperators/TTOrepresentation.h
        #    include/TreeOperators/TensorOperators/TTOcontraction.h
        #    include/TreeOperators/TensorOperators/TTOHoleTree.h
        #    include/TreeOperators/TensorOperators/TTOMatrixTree.h
        #
        #    include/TreeShape/AbstractNode.h
        #    include/TreeShape/Edge.h
        #    include/TreeShape/Leaf.h
        #    include/TreeShape/LeafTypes/DVRBasis.h
        #    include/TreeShape/LeafTypes/FFTGrid.h
        #    include/TreeShape/LeafTypes/HO_Basis.h
        #    include/TreeShape/LeafTypes/LeafInterface.h
        #    include/TreeShape/LeafTypes/LegendrePolynomials.h
        #    include/TreeShape/LeafTypes/SpinGroup.h
        #    include/TreeShape/LinearizedLeaves.h
        #    include/TreeShape/Node.h
        #    include/TreeShape/NodePosition.h
        #    include/TreeShape/Tree.h
        #    include/TreeShape/TreeFactory.h
        #
        #    include/Util/BS_integrator.h
        #    include/Util/FFT.h
        #    include/Util/FFTCooleyTukey.h
        #    include/Util/GradientDescent.h
        #    include/Util/GradientDescent_Implementation.h
        #    include/Util/JacobiRotationFramework.h
        #    include/Util/Lanzcos.h
        #    include/Util/MultiIndex.h
        #    include/Util/QMConstants.h
        #    include/Util/RandomMatrices.h
        #    include/Util/RandomProjector.h
        #    include/Util/RandomProjector_Implementation.h
        #    include/Util/SimultaneousDiagonalization.h
        #    include/Util/string_ext.h
        #    include/Util/Tree.h
        #    include/Util/TreeNode.h
        #    include/Util/WeightedSimultaneousDiagonalization.h
        )

