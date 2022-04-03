# Easily regenerated with find include -name '*.h' | sort
set(QuTree_INCLUDE_FILES
        include/QuTree.h
        include/Check.h
        include/stdafx.h

        # Tensor
        include/Tensor/Tensor
        include/Tensor/Tensor.h
        include/Tensor/TensorBLAS1.h
        include/Tensor/TensorBLAS2.h
        include/Tensor/TensorLapack.h
        include/Tensor/TensorShape.h
        include/Tensor/Tensor_Implementation.h

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
        )

