# Easily regenerate with find src -name '*.cpp' | sort
set(QuTree_SOURCE_FILES

        # Tensor/
        src/Tensor/SVD.cpp
        src/Tensor/TensorShape.cpp
        src/Tensor/TensorBLAS1.cpp
        src/Tensor/TensorBLAS2.cpp
        src/Tensor/TensorLapack.cpp
        src/Tensor/TensorSlow.cpp
        include/Tensor/Tensor_Instantiations.cpp
        src/Tensor/stdafx.cpp

        # Primitve Basis
        src/Tree/PrimitiveBasis/DVR.cpp
        src/Tree/PrimitiveBasis/HarmonicOscillator.cpp
        src/Tree/PrimitiveBasis/LegendrePolynomials.cpp
        src/Tree/PrimitiveBasis/FFTGrid.cpp
        src/Tree/PrimitiveBasis/SpinGroup.cpp
        src/Tree/PrimitiveBasis/BasisAPI.cpp

        # Tree
        src/Tree/NodePosition.cpp
        src/Tree/Node.cpp
        src/Tree/LeafArray.cpp
        src/Tree/Edge.cpp
        src/Tree/SubTree.cpp
        src/Tree/Tree.cpp
        src/Tree/TreeFactory.cpp

        # TensorNetwork
        src/TensorNetwork/contractions.cpp
        src/TensorNetwork/TensorTree.cpp
        src/TensorNetwork/TensorTreeFactory.cpp

        # Operators
        src/Operator/LeafFunction.cpp
        src/Operator/LeafMatrix.cpp
        src/Operator/ProductOperator.cpp
        src/Operator/SumOfProductsOperator.cpp

        # Utilities
        src/Util/GateOperators.cpp
        src/Util/JacobiRotationFramework.cpp
        src/Util/SimultaneousDiagonalization.cpp
        )

