# Easily regenerate with find src -name '*.cpp' | sort
set(QuTree_SOURCE_FILES

        # Tensor/
        src/Tensor/SVD.cpp
        src/Tensor/TensorShape.cpp
        src/Tensor/TensorBLAS1.cpp
        src/Tensor/TensorBLAS2.cpp
        src/Tensor/TensorLapack.cpp
        src/Tensor/TensorSlow.cpp
        src/Tensor/Tensor_Instantiations.cpp
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
        src/TensorNetwork/TensorTree.cpp
        src/TensorNetwork/MatrixTree.cpp
        src/TensorNetwork/DotProduct.cpp

        # Operators
        src/Operator/LeafFunction.cpp
        src/Operator/LeafMatrix.cpp
        src/Operator/ProductOperator.cpp

        #    src/Core/JacobiRotationFramework.cpp
        #    src/Core/LA_lib.f
        #    src/TreeClasses/MatrixTensorTree.cpp
        #    src/TreeClasses/MatrixTensorTreeFunctions.cpp
        #    src/TreeClasses/MatrixTree.cpp
        #    src/TreeClasses/MatrixTreeFunctions.cpp
        #    src/TreeClasses/TreeTransformations.cpp
        #    src/TreeClasses/SparseMatrixTree.cpp
        #    src/TreeClasses/SparseMatrixTreeFunctions.cpp
        #    src/TreeClasses/SparseTree.cpp
        #    src/TreeClasses/SparseTensorTree.cpp
        #    src/TreeClasses/SpectralDecompositionTree.cpp
        #    src/TreeClasses/SymTensorTree.cpp
        #    src/TreeClasses/TensorTree_Instantiation.cpp
        #    src/TreeClasses/TreeIO.cpp

        #    src/TreeOperators/LeafFunction.cpp
        #    src/TreeOperators/LeafMatrix.cpp
        #    src/TreeOperators/MultiLeafOperator.cpp
        #    src/TreeOperators/SumOfProductsOperator.cpp
        #    src/TreeOperators/TreeStructured/TreeSOP.cpp
        #    src/TreeOperators/TensorOperators/TensorTreeOperator.cpp
        #    src/TreeOperators/TensorOperators/TTOMatrixTree.cpp
        #    src/TreeOperators/TensorOperators/contractSOP.cpp
        #    src/TreeOperators/TensorOperators/contractCircuit.cpp
        #    src/TreeOperators/TensorOperators/TTOrepresentation.cpp
        #    src/TreeOperators/TensorOperators/TTNO_instantiations.cpp

        #    src/Util/FFT.cpp
        #    src/Util/GateOperators.cpp
        #    src/Util/SimultaneousDiagonalization.cpp
        #    src/Util/RandomMatrices.cpp
        #    src/Util/string_ext.cpp
        #    src/Util/WeightedSimultaneousDiagonalization.cpp
        )

