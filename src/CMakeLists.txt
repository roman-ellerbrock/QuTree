# Easily regenerate with find src -name '*.cpp' | sort
set(QuTree_SOURCE_FILES
    src/Core/JacobiRotationFramework.cpp
    src/Core/Matrix_Instantiations.cpp
    src/Core/TensorShape.cpp
    src/Core/Tensor_Extension_Instantiations.cpp
    src/Core/Tensor_Instantiations.cpp
    src/Core/Vector_Instantiations.cpp
    src/Core/stdafx.cpp

    src/Core/LA_lib.f

    src/TreeClasses/MatrixTensorTree.cpp
    src/TreeClasses/MatrixTensorTreeFunctions.cpp
    src/TreeClasses/MatrixTree.cpp
    src/TreeClasses/MatrixTreeFunctions.cpp
    src/TreeClasses/TreeTransformations.cpp
    src/TreeClasses/SparseMatrixTree.cpp
    src/TreeClasses/SparseMatrixTreeFunctions.cpp
    src/TreeClasses/SparseTree.cpp
    src/TreeClasses/SparseTensorTree.cpp
    src/TreeClasses/SpectralDecompositionTree.cpp
    src/TreeClasses/SymTensorTree.cpp
    src/TreeClasses/TensorTree_Instantiation.cpp
    src/TreeClasses/TreeIO.cpp

    src/TreeOperators/LeafFunction.cpp
    src/TreeOperators/LeafMatrix.cpp
    src/TreeOperators/MultiLeafOperator.cpp
    src/TreeOperators/SumOfProductsOperator.cpp
    src/TreeOperators/TreeStructured/TreeSOP.cpp
    src/TreeOperators/TensorOperators/MatrixListTree.cpp
    src/TreeOperators/TensorOperators/TensorOperatorTree.cpp
    src/TreeOperators/TensorOperators/TensorOperatorTreeFunctions.cpp

    src/TreeShape/Leaf.cpp
    src/TreeShape/LeafTypes/DVRBasis.cpp
    src/TreeShape/LeafTypes/FFTGrid.cpp
    src/TreeShape/LeafTypes/HO_Basis.cpp
    src/TreeShape/LeafTypes/LegendrePolynomials.cpp
    src/TreeShape/LeafTypes/SpinGroup.cpp
    src/TreeShape/LinearizedLeaves.cpp
    src/TreeShape/Node.cpp
    src/TreeShape/NodePosition.cpp
    src/TreeShape/Tree.cpp
    src/TreeShape/TreeFactory.cpp

    src/Util/FFT.cpp
    src/Util/SimultaneousDiagonalization.cpp
    src/Util/RandomMatrices.cpp
    src/Util/string_ext.cpp
    src/Util/WeightedSimultaneousDiagonalization.cpp
)

