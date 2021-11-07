#pragma once
/**
 * \defgroup QuTree
 * \brief QuTree is a C++ Tensor and Tensor Tree library.
 *
 * This library has classes that are similar to libraries like Lapack,
 * Eigen, Armadillo, etc., but its focus is concentrated on handling
 * high-order Tensors and Tensor Trees.
 * */

/// Core/
#include"Core/Matrix.h"
#include"Core/Matrix_Extension.h"
#include"Core/Matrix_Extension_Implementation.h"
#include"Core/Matrix_Implementation.h"
#include"Core/stdafx.h"
#include"Core/Tensor.h"
#include"Core/Tensor_Functions.h"
#include"Core/Tensor_Functions_Implementation.h"
#include"Core/Tensor_Implementation.h"
#include"Core/TensorShape.h"
#include"Core/Vector.h"
#include"Core/Vector_Implementation.h"

/// TreeClasses/
#include "TreeClasses/EdgeAttribute.h"
#include "TreeClasses/MatrixTensorTree.h"
#include "TreeClasses/MatrixTensorTreeFunctions.h"
#include "TreeClasses/MatrixTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions_Implementation.h"
#include "TreeClasses/NodeAttribute.h"
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/SparseMatrixTree.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions_Implementation.h"
#include "TreeClasses/SparseNodeAttribute.h"
#include "TreeClasses/SparseTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "TreeClasses/SymTensorTree.h"
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/TensorTreeFunctions_Implementation.h"
#include "TreeClasses/TreeIO.h"
#include "TreeClasses/TreeTransformation.h"
#include "TreeClasses/TreeTransformation_Implementation.h"

/// TreeOperators/
#include "TreeOperators/CoordinateTransformation.h"
#include "TreeOperators/LeafFunction.h"
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/LeafOperator.h"
#include "TreeOperators/MultiLeafOperator.h"
#include "TreeOperators/Potential.h"
#include "TreeOperators/PotentialOperator.h"
#include "TreeOperators/SOPVector.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SumOfProductsOperator_Implementation.h"

/// TreeShape/
#include "TreeShape/AbstractNode.h"
#include "TreeShape/Edge.h"
#include "TreeShape/Leaf.h"
#include "TreeShape/LinearizedLeaves.h"
#include "TreeShape/Node.h"
#include "TreeShape/NodePosition.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"

/// Util/
#include"Util/BS_integrator.h"
#include"Util/FFT.h"
#include"Util/FFTCooleyTukey.h"
#include"Util/JacobiRotationFramework.h"
#include"Util/Lanzcos.h"
#include"Util/QMConstants.h"
#include"Util/SimultaneousDiagonalization.h"
#include"../../QuTreeVM/src/Util/long_integer.h"

