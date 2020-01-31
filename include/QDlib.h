#include"Core/BS_integrator.h"
#include"Core/DoubleParticleOperator.h"
#include"Core/FFT.h"
#include"Core/FFTCooleyTukey.h"
#include"Core/JacobiRotationFramework.h"
#include"Core/Lanzcos.h"
#include"Core/Matrix.h"
#include"Core/Matrix_Implementation.h"
#include"Core/MultiIndex.h"
#include"Core/QMConstants.h"
#include"Core/SimultaneousDiagonalization.h"
#include"Core/FactorMatrix.h"
#include"Core/Tensor.h"
#include"Core/TensorDim.h"
#include"Core/TensorDim_Extension.h"
#include"Core/Tensor_Extension.h"
#include"Core/Tensor_Extension_Implementation.h"
#include"Core/Tensor_Implementation.h"
#include"Core/Tree.h"
#include"Core/TreeNode.h"
#include"Core/Vector.h"
#include"Core/Vector_Implementation.h"
#include"Core/long_integer.h"
#include"Core/stdafx.h"

/**
 * \defgroup QDlib
 * \brief QDlib is a C++ Tensor and Tensor Tree library.
 *
 * This library has classes that are similar to libraries like Lapack,
 * Eigen, Armadillo, etc., but its focus is concentrated on handling
 * high-order Tensors and Tensor Trees.
 * */
