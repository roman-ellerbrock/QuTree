//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Extension_Implementation.h"

typedef complex<double> cd;
typedef double doub;

// Tensor-Extension instantiations
template Matrix<cd> Tensor_Extension::OuterProduct(const Tensor<cd>& A, const Tensor<cd>& B);
template void Tensor_Extension::OuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B);
template Matrix<cd> Tensor_Extension::WeightedOuterProduct(const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& m);
template void Tensor_Extension::WeightedOuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& rho);
template Tensor<cd> Tensor_Extension::Merge(Tensor<cd> A, const Tensor<cd>& B);
template Matrix<cd> Tensor_Extension::Map(const Tensor<cd>& A);

