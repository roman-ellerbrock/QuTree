//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Extension_Implementation.h"

typedef complex<double> cd;
typedef double d;

// Tensor-Extension instantiations
template Matrix<cd> Tensor_Extension::outerProduct(const Tensor<cd>& A, const Tensor<cd>& B);
template void Tensor_Extension::OuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B);

template Matrix<cd> Tensor_Extension::weightedOuterProduct(const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& M);
template void Tensor_Extension::weightedOuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& rho);

template Matrix<cd> Tensor_Extension::map(const Tensor<cd>& A);
template Matrix<d> Tensor_Extension::map(const Tensor<d>& A);

template Tensor<cd> Tensor_Extension::regularize(Tensor<cd> A, size_t k, double eps);
template Tensor<d> Tensor_Extension::regularize(Tensor<d> A, size_t k, double eps);


template Tensor<cd> Tensor_Extension::doubleHoleContraction(const Tensor<cd>& A,
	const Tensor<cd>& B, size_t k1, size_t k2);
template Tensor<d> Tensor_Extension::doubleHoleContraction(const Tensor<d>& A,
	const Tensor<d>& B, size_t k1, size_t k2);

////////////////////////////////////////////////////////////
/// Direct Sump
////////////////////////////////////////////////////////////

template Tensor<d> Tensor_Extension::directSum(const Tensor<d>& A, const Tensor<d>& B, bool before, bool last);
template Tensor<cd> Tensor_Extension::directSum(const Tensor<cd>& A, const Tensor<cd>& B, bool before, bool last);

////////////////////////////////////////////////////////////
/// Random
////////////////////////////////////////////////////////////
template void Tensor_Extension::generateNormal(d* A, size_t n, mt19937& gen);
template void Tensor_Extension::generateNormal(cd* A, size_t n, mt19937& gen);
template void Tensor_Extension::generate<d>(Tensor<d>& A, mt19937& gen);
template void Tensor_Extension::generate<cd>(Tensor<cd>& A, mt19937& gen);
template void Tensor_Extension::generate<d>(Matrix<d>& A, mt19937& gen);
template void Tensor_Extension::generate<cd>(Matrix<cd>& A, mt19937& gen);
template void Tensor_Extension::generate<d>(Vector<d>& A, mt19937& gen);
template void Tensor_Extension::generate<cd>(Vector<cd>& A, mt19937& gen);

