//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Extension_Implementation.h"

typedef complex<double> cd;
typedef double d;

// Tensor-Extension instantiations
template Matrix<cd> Tensor_Extension::OuterProduct(const Tensor<cd>& A, const Tensor<cd>& B);
template void Tensor_Extension::OuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B);

template Matrix<cd> Tensor_Extension::WeightedOuterProduct(const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& m);
template void Tensor_Extension::WeightedOuterProductAdd(Matrixcd& M, const Tensor<cd>& A, const Tensor<cd>& B,
	const Matrix<cd>& rho);

template Matrix<cd> Tensor_Extension::Map(const Tensor<cd>& A);
template Matrix<d> Tensor_Extension::Map(const Tensor<d>& A);

template Tensor<cd> Tensor_Extension::DoubleHoleContraction(const Tensor<cd>& A,
	const Tensor<cd>& B, size_t k1, size_t k2);
template Tensor<d> Tensor_Extension::DoubleHoleContraction(const Tensor<d>& A,
	const Tensor<d>& B, size_t k1, size_t k2);

////////////////////////////////////////////////////////////
/// Direct Sump
////////////////////////////////////////////////////////////

template Tensor<d> Tensor_Extension::DirectSum(const Tensor<d>& A, const Tensor<d>& B, bool before, bool last);
template Tensor<cd> Tensor_Extension::DirectSum(const Tensor<cd>& A, const Tensor<cd>& B, bool before, bool last);

////////////////////////////////////////////////////////////
/// Random
////////////////////////////////////////////////////////////
template void Tensor_Extension::Generate_normal(d* A, size_t n, mt19937& gen);
template void Tensor_Extension::Generate_normal(cd* A, size_t n, mt19937& gen);
template void Tensor_Extension::Generate<d>(Tensor<d>& A, mt19937& gen);
template void Tensor_Extension::Generate<cd>(Tensor<cd>& A, mt19937& gen);
template void Tensor_Extension::Generate<d>(Matrix<d>& A, mt19937& gen);
template void Tensor_Extension::Generate<cd>(Matrix<cd>& A, mt19937& gen);
template void Tensor_Extension::Generate<d>(Vector<d>& A, mt19937& gen);
template void Tensor_Extension::Generate<cd>(Vector<cd>& A, mt19937& gen);

