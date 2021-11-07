//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include "Core/Tensor_Functions_Implementation.h"

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

template Tensor<cd> Tensor_Extension::doubleHoleContraction(const Tensor<cd>& A,
	const Tensor<cd>& B, size_t k1, size_t k2);
template Tensor<d> Tensor_Extension::doubleHoleContraction(const Tensor<d>& A,
	const Tensor<d>& B, size_t k1, size_t k2);

/// new ... add double versions
template Tensor<cd> matrixTensor<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template void matrixTensor<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> tMatrixTensor<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode);
template Tensor<cd> multStateAB(const Matrix<cd>& A, const Tensor<cd>& B);
template Tensor<cd> multStateArTB<cd, cd>(const Matrix<cd>& A, const Tensor<cd>& B);
template void multStateArTB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B);
template void multStateAB<cd, cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, bool zero);
template Tensor<cd> project<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> projectOut<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<cd> projectOrthogonal<cd>(const Tensor<cd>& A, const Tensor<cd>& B);
template void multAdd<cd, cd>(Tensor<cd>& A, const Tensor<cd>& B, cd coeff);
template Matrix<cd> toMatrix(const Tensor<cd>& A);
template Matrix<cd> moveToMatrix(Tensor<cd>& A);
template Tensor<cd> toTensor(const Matrix<cd>& B);
template Matrix<cd> toMatrix(const Tensor<cd>& A, size_t mode);
template Tensor<cd> toTensor(const Matrix<cd>& B, const TensorShape& shape, size_t mode);
template Tensor<cd> moveToTensor(Matrix<cd>& B);
template void contraction<cd>(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t before, size_t active1, size_t active2, size_t behind);
template void matrixTensor<cd>(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>&  B, size_t before, size_t activeC, size_t activeB, size_t after, bool zero);
template Matrix<cd> contraction(const Tensor<cd>& A, const Tensor<cd>& B, size_t k);
template void contraction(Matrix<cd>& S, const Tensor<cd>& A, const Tensor<cd>& B, size_t k, bool zero);

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

