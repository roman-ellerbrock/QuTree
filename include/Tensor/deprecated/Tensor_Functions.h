#pragma once
#include "Tensor/Tensor.h"
#include <random>

//! A Dot Product between two Tensors.
/*!
	This function calculates the doc-product between the tensors
	(*this) and A. The result is an overlap matrix sized with
	the number of tensors.
 */
template<typename T>
Matrix<T> dotProduct(const Tensor<T>& A, const Tensor<T>& A);

template<typename T, typename U>
void contraction1(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B,
	size_t before, size_t active1, size_t active2, size_t behind);

template<typename T>
void contraction(Matrix<T>& S, const Tensor<T>& A, const Tensor<T>& B, size_t k, bool zero = true);

template<typename T>
Matrix<T> contraction(const Tensor<T>& A, const Tensor<T>& B, size_t k);

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero = true);

template<typename T, typename U>
void tMatrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B,
	size_t before, size_t activeC, size_t activeB, size_t after, bool zero = true);

template<typename T, typename U>
void matrixTensor(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T, typename U>
Tensor<T> matrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

template<typename T, typename U>
void tensorMatrix(Tensor<T>& C, const Tensor<T>& B, const Matrix<U>& A, size_t mode, bool zero = true);

template<typename T, typename U>
Tensor<T> tensorMatrix(const Tensor<T>& B, const Matrix<U>& A, size_t mode);

template<typename T, typename U>
Tensor<T> tMatrixTensor(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

template<typename T, typename U>
void multStateAB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, bool zero = true);

template<typename T, typename U>
Tensor<T> multStateAB(const Matrix<U>& A, const Tensor<T>& B);

template<typename T, typename U>
void multStateArTB(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B);

template<typename T, typename U>
Tensor<T> multStateArTB(const Matrix<U>& A, const Tensor<T>& B);

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A);

template<typename T>
Matrix<T> toMatrix(const Tensor<T>& A, size_t mode);

template<typename T>
Matrix<T> moveToMatrix(Tensor<T>& A);

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B, const TensorShape& shape, size_t mode);

template<typename T>
Tensor<T> toTensor(const Matrix<T>& B);

template<typename T>
Tensor<T> moveToTensor(Matrix<T>& B);


//Projects B on A
template<typename T>
Tensor<T> project(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Tensor<T> projectOut(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Tensor< complex<double> > projectOrthogonal(const Tensor<complex<double> >& A,
	const Tensor< T >& B);








tuple<Tensorcd, Matrixcd, Vectord> SVD(const Tensorcd& A);

tuple<Matrixcd, Matrixcd, Vectord> SVD(const Matrixcd& A);

Tensorcd normalize(Tensorcd A, size_t k, double eps);

template<typename T>
Matrix<T> map(const Tensor<T>& A);

template<typename T>
Tensor<T> merge(Tensor<T> A, const Tensor<T>& B);

template<typename T>
void OuterProductAdd(Matrix<T>& M,
	const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Matrix<T> outerProduct(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Matrix<T> weightedOuterProduct(const Tensor<T>& A, const Tensor<T>& B,
	const Matrix<T>& M);

template<typename T>
void weightedOuterProductAdd(Matrix<T>& M, const Tensor<T>& A,
	const Tensor<T>& B, const Matrix<T>& rho);

template<typename T>
Tensor<T> doubleHoleContraction(const Tensor<T>& A, const Tensor<T>& B,
	size_t k1, size_t k2);

//////////////////////////////////////////////////////////////////////
/// Direct Sum + Product
//////////////////////////////////////////////////////////////////////

void shiftIndices(vector<size_t>& Ibreak, const TensorShape& shift,
	bool beforeLast, bool last);

TensorShape directSum(const TensorShape& A, const TensorShape& B,
	bool before, bool last);

template<typename T>
Tensor<T> directSum(const Tensor<T>& A, const Tensor<T>& B,
	bool before, bool last);

template<typename T>
Tensor<T> directProduct(const Tensor<T>& A, const Tensor<T>& B);

//////////////////////////////////////////////////////////////////////
/// Random number routines for tensors
//////////////////////////////////////////////////////////////////////

// Generate randomly occupied Tensor
template<typename T>
void generate(Tensor<T>& A, mt19937& gen);
// Generate randomly occupied Matrix
template<typename T>
void generate(Matrix<T>& A, mt19937& gen);
// Generate randomly occupied Vector
template<typename T>
void generate(Vector<T>& A, mt19937& gen);

template<typename T>
void generateNormal(T *A, size_t n, mt19937& gen);

