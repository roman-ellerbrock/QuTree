#pragma once
#include "Tensor.h"
#include <random>

namespace Tensor_Extension {

	tuple<Tensorcd, Matrixcd, Vectord> SVD(const Tensorcd& A);

	tuple<Matrixcd, Matrixcd, Vectord> SVD(const Matrixcd& A);

	template <typename T>
	Tensor<T> regularize(Tensor<T> A, size_t k, double eps);

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
	void generateNormal(T* A, size_t n, mt19937& gen);

	template<typename T>
	Matrix<T> oldStateAveragedHoleProduct(const Tensor<T>& A, const Tensor<T>& B, size_t k);

	template<typename T, typename U>
	Tensor<T> oldMultAB(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

	template<typename T, typename U>
	Tensor<T> oldmultStateAb(const Matrix<U>& A, const Tensor<T>& B);
}
