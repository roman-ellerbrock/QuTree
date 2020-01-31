#pragma once
#include "Tensor.h"
#include <random>

namespace Tensor_Extension {

	tuple<Tensorcd, Matrixcd, Vectord> SVD(const Tensorcd& A);

	tuple<Matrixcd, Matrixcd, Vectord> SVD(const Matrixcd& A);

	template<typename T>
	Matrix<T> Map(const Tensor<T>& A);

	template<typename T>
	Tensor<T> Merge(Tensor<T> A, const Tensor<T>& B);

	template<typename T>
	void OuterProductAdd(Matrix<T>& M,
		const Tensor<T>& A, const Tensor<T>& B);

	template<typename T>
	Matrix<T> OuterProduct(const Tensor<T>& A, const Tensor<T>& B);

	template<typename T>
	Matrix<T> WeightedOuterProduct(const Tensor<T>& A, const Tensor<T>& B,
		const Matrix<T>& M);

	template<typename T>
	void WeightedOuterProductAdd(Matrix<T>& M, const Tensor<T>& A,
		const Tensor<T>& B, const Matrix<T>& rho);

	// Generate randomly occupied Tensor
	template<typename T>
	void Generate(Tensor<T>& A, mt19937& gen);
	// Generate randomly occupied Matrix
	template<typename T>
	void Generate(Matrix<T>& A, mt19937& gen);
	// Generate randomly occupied Vector
	template<typename T>
	void Generate(Vector<T>& A, mt19937& gen);

	template<typename T>
	void Generate_normal(T* A, size_t n, mt19937& gen);

	template<typename T>
	Matrix<T> OldStateAveragedHoleProduct(const Tensor<T>& A, const Tensor<T>& B, size_t k);

	template<typename T, typename U>
	Tensor<T> OldmultAB(const Matrix<U>& A, const Tensor<T>& B, size_t mode);

	template<typename T, typename U>
	Tensor<T> OldmultStateAB(const Matrix<U>& A, const Tensor<T>& B);
}
