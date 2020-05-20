//
// Created by Roman Ellerbrock on 5/13/20.
//

#ifndef MATRIX_EXTENSION_IMPLEMENTATION_H
#define MATRIX_EXTENSION_IMPLEMENTATION_H
#include "Matrix_Extension.h"


template <typename T>
SpectralDecomposition<T> reduceRank(const SpectralDecomposition<T>& x,
	size_t rank) {
	size_t dim = x.second.Dim();
	assert(rank <= dim);
	Matrix<T> U(dim, rank);
	Vector<double> ew(rank);
	const Matrix<T>& evec = x.first;
	const Vectord& eval = x.second;
	size_t start = dim - rank;
	for (size_t i = 0; i < rank; ++i) {
		size_t shifted = i + start;
		for (size_t j = 0; j < dim; ++j) {
			U(j, i) = evec(j, shifted);
		}
		ew(i) = eval(shifted);
	}
	eval.print();
	ew.print();
	return {U, ew};
}

template<typename T>
Matrix<T> Merge(const Matrix<T>& A, const Matrix<T>& B,
	const Matrix<T>& AB) {
	// Merge Block matrices into one matrix
	// 	C =	(	A	AB	)
	// 		(	AB	B	)
	size_t n = A.Dim1();
	size_t m = B.Dim1();
	assert(n == AB.Dim1());
	assert(m == AB.Dim2());
	size_t N = m + n;

	Matrix<T> C(N, N);
	for (size_t i = 0; i < n; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C(j, i) = A(j, i);
		}
	}
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < m; ++j) {
			C(j + n, i + n) = B(j, i);
		}
	}
	// Off diagonals
	for (size_t i = 0; i < m; ++i) {
		for (size_t j = 0; j < n; ++j) {
			C(j, i + n) = AB(j, i);
			C(i + n, j) = conj(AB(j, i));
		}
	}
	return C;
}

#endif //MATRIX_EXTENSION_IMPLEMENTATION_H
