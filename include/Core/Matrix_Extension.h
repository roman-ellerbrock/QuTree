//
// Created by Roman Ellerbrock on 5/13/20.
//

#ifndef MATRIX_EXTENSION_H
#define MATRIX_EXTENSION_H
#include "Matrix.h"

template <typename T>
SpectralDecomposition<T> reduceRank(const SpectralDecomposition<T>& x, size_t rank);

template <typename T>
Matrix<T> merge(const Matrix<T>& A, const Matrix<T>& B,
	const Matrix<T>& AB);

#endif //MATRIX_EXTENSION_H
