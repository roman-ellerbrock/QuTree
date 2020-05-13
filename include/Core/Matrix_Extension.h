//
// Created by Roman Ellerbrock on 5/13/20.
//

#ifndef MATRIX_EXTENSION_H
#define MATRIX_EXTENSION_H
#include "Matrix.h"

template <typename T>
SpectralDecomposition<T> reduceRank(const SpectralDecomposition<T>& x, size_t rank) {
	size_t dim = x.second.Dim();
	Matrix<T> U
}

#endif //MATRIX_EXTENSION_H
