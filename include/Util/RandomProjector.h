//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_H
#define RANDOMPROJECTOR_H
#include <random>
#include "Core/Matrix.h"
#include "Core/Tensor.h"

namespace Random {

	template <typename T, class LinearOperator>
	Matrix<T> RandomQ(const LinearOperator& A, size_t k_plus_p,
		mt19937& gen);

	template <typename T, class LinearOperator>
	SpectralDecomposition<T> DiagonalizeRandom(const LinearOperator& A,
		size_t rank, size_t pow, mt19937& gen);
}

#endif //RANDOMPROJECTOR_H
