//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef RANDOMMATRICES_H
#define RANDOMMATRICES_H
#include "Core/Matrix.h"
#include <random>

namespace RandomMatrices {
	Matrixcd GUE(size_t dim, mt19937& gen);
	Matrixd GOE(size_t dim, mt19937& gen);
	Matrixcd BuildQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen);
	Matrixcd RandomProjection(const Matrixcd A,
		size_t stat_dim, mt19937& gen);
}

#endif //RANDOMMATRICES_H
