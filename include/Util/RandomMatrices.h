//
// Created by Roman Ellerbrock on 2/2/20.
//

#ifndef RANDOMMATRICES_H
#define RANDOMMATRICES_H
#include "Core/Matrix.h"
#include <random>

namespace RandomMatrices {

	Matrixcd RandomRealGauss(size_t dim1, size_t dim2, mt19937& gen);
	Matrixcd RandomGauss(size_t dim1, size_t dim2, mt19937& gen);

	Matrixcd RandomProjector(size_t dim1, size_t dim2, mt19937& gen);

	Matrixcd GUE(size_t dim, mt19937& gen);

	Matrixd GOE(size_t dim, mt19937& gen);

	Matrixcd GUEProjector(size_t dim1, size_t dim2, mt19937& gen);

	Matrixcd RandomProjection(const Matrixcd& A,
		size_t stat_dim, mt19937& gen);


	/// randomized Factorizations
	SpectralDecompositioncd DiagonalizeRandom(const Matrixcd& A,
		size_t rank, size_t p, mt19937& gen);

	Matrixcd BuildQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen);

	SVDcd svdRandom(const Matrixcd& A, size_t rank, mt19937& gen);


	/// Probability distribution framework
	/// Take absolute square of diagonals(A)
	Vectord probabilitiyDist(const Matrixcd& A);

	double entropy(const Vectord& p);

	double crossEntropy(const Vectord& p, const Vectord& q);
}

#endif //RANDOMMATRICES_H
