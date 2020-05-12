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

	Matrixcd RandomQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen);

	/// randomized Factorizations
	SpectralDecompositioncd DiagonalizeRandom(const Matrixcd& A,
		size_t rank, size_t p, mt19937& gen);

	Matrixcd BuildQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen);

	SVDcd svdRandom(const Matrixcd& A, size_t rank, mt19937& gen);

	Matrixcd randomSparse(size_t dim1, size_t dim2, mt19937& gen);

	/// Probability distribution framework
	/// Take absolute square of diagonals(A)
	Vectord probabilitiyDist(const Matrixcd& A);

	double entropy(const Vectord& p);

	double crossEntropy(const Vectord& p, const Vectord& q);

	double entropy(const Matrixcd& p);
	double crossEntropy(const Matrixcd& p, const Matrixcd& q);
	double crossEntropyDifference(const Matrixcd& p, const Matrixcd& q);


	/// Krylov subspace stuff
	vector<Vectorcd> BuildKrylovSpace(Vectorcd x,
	const Matrixcd& A, size_t dim_subspace);

	template <typename T>
	void GramSchmidt(Vector<T>& v, const vector<Vector<T>>& es);

	Matrixcd toMatrix(const vector<Vectorcd>& x);
}

#endif //RANDOMMATRICES_H
