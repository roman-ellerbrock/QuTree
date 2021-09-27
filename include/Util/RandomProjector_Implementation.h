//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_IMPLEMENTATION_H
#define RANDOMPROJECTOR_IMPLEMENTATION_H
#include "Util/RandomProjector.h"
#include "Core/MatrixBLAS.h"

namespace Random {

	/**
	 * \brief Draw a random matrix from the gaussian unitary ensemble.
	 * @tparam T base type
	 * @param dim1 number of rows of the matrix
	 * @param dim2 number of columns of the matrix
	 * @param gen random number generator
	 * @return Matrix with random gaussian unitary entries.
	 */
	template<typename T>
	Matrix<T> randomGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim1 * dim2; ++i) {
			complex<double> x(dist(gen), dist(gen));
			g[i] = x;
		}
		return g;
	}

	/**
	 * @tparam T base type
	 * @param dim dimension of square matrix
	 * @param gen random number generator
	 * @return GUE matrix created by symmetrizing a random Gaussian matrix.
	 */
	template<typename T>
	Matrix<T> gue(size_t dim, mt19937& gen) {
		Matrix<T> r = randomGauss<T>(dim, dim, gen);
		return 0.5 * (r + r.adjoint());
	}

	template<typename T>
	Matrix<T> product(const Matrix<T>& x, const Matrix<T>& y, void *mem) {
		return x * y;
	}

	/**
	 * Create matrix Q as required by random algorithms in Ref. [1].
	 *
	 * [1] SIAM Rev., 53(2), 217–288. (72 pages)
	 *
	 * @tparam Mem container for memory attached to linear operator
	 * @param A the linear operator used for preconditioning the random matrix
	 * @param k_plus_p dimension of the random space. Refers to notation in Ref. [1]
	 * @return Random Matrix Q as required by Ref. [1].
	 */
	template<typename T, class LinearOperator, class Mem>
	Matrix<T> randomQ(const LinearOperator& A, size_t k_plus_p, mt19937& gen,
		Mem *mem) {
		assert(k_plus_p <= A.dim2());
		Matrix<T> Omega = randomGauss<T>(k_plus_p, A.dim1(), gen);

		Matrix<T> Y = product(A, Omega.adjoint(), mem);
		/// Y = QR
		/// YY^ = QRR^Q^
//		auto Q2 = qr(Y); // for non-BLAS version of QuTree
		auto Q2(Y);
		qrBLAS(Q2, Y);

		auto Q = subMatrix(Q2, Y.dim1(), Y.dim2());
		return Q;
	}

	template<typename T, class LinearOperator, class Mem>
	SpectralDecomposition<T> diagonalizeRandom(const LinearOperator& A,
		size_t rank, size_t power, mt19937& gen, Mem *mem) {
		/**
		 * \brief Diagonalize using random projection
		 *
		 * For a detailed description see algorithm 5.3 in Ref. [1].
		 *
		 * [1] SIAM Rev., 53(2), 217–288. (72 pages)
		 *
		 * @param A Matrix that should be diagonalized
		 * @param rank target rank of the matrix
		 * @param p oversampling parameter (extra dimension)
		 * @param power exponent for power iteration
		 * @param gen A random matrix generator
		 * @return Rectangular transformation matrix and eigenvalues
		 */

		auto Q = randomQ<T, LinearOperator, Mem>(A, rank, gen, mem);
		for (size_t i = 0; i < power; ++i) {
			Matrix<T> Y = product(A, Q, mem);
			auto Q2(Y);
			qrBLAS(Q2, Y);
			Q = subMatrix(Q2, Y.dim1(), Y.dim2());
		}

		Matrix<T> Y = product(A, Q, mem);

		/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
		auto B = Q.adjoint() * Y;
		auto x = diagonalize(B);
		const Matrix<T>& V = x.first;
		auto& ew = x.second;
		/// Build eigenvector of full matrix
		auto U = Q * V;
		return {U, ew};
	}
}


#endif //RANDOMPROJECTOR_IMPLEMENTATION_H
