//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_IMPLEMENTATION_H
#define RANDOMPROJECTOR_IMPLEMENTATION_H
#include "Util/RandomProjector.h"

namespace Random {

	template <typename T>
	Matrix<T> RandomGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				complex<double> x(dist(gen), dist(gen));
				g(j, i) = x;
			}
		}
		return g;
	}

	template <typename T>
	Matrix<T> GUE(size_t dim, mt19937& gen) {
		Matrix<T> r = RandomGauss<T>(dim, dim, gen);
		return 0.5 * (r + r.Adjoint());
	}

	template <typename T, class LinearOperator>
	Matrix<T> RandomQ(const LinearOperator& A, size_t k_plus_p, mt19937& gen) {
		assert(k_plus_p <= A.Dim2());
		Matrix<T> Omega = RandomGauss<T>(k_plus_p, A.Dim1(), gen);
//		Matrix<T> Omega = GUE<T>(k_plus_p, A.Dim1(), gen);
		Matrix<T> Y = A * Omega.Adjoint();
		/// Y = QR
		/// YY^ = QRR^Q^
		auto Q2 = QR(Y);

		auto Q = Submatrix(Q2, Y.Dim1(), Y.Dim2());
		return Q;
	}

/*	template <typename T, class LinearOperator>
	SpectralDecomposition<T> DiagonalizeRandom(const LinearOperator& A,
		size_t rank, size_t p, mt19937& gen) {
		/ **
		 * \brief Diagonalize using random projection
		 *
		 * For a detailed description see algorithm 5.3 in Ref. [1].
		 *
		 * [1] SIAM Rev., 53(2), 217–288. (72 pages)
		 *
		 * @param A Matrix that should be diagonalized
		 * @param rank target rank of the matrix
		 * @param p oversampling parameter (extra dimension)
		 * @param gen A random matrix generator
		 * @return Rectangular transformation matrix and eigenvalues
		 * /

		auto Q = RandomQ<T, LinearOperator>(A, rank + p, gen);

		/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
		auto B = A * Q;
		auto Aprime = Q.Adjoint() * B;
		auto x = Diagonalize(Aprime);
		const Matrixcd& V = x.first;
		const Vectord& ew = x.second;

		/// Build eigenvector of full matrix
		auto U = Q * V;

		return {U, ew};
	}
	*/

// dimensions: 6000*160000
// rank: 1500 or 3000

	template <typename T, class LinearOperator>
	SpectralDecomposition<T> DiagonalizeRandom(const LinearOperator& A,
		size_t rank, size_t power, mt19937& gen) {
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

		auto Q = RandomQ<T, LinearOperator>(A, rank, gen);
//		auto Q = RandomGauss<T>(rank, A.Dim1(), gen);
		for (size_t i = 0; i < power; ++i) {
			Matrix<T> Y = A * Q;
			auto Q2 = QR(Y);
			Q = Submatrix(Q2, Y.Dim1(), Y.Dim2());
		}

		auto Y = A * Q;
/*		for (size_t i = 1; i < power; ++i) {
			Y = A * Y;
		}*/

		/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
		auto B = Q.Adjoint() * Y;
		auto x = Diagonalize(B);
		const Matrix<T>& V = x.first;
		auto& ew = x.second;
/*		if (power != 1) {
			for (size_t i = 0; i < ew.Dim(); ++i) {
				ew(i) = pow(ew(i), 1. / ((double) power));
			}
		}
*/
		/// Build eigenvector of full matrix
		auto U = Q * V;

		return {U, ew};
	}
}


#endif //RANDOMPROJECTOR_IMPLEMENTATION_H
