//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_IMPLEMENTATION_H
#define RANDOMPROJECTOR_IMPLEMENTATION_H
#include "Util/RandomProjector.h"
//#include "Core/MatrixBLAS.h"

namespace Random {

	template<typename T>
	Matrix<T> randomGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim1 * dim2; ++i) {
			complex<double> x(dist(gen), dist(gen));
			g[i] = x;
		}
/*		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				complex<double> x(dist(gen), dist(gen));
				g(j, i) = x;
			}
		}*/
		return g;
	}

	template <typename T>
	Matrix<T> gue(size_t dim, mt19937& gen) {
		Matrix<T> r = randomGauss<T>(dim, dim, gen);
		return 0.5 * (r + r.adjoint());
	}

	template<typename T>
	Matrix<T> product(const Matrix<T>& x, const Matrix<T>& y, void *mem) {
		return x * y;
	}

	template <typename T, class LinearOperator, class Mem>
	Matrix<T> randomQ(const LinearOperator& A, size_t k_plus_p, mt19937& gen,
		Mem* mem) {
		assert(k_plus_p <= A.dim2());
		Matrix<T> Omega = randomGauss<T>(k_plus_p, A.dim1(), gen);
//		Matrix<T> Omega = GUE<T>(k_plus_p, A.Dim1(), gen);

//		Matrix<T> Y = A * Omega.adjoint();
		Matrix<T> Y = product(A, Omega.adjoint(), mem);
		/// Y = QR
		/// YY^ = QRR^Q^
		auto Q2 = qr(Y);
//		auto Q2(Y);
//		qrBLAS(Q2, Y);

		auto Q = subMatrix(Q2, Y.dim1(), Y.dim2());
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
//		auto Q = RandomGauss<T>(rank, A.Dim1(), gen);
			for (size_t i = 0; i < power; ++i) {
//			Matrix<T> Y = A * Q;
				Matrix<T> Y = product(A, Q, mem);
				auto Q2 = qr(Y);
//				auto Q2(Y);
//				qrBLAS(Q2, Y);
				Q = subMatrix(Q2, Y.dim1(), Y.dim2());
			}

//			auto Y = A * Q;
			Matrix<T> Y = product(A, Q, mem);

			/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
			auto B = Q.adjoint() * Y;
			auto x = diagonalize(B);
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
//		SpectralDecomposition<T> res(U, ew);
//		return res;
			return {U, ew};
		}
	}


#endif //RANDOMPROJECTOR_IMPLEMENTATION_H
