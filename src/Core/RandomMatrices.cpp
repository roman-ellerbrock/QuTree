//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "Util/RandomMatrices.h"

namespace RandomMatrices {

	Matrixcd RandomRealGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				g(j, i) = dist(gen);
			}
		}
		return g;
	}

	Matrixcd RandomGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				g(j, i) = complex<double>(dist(gen), dist(gen));
			}
		}
		return g;
	}

	Matrixcd GUE(size_t dim, mt19937& gen) {
		Matrixcd r = RandomGauss(dim, dim, gen);
		return 0.5 * (r + r.Adjoint());
	}

	SpectralDecompositioncd GUE_diag(size_t dim, mt19937& gen) {
		auto A = GUE(dim, gen);
		return Diagonalize(A);
	}

	Matrixd GOE(size_t dim, mt19937& gen) {
		Matrixd r(dim, dim);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				r(j, i) = dist(gen);
			}
		}
		return 0.5 * (r + r.Adjoint());
	}

	SpectralDecompositiond GOE_diag(size_t dim, mt19937& gen) {
		auto A = GOE(dim, gen);
		return Diagonalize(A);
	}

	Matrixcd GUE(size_t dim1, size_t dim2, mt19937& gen) {
		size_t dim = max(dim1, dim2);
		auto G = GUE(dim, gen);
		auto Grect = Submatrix(G, dim1, dim2);
		return Grect;
	}

	Matrixcd RandomProjector(size_t dim1, size_t dim2, mt19937& gen) {
		/// Right hand side projector
		auto P = RandomRealGauss(dim1, dim2, gen);
//		P /= sqrt((double) dim2);
		for (size_t i = 0; i < P.Dim1(); ++i) {
			double norm = P.row(i).Norm();
			for (size_t j = 0; j < P.Dim2(); ++j) {
				P(i, j) /= norm;
			}
		}
		return P;
	}

	/// Build AP
	Matrixcd GUEProjector(size_t dim1, size_t dim2, mt19937& gen) {
		auto G = GUE(dim1, dim2, gen);
		auto Q = QR(G);
		return Submatrix(Q, dim1, dim2);
	}

	Matrixcd RandomQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen) {
		assert(k_plus_p <= A.Dim2());
		Matrixcd Omega = GUE(k_plus_p, A.Dim1(), gen);
		Matrixcd Y = A * Omega.Adjoint();
		/// Y = QR
		/// YY^ = QRR^Q^
		auto Q2 = QR(Y);
		auto Q = Submatrix(Q2, Y.Dim1(), Y.Dim1());
		return Q;
	}

	Matrixcd RandomProjection(const Matrixcd& A,
		size_t rdim, mt19937& gen) {
		assert(A.Dim1() == A.Dim2());
		Matrixcd Q = RandomQ(A, rdim, gen);
		return Q.Adjoint() * A * Q;
	}

	SpectralDecompositioncd DiagonalizeRandom(const Matrixcd& A,
		size_t rank, size_t p, mt19937& gen) {
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
		 * @param gen A random matrix generator
		 * @return Rectangular transformation matrix and eigenvalues
		 */

		Matrixcd Q = RandomQ(A, rank + p, gen);

		/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
		auto Aprime = Q.Adjoint() * A * Q;
		auto x = Diagonalize(Aprime);
		const Matrixcd& V = x.first;
		const Vectord& ew = x.second;

		/// Build eigenvector of full matrix
		auto U = Q * V;

		return {U, ew};
	}

	SVDcd svdRandom(const Matrixcd& A,
		size_t rank, mt19937& gen) {
		/**
		 * \brief Perform a randomized SVD
		 *
		 * For a detailed description see algorithm 5.1 in Ref. [1].
		 *
		 * [1] SIAM Rev., 53(2), 217–288. (72 pages)
		 *
		 * */
		auto Q = RandomQ(A, rank, gen);
		auto B = Q.Adjoint() * A;
		SVDcd Bsvd = svd(B);
		auto& U = get<0>(Bsvd);
		U = Q * U;
		return Bsvd;
	}

	Vectord probabilitiyDist(const Matrixcd& A) {
		auto tmp = A.diag();
		size_t dim = min(A.Dim1(), A.Dim2());
		Vectord p(dim);
		for (size_t i = 0; i < dim; ++i) {
			p(i) = pow(abs(tmp(i)), 2);
		}
		p /= sum(p);
		return p;
	}

	double entropy(const Vectord& p) {
		double S = 0;
		for (size_t i = 0; i < p.Dim(); ++i) {
			S -= p(i) * log(p(i));
		}
		return S;
	}

	double crossEntropy(const Vectord& p, const Vectord& q) {
		double H = 0;
		for (size_t i = 0; i < p.Dim(); ++i) {
			H -= p(i) * log(q(i));
		}
		return H;
	}
}

