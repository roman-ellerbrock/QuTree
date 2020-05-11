//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "Util/RandomMatrices.h"

namespace RandomMatrices {
	Matrixcd GUE(size_t dim, mt19937& gen) {
		Matrixcd r(dim, dim);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				r(j, i) = complex<double>(dist(gen), dist(gen));
			}
		}
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
		Matrixcd Grect(dim1, dim2);
		for (size_t c = 0; c < dim2; ++c) {
			for (size_t r = 0; r < dim1; ++r) {
				Grect(r, c) = G(r, c);
			}
		}
		return Grect;
	}

	Matrixcd QR_viaDiag(const Matrixcd Y) {
		/// Do Y = QR and return Q

		auto YY = Y * Y.Adjoint();
		auto x = Diagonalize(YY);
		Matrixcd Q(Y.Dim1(), Y.Dim2());
		size_t last = YY.Dim2() - 1;
		for (size_t r = 0; r < Y.Dim1(); ++r) {
			for (size_t c = 0; c < Y.Dim2(); ++c) {
				Q(r, c) = x.first(r, last - c);
			}
		}
		return Q;
	}

	Matrixcd RandomQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen) {
		assert(k_plus_p <= A.Dim2());
		Matrixcd Omega = GUE(k_plus_p, A.Dim1(), gen);
		Matrixcd Y = A * Omega.Adjoint();
//		Matrixcd Q = Y;
//		GramSchmidt(Q);
		/// Y = QR
		/// YY^ = QRR^Q^
//		auto Q = QR_viaDiag(Y);
		auto Q2 = QR(Y);
//		Q.print();
/*		cout << "Q2:\n";
		Q2.print();
		auto R = Q2.Adjoint() * Y;
		cout << "R:\n";
		R.print();*/

		Matrixcd Q(Y.Dim1(), Y.Dim2());
		for (size_t j = 0; j < Y.Dim2(); ++j) {
			for (size_t i = 0; i < Y.Dim1(); ++i) {
				Q(i, j) = Q2(i, j);
			}
		}
		return Q;
	}

	Matrixcd RandomProjection(const Matrixcd& A,
		size_t rdim, mt19937& gen) {
		Matrixcd Q = RandomQ(A, rdim, gen);
		return Q.Adjoint() * A * Q;
	}

	SpectralDecompositioncd DiagonalizeRandom(const Matrixcd& A,
		size_t rank, size_t p, mt19937& gen) {
		/**
		 * \brief Diagonalize using random projection
		 *
		 * For a detailed description see Ref. [1].
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
}

