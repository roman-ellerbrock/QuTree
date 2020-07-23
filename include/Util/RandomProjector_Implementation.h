//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_IMPLEMENTATION_H
#define RANDOMPROJECTOR_IMPLEMENTATION_H

namespace Random {

	template <typename T, class LinearOperator A>
	Matrix<T> RandomQ(const LinearOperator& A, size_t k_plus_p, mt19937& gen) {
		assert(k_plus_p <= A.Dim2());
		Matrixcd Omega = GUE(k_plus_p, A.Dim1(), gen);
		Matrixcd Y = A * Omega.Adjoint();
		/// Y = QR
		/// YY^ = QRR^Q^
		auto Q2 = QR(Y);

		auto Q = Submatrix(Q2, Y.Dim1(), Y.Dim2());
		return Q;
	}

	template <typename T, class LinearOperator A>
	SpectralDecomposition<T> DiagonalizeRandom(const LinearOperator& A,
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

		Matrix<T> Q = RandomQ(A, rank + p, gen);

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

}


#endif //RANDOMPROJECTOR_IMPLEMENTATION_H
