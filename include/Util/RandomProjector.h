//
// Created by Roman Ellerbrock on 7/14/20.
//

#ifndef RANDOMPROJECTOR_H
#define RANDOMPROJECTOR_H

namespace Random {

	template <typename T>
	SpectralDecompositioncd DiagonalizeRandom(const Matrix<T>& A,
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
		auto Aprime = Q.Adjoint() * A * Q;
		auto x = Diagonalize(Aprime);
		const Matrixcd& V = x.first;
		const Vectord& ew = x.second;

		/// Build eigenvector of full matrix
		auto U = Q * V;

		return {U, ew};
	}

}

#endif //RANDOMPROJECTOR_H
