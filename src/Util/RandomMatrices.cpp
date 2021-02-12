//
// Created by Roman Ellerbrock on 2/2/20.
//

#include "Util/RandomMatrices.h"

namespace RandomMatrices {

	Matrixcd randomRealGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				g(j, i) = dist(gen);
			}
		}
		return g;
	}

	Matrixcd randomGauss(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd g(dim1, dim2);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				g(j, i) = complex<double>(dist(gen), dist(gen));
			}
		}
		return g;
	}

	Vectorcd gaussVector(size_t dim, mt19937& gen) {
		normal_distribution<double> dist;
		Vectorcd r(dim);
		for (size_t i = 0; i < dim; ++i) {
			r(i) = dist(gen);
		}
		return r;
	}

	Matrixcd randomSparse(size_t dim1, size_t dim2, mt19937& gen) {
		Matrixcd ran(dim1, dim2);
		uniform_real_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim2; ++i) {
			for (size_t j = 0; j < dim1; ++j) {
				double r = dist(gen);
				if (r < 1./3.) {
					ran(j, i) = 1.;
				} else if (r < 2./3.) {
					ran(j, i) = 0.;
				} else {
					ran(j, i) = -1.;
				}
			}
		}
		return ran;
	}

	Matrixcd gue(size_t dim, mt19937& gen) {
		Matrixcd r = randomGauss(dim, dim, gen);
		return 0.5 * (r + r.adjoint());
	}

	SpectralDecompositioncd gueDiag(size_t dim, mt19937& gen) {
		auto A = gue(dim, gen);
		return diagonalize(A);
	}

	Matrixd goe(size_t dim, mt19937& gen) {
		Matrixd r(dim, dim);
		normal_distribution<double> dist(0., 1.);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				r(j, i) = dist(gen);
			}
		}
		return 0.5 * (r + r.adjoint());
	}

	SpectralDecompositiond goeDiag(size_t dim, mt19937& gen) {
		auto A = goe(dim, gen);
		return diagonalize(A);
	}

	Matrixcd gue(size_t dim1, size_t dim2, mt19937& gen) {
		size_t dim = max(dim1, dim2);
		auto G = gue(dim, gen);
		auto Grect = subMatrix(G, dim1, dim2);
		return Grect;
	}

	Matrixcd randomProjector(size_t dim1, size_t dim2, mt19937& gen) {
		/// Right hand side projector
		auto P = randomRealGauss(dim1, dim2, gen);
//		P /= sqrt((double) dim2);
		for (size_t i = 0; i < P.dim1(); ++i) {
			double norm = P.row(i).norm();
			for (size_t j = 0; j < P.dim2(); ++j) {
				P(i, j) /= norm;
			}
		}
		return P;
	}

	/// Build AP
	Matrixcd gueProjector(size_t dim1, size_t dim2, mt19937& gen) {
		auto G = gue(dim1, dim2, gen);
		auto Q = qr(G);
		return subMatrix(Q, dim1, dim2);
	}

	Matrixcd randomQ(const Matrixcd& A, size_t k_plus_p, mt19937& gen) {
		assert(k_plus_p <= A.dim2());
		Matrixcd Omega = randomSparse(k_plus_p, A.dim1(), gen);
//		Matrixcd Omega = GUE(k_plus_p, A.Dim1(), gen);
		Matrixcd Y = A * Omega.adjoint();
		/// Y = QR
		/// YY^ = QRR^Q^
		auto Q2 = qr(Y);

		auto Q = subMatrix(Q2, Y.dim1(), Y.dim2());
		return Q;
	}

	Matrixcd randomProjection(const Matrixcd& A,
		size_t rdim, mt19937& gen) {
		assert(A.dim1() == A.dim2());
		Matrixcd Q = randomQ(A, rdim, gen);
		return Q.adjoint() * A * Q;
	}

	SpectralDecompositioncd diagonalizeRandom(const Matrixcd& A,
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

		Matrixcd Q = randomQ(A, rank + p, gen);

		/// Build and Diagonalize Aprime = Q^* A Q = V ew V^*
		auto Aprime = Q.adjoint() * A * Q;
		auto x = diagonalize(Aprime);
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
		auto Q = randomQ(A, rank, gen);
		auto B = Q.adjoint() * A;
		SVDcd Bsvd = svd(B);
		auto& U = get<0>(Bsvd);
		U = Q * U;
		return Bsvd;
	}

//////////////////////////////////////////////////////////////
/// entropy & cross entropy
//////////////////////////////////////////////////////////////

	Vectord probabilitiyDist(const Matrixcd& A) {
		auto tmp = A.diag();
		size_t dim = min(A.dim1(), A.dim2());
		Vectord p(dim);
		for (size_t i = 0; i < dim; ++i) {
			p(i) = pow(abs(tmp(i)), 2);
		}
		p /= sum(p);
		return p;
	}

	double entropy(const Vectord& p) {
		double S = 0;
		for (size_t i = 0; i < p.dim(); ++i) {
			S -= p(i) * log(p(i));
		}
		return S;
	}

	double crossEntropy(const Vectord& p, const Vectord& q) {
		double H = 0;
		for (size_t i = 0; i < p.dim(); ++i) {
			H -= p(i) * log(q(i));
		}
		return H;
	}

	double entropy(const Matrixcd& A) {
		return entropy(probabilitiyDist(A));
	}

	double crossEntropy(const Matrixcd& p, const Matrixcd& q) {
		return crossEntropy(probabilitiyDist(p), probabilitiyDist(q));
	}

	double crossEntropyDifference(const Matrixcd& p, const Matrixcd& q) {
		double H = crossEntropy(p, q);
		double S = entropy(q);
		return (H - S);
	}

	template <typename T>
	void gramSchmidt(Vector<T>& v, const vector<Vector<T>>& es) {
		double conv = 0.;
		for (const auto& e : es) {
			T c = e * v;
			v -= e * c;
		}
		normalize(v);
	}

	vector<Vectorcd> buildKrylovSpace(Vectorcd x,
		const Matrixcd& A, size_t dim_subspace) {
		normalize(x);
		vector<Vectorcd> vec({x});
		for (size_t i = 0; i < dim_subspace; ++i) {
			Vectorcd Av = A * vec[i];
			gramSchmidt(Av, vec);
			vec.emplace_back(Av);
		}
		return vec;
	}

	Matrixcd toMatrix(const vector<Vectorcd>& x) {
		assert(!x.empty());
		size_t dim1 = x.front().dim();
		size_t dim2 = x.size();
		Matrixcd y(dim1, dim2);
		for (size_t j = 0; j < dim2; ++j) {
			const auto& a = x[j];
			for (size_t i = 0; i < dim1; ++i) {
				y(i, j) = a(i);
			}
		}
		return y;
	}
}

