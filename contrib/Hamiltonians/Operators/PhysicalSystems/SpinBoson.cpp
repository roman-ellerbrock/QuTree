//
// Created by Roman Ellerbrock on 8/8/20.
//
#include "SpinBoson.h"

namespace Operator {

	Matrixd Coefficients(const string& filename, const Tree& tree) {
		size_t n_fragment = tree.nLeaves();
		size_t dim = 0;
		for (size_t n = 0; n < n_fragment; ++n) {
			const Leaf& leaf = tree.getLeaf(n);
			size_t n_states = leaf.dim();
			dim += n_states;
		}
		Matrixd C(dim, dim);
		ifstream is(filename);
		for (size_t i = 0; i < dim; ++i) {
			for (size_t j = 0; j < dim; ++j) {
				is >> C(j, i);
			}
		}
		return C;
	}

	Matrixd Submatrix(const Matrixd& C, size_t m, size_t n,
		const Tree& tree) {

		/// Copy Matrix from supermatrix
		size_t n_fragment = tree.nLeaves();
		assert(n < n_fragment && m < n_fragment);
		const Leaf& leaf1 = tree.getLeaf(m);
		const Leaf& leaf2 = tree.getLeaf(n);
		size_t dim1 = leaf1.dim();
		size_t dim2 = leaf2.dim();

		size_t start1 = 0;
		for (size_t l = 0; l < m; ++l) {
			const Leaf& leaf = tree.getLeaf(l);
			size_t dim = leaf.dim();
			start1 += dim;
		}
		size_t start2 = 0;
		for (size_t l = 0; l < n; ++l) {
			const Leaf& leaf = tree.getLeaf(l);
			size_t dim = leaf.dim();
			start2 += dim;
		}

		cout << "dims: " << dim1 << " " << dim2 << endl;
		cout << "starts: " << start1 << " " << start2 << endl;
		Matrixd sub(dim1, dim2);
		for (size_t i = 0; i < dim1; ++i) {
			for (size_t j = 0; j < dim2; ++j) {
				sub(j, i) = C(start1 + j, start2 + i);
			}
		}
		return sub;
	}

	SOPcd Exciton(const string& filename, const Tree& tree) {
		size_t n_fragments = tree.nLeaves();

		/// Get coefficients from file
		Matrixd C = Coefficients(filename, tree);
		auto x = diagonalize(C);
		cout << "Eigenvalues (Singles approx.): " << endl;
		x.second *= 27.2114;
		x.second.print();

		/// Append terms
		SOPcd H;
/*
		for (size_t n = 0; n < n_fragments; ++n) {
			for (size_t m = 0; m < n_fragments; ++m) {
//				if (m != n ) {continue;}
					Matrixd lambda = Submatrix(C, m, n, tree);

					auto x = svd(lambda);
					const Matrixd U = get<0>(x);
					const Matrixd V = get<1>(x);
					const Vectord alpha = get<2>(x);

					if (m != n) {
						lambda *= 0.5;
					}
					cout << "Coupling Framents: " << m << " and " << n << endl;
					lambda.print();
					cout << "singular values:" << endl;
					alpha.print();

					for (size_t k = 0; k < alpha.Dim(); ++k) {
						MLOcd M;
						Vectord u = U.col(k);
						Vectord v = V.col(k);
						cout << "k = " << k << endl;
						u.print();
						v.print();
						shared_ptr<LeafOperatorcd> hu = make_shared<VectorProjector>(u);
						shared_ptr<LeafOperatorcd> hv = make_shared<VectorProjector>(v);
						M.push_back(hu, m);
						M.push_back(hv, n);
						H.push_back(M, alpha(k));
					}
				}
			}

		*/

		/*
		for (size_t i = 0; i < lambda.Dim1(); ++i) {
			shared_ptr<LeafOperatorcd> Pi = make_shared<PrimitiveProjector>(i);
			for (size_t j = 0; j < lambda.Dim2(); ++j) {
				shared_ptr<LeafOperatorcd> Pj = make_shared<PrimitiveProjector>(j);
				if (abs(lambda(i, j)) > 1e-10) {
					MLOcd M;
					M.push_back(Pi, m);
					M.push_back(Pj, n);
					H.push_back(M, lambda(i, j));
				}
			}
		}
		 */

		for (size_t n = 0; n < n_fragments; ++n) {
			for (size_t m = n; m < n_fragments; ++m) {
				Matrixd lambda = Submatrix(C, m, n, tree);

				if (n == m) {
					cout << "Coupling Framents: " << m << " and " << n << endl;
//					lambda.print();
					Vectord diag = lambda.diag();
					diag.print();
					MLOcd M;
					shared_ptr<LeafOperatorcd> e = make_shared<VectorProjector>(diag);
					M.push_back(e, m);
					H.push_back(M, 1.);

/*					for (size_t i = 0; i < lambda.Dim1(); ++i) {
						shared_ptr<LeafOperatorcd> Xi = make_shared<ExtendedX>(i, 0);
						for (size_t j = 0; j < lambda.Dim2(); ++j) {
							shared_ptr<LeafOperatorcd> Xj = make_shared<ExtendedX>(j, 0);
							if (abs(lambda(i, j)) > 1e-10) {
								MLOcd M;
								M.push_back(Xi, m);
								M.push_back(Xj, n);
								H.push_back(M, lambda(i, j));
							}
						}
					}
					*/
				} else {

					for (size_t i = 0; i < lambda.dim1(); ++i) {
						shared_ptr<LeafOperatorcd> Xi = make_shared<Excite>(i, 0);
						Vectord alpha = lambda.col(i);
						shared_ptr<LeafOperatorcd> Xjs = make_shared<VectorExcite>(0, alpha);
						MLOcd M;
						M.push_back(Xi, m);
						M.push_back(Xjs, n);
						H.push_back(M, 1.);
					}

					for (size_t j = 0; j < lambda.dim2(); ++j) {
						shared_ptr<LeafOperatorcd> Xj = make_shared<Excite>(j, 0);
						Vectord beta = lambda.row(j);
						shared_ptr<LeafOperatorcd> Xis = make_shared<VectorExcite>(0, beta);
						MLOcd M;
						M.push_back(Xj, n);
						M.push_back(Xis, m);
						H.push_back(M, 1.);
					}
				}
/*
//					lambda *= 0.5;
					for (size_t i = 0; i < lambda.Dim1(); ++i) {
						shared_ptr<LeafOperatorcd> Xi = make_shared<Excite>(i, 0);
						for (size_t j = 0; j < lambda.Dim2(); ++j) {
							shared_ptr<LeafOperatorcd> Xj = make_shared<Excite>(0, j);
							if (abs(lambda(i, j)) > 1e-7) {
								MLOcd M;
								M.push_back(Xi, m);
								M.push_back(Xj, n);
								H.push_back(M, lambda(i, j));
							}
						}
					}
					//
					for (size_t i = 0; i < lambda.Dim1(); ++i) {
						shared_ptr<LeafOperatorcd> Xi = make_shared<Excite>(0, i);
						for (size_t j = 0; j < lambda.Dim2(); ++j) {
							shared_ptr<LeafOperatorcd> Xj = make_shared<Excite>(j, 0);
							if (abs(lambda(i, j)) > 1e-7) {
								MLOcd M;
								M.push_back(Xi, m);
								M.push_back(Xj, n);
								H.push_back(M, lambda(i, j));
							}
						}
					}

				}
			}
			*/
			}
		}
		cout << "Hamiltonian size: " << H.size() << endl;

		return H;
	}
}

