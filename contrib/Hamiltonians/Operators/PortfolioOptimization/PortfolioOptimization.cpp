//
// Created by Roman Ellerbrock on 5/3/22.
//

#include "PortfolioOptimization.h"
#include "TreeOperators/SumOfProductsOperator.h"

Tensord readAssets(const string& name, size_t n_assets, size_t n_time) {
	TensorShape shape({n_time, n_assets});
	ifstream is(name);
	if (is.fail()) {
		cerr << "Error opening ticker file.\n";
		exit(1);
	}
	Tensord assets(shape);
	for (size_t j = 0; j < n_time; ++j) {
		for (size_t i = 0; i < n_assets; ++i) {
			is >> assets(j, i);
		}
	}
	return assets;
}

Tensord readAssets(const vector<string>& names, size_t m) {
	size_t N = names.size();

	Tensord assets({m, N});
	size_t j = 0;
	for (const string& name: names) {
		ifstream is(name);
		if (is.fail()) {
			cerr << "Problem opening asset file.\n";
			exit(0);
		}
		for (size_t i = 0; i < m; ++i) {
			is >> assets(i, j);
		}
		j++;
	}
	return assets;
}

template<typename T>
Tensor<T> log(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = log(A(i));
	}
	return A;
}

template<typename T>
Tensor<T> exp(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = exp(A(i));
	}
	return A;
}

Tensord log_returns(const Tensord& A) {
	size_t m = A.shape_.lastBefore() - 1;
	size_t N = A.shape_.lastDimension();
	TensorShape shape({m, N});
	Tensord mu(shape);
	for (size_t j = 0; j < N; ++j) {
		for (size_t i = 0; i < m; ++i) {
			mu(i, j) = log(A(i + 1, j) / A(i, j));
		}
	}
	return mu;
}

Tensord bare_returns(const Tensord& A) {
	const TensorShape& shape = A.shape_;
	size_t N = shape.lastDimension(); // n_assets
	size_t M = shape.lastBefore(); // n_assets
	Tensord mu({M, N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < M; ++m) {
			mu(m, n) = (A(m + 1, n) - A(m, n)) / A(m, n);
		}
	}
	return mu;
}

Tensord moving_avg(const Tensord& A, size_t range_back) {
	size_t M = A.shape_.lastBefore();
	size_t N = A.shape_.lastDimension();
	size_t L = M - range_back;
	Tensord av({M, N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t l = 0; l < L; ++l) {
			for (size_t k = 0; k < range_back; ++k) {
				av(l, n) += A(l + k, n);
			}
			av(l, n) /= (double) range_back;
		}
	}
	return av;
}

Tensord moving_covariance(const Tensord& mu, const Tensord& av, size_t range_back) {
	size_t M = mu.shape_.lastBefore();
	size_t N = mu.shape_.lastDimension();
	TensorShape shape({M, N, N});
	Tensord cov(shape);

	Matrixd sig(N, N);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			for (size_t m = 0; m < M; ++m) {
				cov({m, i, j}) = (mu(m, i) - av(m, i)) * (mu(m, j) - av(m, j));
				sig(i, j) += cov({m, i, j});
			}
		}
	}
	auto x = diagonalize(sig);
	cout << "avg(cov):\n";
	sig.print();
	cout << "lambda:\n";
	x.second.print();
	cout << "U:\n";
	x.first.print();
	return cov;
}

Tensord avg(const Tensord& A) {
	size_t M = A.shape_.lastBefore();
	size_t N = A.shape_.lastDimension();
	Tensord av({N});
	for (size_t n = 0; n < N; ++n) {
		for (size_t m = 0; m < M; ++m) {
			av(n) += A(m, n);
		}
		av(n) /= (double) M;
	}
	return av;
}

Tensord covariance(const Tensord& mu, const Tensord& av) {
	size_t M = mu.shape_.lastBefore();
	size_t N = mu.shape_.lastDimension();
	TensorShape shape({M, N, N});
	Tensord cov(shape);

	Matrixd sig(N, N);
	for (size_t i = 0; i < N; ++i) {
		for (size_t j = 0; j < N; ++j) {
			for (size_t m = 0; m < M; ++m) {
				cov({m, i, j}) = (mu(m, i) - av(i)) * (mu(m, j) - av(j));
				sig(i, j) += cov({m, i, j});
			}
		}
	}
	auto x = diagonalize(sig);
/*	cout << "avg(cov):\n";
	sig.print();
	cout << "lambda:\n";
	x.second.print();
	cout << "U:\n";
	x.first.print();*/
	return cov;
}

Tensord filter(const Tensord& A, size_t Delta) {
	const TensorShape& shape = A.shape_;
	size_t M = shape.lastBefore();
	size_t N = shape.lastDimension();
	size_t O = M / Delta;
	TensorShape sh({O, N});
	Tensord B(sh);
	cout << "filter:\n";
	shape.print();
	sh.print();
	for (size_t n = 0; n < N; ++n) {
		for (size_t o = 0; o < O; ++o) {
			size_t i = o * Delta;
			B(o, n) = A(i, n);
		}
	}
	return B;
}

Tensord selectAssets(const Tensord& A, size_t Na) {
	const TensorShape& shape = A.shape_;
	size_t M = shape.lastBefore();
	size_t N = shape.lastDimension();
	TensorShape sh({M, Na});
	Tensord B(sh);
	for (size_t n = 0; n < Na; ++n) {
		for (size_t i = 0; i < M; ++i) {
			B(i, n) = A(i, n);
		}
	}
	return B;
}

void diagonalizeCovariance(const Tensord& cov) {
	const TensorShape& shape = cov.shape_;
	Matrixd sig(shape[1], shape[2]);
	Matrixd sig_av(shape[1], shape[2]);
	for (size_t m = 0; m < shape[0]; ++m) {
		for (size_t i = 0; i < sig.dim1(); ++i) {
			for (size_t j = 0; j < sig.dim2(); ++j) {
				sig(i, j) = cov({m, i, j});
			}
		}
		sig_av += sig;
	}

	sig_av /= (double) shape[0];
	auto y = diagonalize(sig_av);
/*	cout << "cov_av:\n";
	sig_av.print();
	cout << "cov eigenvalues:\n";
	y.second.print();*/
}

Matrixcd set1() {
	Matrixcd s(2, 2);
	s(1, 1) = 1.;
	return s;
}

double weight(size_t k, size_t Nq) {
	return (double) pow(2, k) / (double) (pow(2, Nq) - 1);
}

void addReturns(SOPcd& H, const Tensord& mu, double alpha, size_t Nt, size_t Na, size_t Nq) {
	TensorShape shape_q({Nq, Na, Nt});
//	cout << "Nt = " << Nt << endl;
	for (size_t t = 0; t < Nt; ++t) {
//		cout << "t = " << t << endl;
		for (size_t i = 0; i < Na; ++i) {
			size_t muidx = indexMapping({t, i}, mu.shape_);
			for (size_t k = 0; k < Nq; ++k) {
				/// total weight including negative returns
				double c = -alpha * mu(muidx) * weight(k, Nq);
//				cout << i << "\t" << c << endl;

				size_t qidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(set1(), qidx);

				H.push_back(M, c);
			}
		}
	}
}

void addCovariance(SOPcd& H, const Tensord& cov, double gamma, size_t Nt, size_t Na, size_t Nq) {
	TensorShape shape_q({Nq, Na, Nt});
	Matrixcd mat(Na, Na);
	for (size_t t = 0; t < Nt; ++t) { /// time
		/// diagonals
		for (size_t i = 0; i < Na; ++i) { /// asset i
			size_t sigidx = indexMapping({t, i, i}, cov.shape_);
			for (size_t k = 0; k < Nq; ++k) {
				for (size_t l = 0; l < Nq; ++l) {
					double c = 0.5 * gamma * cov(sigidx) * weight(k, Nq) * weight(l, Nq);

					size_t kidx = indexMapping({k, i, t}, shape_q);

					assert(kidx < shape_q.totalDimension());

					MLOcd M(set1(), kidx);

					H.push_back(M, c);
					mat(i, i) = c;
				}
			}
		}

		/// offdiagonals
		for (size_t i = 0; i < Na; ++i) { /// asset i
			for (size_t j = i + 1; j < Na; ++j) { /// asset j
				/// qubit indices
				size_t sigidx = indexMapping({t, i, j}, cov.shape_);
				for (size_t k = 0; k < Nq; ++k) {
					for (size_t l = 0; l < Nq; ++l) {
						double c = 0.5 * gamma * cov(sigidx) * weight(k, Nq) * weight(l, Nq);
						c *= 2.; // only upper triangle of cov used

						size_t kidx = indexMapping({k, i, t}, shape_q);
						size_t lidx = indexMapping({l, j, t}, shape_q);

						assert(kidx < shape_q.totalDimension());
						assert(lidx < shape_q.totalDimension());

						MLOcd M(set1(), kidx);
						M.push_back(set1(), lidx);

						H.push_back(M, c);
//						cout << i << " " << j << " | " << kidx << " " << lidx << " | " << c << endl;
						mat(i, j) = c;
					}
				}
			}
		}
	}
}

void addConstraint(SOPcd& H, const double rho, size_t Nt, size_t Na, size_t Nq, double K) {
	TensorShape shape_q({Nq, Na, Nt});

	for (size_t t = 0; t < Nt; ++t) {
		/// diagonals
		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t kidx = indexMapping({k, i, t}, shape_q);
				for (size_t l = 0; l < Nq; ++l) {
					size_t lidx = indexMapping({l, i, t}, shape_q);
					MLOcd M(set1(), kidx);
					M.push_back(set1(), lidx);
					double c = rho * weight(k, Nq) * weight(l, Nq);
					H.push_back(M, c);
				}
			}
		}

		/// off-diagonals
		for (size_t i = 0; i < Na; ++i) {
			for (size_t j = i + 1; j < Na; ++j) {
				for (size_t k = 0; k < Nq; ++k) {
					size_t kidx = indexMapping({k, i, t}, shape_q);
					for (size_t l = 0; l < Nq; ++l) {
						size_t lidx = indexMapping({l, j, t}, shape_q);
						MLOcd M(set1(), kidx);
						M.push_back(set1(), lidx);
						double c = rho * weight(k, Nq) * weight(l, Nq);
						c *= 2.; // only upper triangle of cov used
						H.push_back(M, c);
					}
				}
			}
		}
	}

	///
	for (size_t t = 0; t < Nt; ++t) {
		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t qkidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(set1(), qkidx);
				H.push_back(M, -2. * K * rho * weight(k, Nq));
			}
		}
	}

	///
	for (size_t t = 0; t < Nt; ++t) {
		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t kidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(identityMatrixcd(2), kidx);
				H.push_back(M, K * K / ((double) (Na * Nt * Nq)));
			}
		}
	}
}

vector<string> read_tickers(const string& filename) {
	ifstream is(filename);
	vector<string> names;
	for (std::string line; std::getline(is, line);) {
		names.push_back(line + ".csv");
		cout << line << endl;
	}
	return names;
}

SOPcd meanVarianceAnalysis(const Tensord& mu, const Tensord& cov,
	size_t Na, size_t Nt, size_t Nq,
	double alpha, double gamma, double rho, double K) {

	SOPcd H;
	/// Add mu
	addReturns(H, mu, alpha, Nt, Na, Nq);
	cout << "# |H| after returns: " << H.size() << " | expected: " << Na * Nt * Nq << endl;

	/// Add sigma
	addCovariance(H, cov, gamma, Nt, Na, Nq);

	size_t expect_sig = Na * Nt * Nq + Na * Na * Nt * Nq * Nq;
	cout << "# |H| after covariance: " << H.size() << " | expected: " << expect_sig << endl;

	/// Add constraint
	addConstraint(H, rho, Nt, Na, Nq, K);
	size_t expect_rho = 3 * Na * Nt * Nq + 2 * Na * Na * Nt * Nq * Nq;
	cout << "# |H| after constraint: " << H.size() << " | expected: " << expect_rho << endl;

	return H;
}

SOPcd meanVarianceAnalysis(string tickers,
	size_t Na, size_t Nt, size_t NaTot, size_t NtTot, size_t Nq,
	double alpha, double gamma, double rho, double K) {

	tickers.erase(0, tickers.find_first_not_of("\n/ "));
	cout << "# Assets: " << Na << " / " << NaTot << endl;
	cout << "# Time steps: " << Nt << " / " << NtTot << endl;
	cout << "# Qubits per asset: " << Nq << endl;
	cout << "# H = -" << alpha << " mu * omega + ";
	cout << gamma << " / 2 * omega^T* sigma * omega + ";
	cout << rho << " * (u^T omega - " << K << ")^2" << endl;

	/// read & select assets
	auto A = readAssets(tickers, NaTot, NtTot);
	A = selectAssets(A, Na);

	/// returns
	auto mu = log_returns(A);

	/// average returns
	Tensord mu_avg = avg(mu);

	auto cov = covariance(mu, mu_avg);

	return meanVarianceAnalysis(mu, cov, Na, Nt, Nq, alpha, gamma, rho, K);
}
