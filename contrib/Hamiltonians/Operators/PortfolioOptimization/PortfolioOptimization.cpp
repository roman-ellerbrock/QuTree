//
// Created by Roman Ellerbrock on 5/3/22.
//

#include "PortfolioOptimization.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeClasses/Discrete/SymmetricSCF.h"
#include "Util/qutree_rng.h"
#include "TreeShape/TreeFactory.h"

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

double weight(size_t k, size_t Nq) {
	return (double) pow(2, k) / (double) (pow(2, Nq) - 1);
}

void addReturns(SOPcd& H, const Tensord& mu,
	double alpha, size_t Nt, size_t Na, size_t Nq) {

	LeafFuncd x = &LeafInterface::applyX;

	TensorShape shape_q({Nq, Na, Nt});

	for (size_t t = 0; t < Nt; ++t) {
		for (size_t i = 0; i < Na; ++i) {
			size_t muidx = indexMapping({t, i}, mu.shape_);
			for (size_t k = 0; k < Nq; ++k) {
				double c = -alpha * mu(muidx) * weight(k, Nq);
				cout << mu(muidx) << " ";

				size_t qidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(x, qidx);

				H.push_back(M, c);
			}
		}
	}
	cout << endl;
	getchar();
}

void addCovariance(SOPcd& H, const Tensord& cov, double gamma,
	size_t Nt, size_t Na, size_t Nq) {
	LeafFuncd x = &LeafInterface::applyX;

	TensorShape shape_q({Nq, Na, Nt});
	Matrixcd mat(Na, Na);
	for (size_t t = 0; t < Nt; ++t) { /// time
		/// diagonals
		for (size_t i = 0; i < Na; ++i) { /// asset i
			size_t sigidx = indexMapping({t, i, i}, cov.shape_);
			for (size_t k = 0; k < Nq; ++k) {
				for (size_t l = 0; l < Nq; ++l) {
					double c = 0.5 * gamma * cov(sigidx) * weight(k, Nq) * weight(l, Nq);
					cout << i << " - " << cov(sigidx) << "\n";

					size_t kidx = indexMapping({k, i, t}, shape_q);

					assert(kidx < shape_q.totalDimension());

					MLOcd M(x, kidx);

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
						cout << i << " " << j << " - " << cov(sigidx) << "\n";
						c *= 2.; // only upper triangle of cov used

						size_t kidx = indexMapping({k, i, t}, shape_q);
						size_t lidx = indexMapping({l, j, t}, shape_q);

						assert(kidx < shape_q.totalDimension());
						assert(lidx < shape_q.totalDimension());

						MLOcd M(x, kidx);
						M.push_back(x, lidx);

						H.push_back(M, c);
//						cout << i << " " << j << " | " << kidx << " " << lidx << " | " << c << endl;
						mat(i, j) = c;
					}
				}
			}
		}
	}
}

void addConstraint(SOPcd& H, const double rho, size_t Nt, size_t Na,
	size_t Nq, double K) {
	LeafFuncd x = &LeafInterface::applyX;
	size_t dim = pow(2, Nq);
	TensorShape shape_q({Nq, Na, Nt});

	for (size_t t = 0; t < Nt; ++t) {
		/// diagonals
/*		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t kidx = indexMapping({k, i, t}, shape_q);
				for (size_t l = 0; l < Nq; ++l) {
					size_t lidx = indexMapping({l, i, t}, shape_q);
					MLOcd M;
					M.push_back(x, kidx);
					M.push_back(x, lidx);
					double c = rho * weight(k, Nq) * weight(l, Nq);
					H.push_back(M, c);
				}
			}
		}*/

		/// off-diagonals
		for (size_t i = 0; i < Na; ++i) {
			for (size_t j = 0; j < Na; ++j) {
//			for (size_t j = i + 1; j < Na; ++j) {
				for (size_t k = 0; k < Nq; ++k) {
					size_t kidx = indexMapping({k, i, t}, shape_q);
					for (size_t l = 0; l < Nq; ++l) {
						size_t lidx = indexMapping({l, j, t}, shape_q);
						MLOcd M(x, kidx);
						M.push_back(x, lidx);
						double c = rho * weight(k, Nq) * weight(l, Nq);
//						c *= 2.; // only upper triangle of cov used
						H.push_back(M, c);
					}
				}
			}
		}

		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t qkidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(x, qkidx);
				H.push_back(M, -2. * K * rho * weight(k, Nq));
			}
		}

		for (size_t i = 0; i < Na; ++i) {
			for (size_t k = 0; k < Nq; ++k) {
				size_t kidx = indexMapping({k, i, t}, shape_q);
				MLOcd M(identityMatrixcd(dim), kidx);
				H.push_back(M, K * K / ((double) (Na * Nq)));
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
//	addCovariance(H, cov, gamma, Nt, Na, Nq);

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

Tensord to_Tensord(const vector<double>& a) {
	TensorShape shape({a.size(), 1});
	Tensord A(shape);
	for (size_t i = 0; i < a.size(); ++i) {
		A(i) = a[i];
	}
	return A;
}

Tensord selectDayReturns(const Tensord& mu, size_t day) {
	Tensord m({mu.shape_[1], 1});
	for (size_t i = 0; i < m.shape_[0]; ++i) {
		m(i) = mu(day, i);
	}
	return m;
}

auto selectDayCov(const Tensord& cov, size_t day) {
	const TensorShape& shape = cov.shape_;
	size_t dim = shape[1];
	Matrixd c(dim, dim);
	for (size_t i = 0; i < dim; ++i) {
		for (size_t j = 0; j < dim; ++j) {
			c(j, i) = cov({day, j, i});
		}
	}
	return c;
}

void meanVarianceAnalysisOptimization(string tickers,
	size_t Na, size_t Nt, size_t NaTot, size_t NtTot, size_t Nq,
	double alpha, double gamma, double rho, double K, const Tree& tree) {

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

	/// covariance
	Tensord cov_t = covariance(mu, mu_avg);

	/// select a day
	mu = selectDayReturns(mu, 0);

	Matrixd cov = selectDayCov(cov_t, 0);

	/**
	 * Quality of results:
	 * - compare to Gurobi
	 *
	 * Error of model:
	 * -
	 *
	 * Data:
	 * - NASDAQ
	 * - S&P 500
	 *
	 * Constraints:
	 * - spending of whole budget (done)
	 * - Portfolio must not be greater than a specific volatility
	 * - concentration of allocation (40-5-10)
	 * - transaction cost
	 * - permanent price shift
	 *
	 * Risk metrics:
	 * - covariance
	 * - non-convex
	 *
	 * Initialization:
	 * - initialize from convex solver
	 * - use only small intervall around convex solver (1 qubit/asset?)
	 *
	 * Problems:
	 * - constraint dominates optimization or is irrelevant
	 */

	auto f = [mu, cov, alpha, gamma, K, rho](const Configuration<>& c) {

		size_t Na = mu.shape_[0];
		size_t Nq = c.size() / Na;
		static vector<double> xs;
		static vector<size_t> x_int;
		static vector<size_t> tmp;
		if (xs.size() != Na) { xs.resize(Na); }
		if (x_int.size() != Na) { x_int.resize(Na); }
		if (tmp.size() != Nq) { tmp.resize(Nq); }
		split_doubles(xs, x_int, tmp, c, Na);
		auto x = to_Tensord(xs);

		double sum = 0;
		for (const auto& z : xs) { sum += z; }
		if (abs(sum - K) >= 0.999) { return rho * abs(sum - K); }

		double returns = - alpha * mu.dotProduct(x)[0];

		Tensord Cx = matrixTensor(cov, x, 0);
		Matrixd risk_mat = 0.5 * gamma * x.dotProduct(Cx);
		double risk = risk_mat[0];
//		if (risk > max_value_of_risk) { return large_penalty; }

		Tensord u = ones<double>({Na, 1});
		double inv = (u.dotProduct(x))[0];
		double investment = rho * abs(inv - K);

		double H = returns + risk + investment;
		return (double) H;
	};

	qutree::rng = mt19937(time(NULL));
	size_t run_max = 5;
	for (size_t run = 0; run < run_max; ++run) {
		auto Psi = randomConfigurationTree(tree, qutree::rng);
		auto c = optimize(Psi, f, tree, 5, 0);

		auto xs = split_doubles(c, Na);
		auto x = to_Tensord(xs);
		cout << "Final result:\n";
		double sum = 0.;
		for (const auto& z: xs) { sum += z; }
		cout << "Investment: " << sum << endl;
		cout << "Portfolio:\n";
		x.print();
		cout << "Normalized Portfolio:\n";
		x /= K;
		x.print();
		cout << "E = " << f(c) << endl;
	}
}

