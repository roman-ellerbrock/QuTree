//
// Created by Roman Ellerbrock on 9/22/22.
//

#include "NumberTheory.h"
#include <iostream>

namespace number_theory {
	void add(map<cpp_int, size_t>& x, cpp_int p) {
		if (x.find(p) == x.end()) {
			x[p] = 1;
		} else {
			x[p] += 1;
		}
	}

	PrimeFactorization factorize(cpp_int N) {
		PrimeFactorization x;
		while (N != 1) {
			for (cpp_int p = 2; p <= N; ++p) {
				if (N % p == 0) {
					add(x, p);
					N /= p;
					break;
				}
			}
		}
		return x;
	}

	vector<PrimeFactorization> factorize(const vector<cpp_int>& Ns) {
		vector<PrimeFactorization> fs;
		for (const auto& N: Ns) {
			fs.push_back(factorize(N));
		}
		return fs;
	}

	bool is_square_free(const PrimeFactorization& f) {
		return std::all_of(f.begin(), f.end(), [](const pair<cpp_int, size_t>& x) { return (x.second < 2); });
	}

	vector<cpp_int> coprimes(const cpp_int& N) {
		vector<cpp_int> z;
		for (cpp_int i = 1; i < N; ++i) {
			if (gcd(i, N) == 1) {
				z.push_back(N);
			}
		}
		return z;
	}

	cpp_int find_period(cpp_int a, cpp_int N) {
		cpp_int aa = a;
		for (cpp_int x = 1; x < N; ++x) {
			aa = aa % N;
			if (aa == 1) {
				return x;
			}
			aa *= a;
		}
		return 0;
	}

	PrimeFactorization lcm(PrimeFactorization f1, const PrimeFactorization& f2) {
		for (const auto& g: f2) {
			auto factor = g.first;
			auto exponent2 = g.second;
			if (f1.find(factor) == f1.end()) {
				f1[factor] = exponent2;
			} else {
				auto exponent1 = f1.at(factor);
				f1[factor] = max(exponent1, exponent2);
			}
		}
		return f1;
	}

	PrimeFactorization lcm(const vector<PrimeFactorization>& fs) {
		PrimeFactorization f(fs.front());
		for (const PrimeFactorization& f2: fs) {
			f = lcm(f, f2);
		}
		return f;
	}

	PrimeFactorization gcd(PrimeFactorization f1, const PrimeFactorization& f2) {
		for (const auto& g: f2) {
			auto factor = g.first;
			auto exponent2 = g.second;
			if (f1.find(factor) != f1.end()) {
				auto exponent1 = f1.at(factor);
				f1[factor] = min(exponent1, exponent2);
				if (f1[factor] == 0) { f1.erase(factor); }
			}
		}
		vector<cpp_int> deletables;
		for (const auto& g: f1) {
			auto factor = g.first;
			if (f2.find(factor) == f2.end()) {
				deletables.emplace_back(factor);
			}
		}
		for (const auto& x : deletables) {
			f1.erase(x);
		}
		return f1;
	}

	PrimeFactorization gcd(const vector<PrimeFactorization>& fs) {
		PrimeFactorization g = fs.front();
		for (const auto& f : fs) {
			g = gcd(f, g);
		}
		return g;
	}

	PrimeFactorization lambda(const PrimeFactorization& f) {
		vector<PrimeFactorization> rf;
		for (const auto& pp: f) {
			auto p = pp.first;
			auto e = pp.second;
			if (e > 1) {
				cout << "not sure if this is still true...\n";
				exit(1);
			}
			rf.push_back(factorize(p - 1));
		}
		return lcm(rf);
	}

	cpp_int product(const PrimeFactorization& f) {
		cpp_int N = 1;
		for (const auto& pp: f) {
			N *= pow(pp.first, pp.second);
		}
		return N;
	}

	PrimeFactorization operator/(PrimeFactorization Q, const PrimeFactorization& d) {
		/// return Q / d
		for (auto pp: d) {
			cpp_int p = pp.first;
			size_t e = pp.second;
			Q[p] -= e;
			if (Q[p] == 0) { Q.erase(p); }
		}
		return Q;
	}

	cpp_int number_of_divisors(const PrimeFactorization& f) {
		/// returns the number of divisors d of f, d = prod_k e_k
		cpp_int n = 1;
		for (const auto& pp: f) {
			n *= pp.second + 1;
		}
		return n;
	}

	/// create the i-th divisor of f
	PrimeFactorization divisor(const PrimeFactorization& f, cpp_int I) {
		PrimeFactorization d;
		for (const auto& pp: f) {
			const auto& p = pp.first;
			const auto& e = pp.second;
			auto de = I % (e + 1);
			if (de > 0) {
				d[p] = (size_t) de;
			}
			I /= (e + 1);
		}
		return d;
	}

	cpp_int mu(const PrimeFactorization& f) {
		if (!is_square_free(f)) {
			return 0;
		} else {
			size_t k = f.size();
			return (cpp_int) pow(-1, k);
		}
	}

	cpp_int number_of_elements_with_order(const vector<PrimeFactorization>& G, const PrimeFactorization& m) {
		cpp_int D = number_of_divisors(m);
		cpp_int res = 0;
		for (size_t I = 0; I < D; ++I) {
			PrimeFactorization l = divisor(m, I);
//			cout << "l = " << product(l) << endl;
//			cout << "m / l = " << m / l << endl;
			cpp_int muu = mu(m / l);
//			cout << "mu(m / l) = " << muu << endl;
//			cout << " + (" << muu ;
			if (muu != 0) {
				cpp_int x = muu;
				for (const PrimeFactorization& g: G) {
//					cout << " * ";
//					cout << "gcd(" << product(l) << ", " << product(g) << ") = ";
					PrimeFactorization xx = gcd(l, g);
					x *= product(xx);
//					cout << product(xx);
				}
				res += x;
			}
//			cout << ")";
		}
//		cout << " = " << res << endl;
		return res;
	}

	vector<PrimeFactorization> factorizePeriods(const PrimeFactorization& f) {
		vector<PrimeFactorization> rfs;
		for (const auto& p: f) {
			PrimeFactorization rf = number_theory::factorize(p.first - 1);
			rfs.push_back(rf);
		}
		return rfs;
	}
}

std::ostream& operator<<(std::ostream& os, const number_theory::PrimeFactorization& f) {
	bool beg = true;
	if (f.empty()) {
		os << "1";
	}
	for (const auto& x: f) {
		if (!beg) {
			os << " * ";
		}
		beg = false;
		os << x.first << "^" << x.second;
	}
	return os;
}
