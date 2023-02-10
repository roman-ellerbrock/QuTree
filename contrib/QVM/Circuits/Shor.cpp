//
// Created by Roman Ellerbrock on 8/31/20.
//

#include "Shor.h"
#include <boost/multiprecision/cpp_int.hpp>
using namespace boost::multiprecision;

namespace Shor {

	pair<vector<size_t>, vector<size_t>> denominator(double x, size_t N)  {
		vector<size_t> seq;
		if (x < 1) {
			x = 1. / x;
		}
		size_t dim = 10;
		size_t k = dim;
		for (size_t i = 1; i <= dim; ++i) {
			seq.push_back((size_t) x);
			if (abs(x - (size_t) x) < 1e-9) {
				k = i + 1;
				break;
			}
			x = 1. / (x - (size_t) x);
		}

		vector<size_t> ans, ans2;
		for (size_t i = 0; i < k; ++i) {
			size_t up = 1;
			size_t down = seq[i];
			for (size_t j = 0; j < i; ++j) {
				size_t t = seq[i - 1 - j];
				size_t d = down;
				down = down * t + up;
				up = d;
			}
			ans.push_back((size_t) down);
			ans2.push_back((size_t) up);
		}
		return {ans, ans2};
	}

	int power(int x, unsigned int y, int p) {
		int res = 1;     // Initialize result

		x = x % p; // Update x if it is more than or
		// equal to p

		if (x == 0) return 0; // In case x is divisible by p;

		while (y > 0) {
			// If y is odd, multiply x with result
			if (y & 1) {
				res = (res * x) % p;
			}

			// y must be even now
			y = y >> 1; // y = y/2
			x = (x * x) % p;
		}
		return res;
	}

	cpp_int pow_mod(cpp_int a, size_t r, cpp_int N) {
		cpp_int x = 1;
		for (size_t i = 0; i < r; ++i) {
			x *= a;
			x = x % N;
		}
		return x;
	}

	void check(size_t a, size_t N, size_t m, const Register& x) {

		size_t Q = pow(2, 2 * x.size());
		double theta = (double) m / (double) Q;
		auto xs = Shor::denominator(theta, N);
		auto s = xs.first;
		auto d = xs.second;
		auto rs = Shor::check_criteria(s, d, N, theta, Q);

		cout << "a: " << a << endl;
		cout << "N: " << N << endl;
		cout << "m: " << m << endl;
		cout << "Q: " << Q << endl;
		cout << "pruned s:" << endl;
		cout << "theta: " << theta << endl;
		cout << "pruned estimate s: [";
		for (auto x: rs) {
			cout << x << " ";
		}
		cout << "]" << endl;

		cpp_int aa(a);
		cpp_int NN(N);
		for (size_t mult = 1; mult < N; ++mult) {
			vector<size_t> r2 = rs;
			for (auto& r: r2) {
				r *= mult;
			}

			for (size_t i: r2) {
				if (i < N) {
					size_t mm = power(a, i, N);
					if (mm == 1) {
						cout << "found period: r = " << i << endl;
						if (i % 2 == 0) {
							auto zz = pow_mod(aa, i / 2, NN);
							//auto z = zz.convert();
							cout << "a^(r/2) % N = " << zz << endl;
							auto f1 = gcd(zz + 1, N);
							auto f2 = gcd(zz - 1, N);
							cout << "Factors: " << f1 << ", " << f2 << endl;
						}
						break;
					}
				}
			}
		}
	}

	cpp_int find_period(cpp_int a, cpp_int N) {
		cpp_int aa = a;
		for (cpp_int x = 1; x < N; ++x) {
//			cout << x << "\t" << aa << endl;
			aa = aa % N;
			if (aa == 1) {
				return x;
			}
//			if (x % 10000 == 0) {
//				cout << x << endl;
//			}
			aa *= a;
		}
//		cerr << "could not find period.\n";
//		exit(1);
		return 0;
	}

	void factorize(cpp_int a, cpp_int r, cpp_int N) {
		cout << "r = " << r << ", ";
		if (r % 2) {
			cout << "period not even. Try again.\n";
		} else {
			cpp_int t1 = pow_mod(a, (size_t) r / 2, N) + 1;
//			cout << "t1: " << t1 << endl;
//			cpp_int t1 = pow(a, (size_t) r / 2) + 1;
			auto f1 = boost::math::gcd(t1, N);
			cpp_int t2 = pow_mod(a, (size_t) r / 2, N) - 1;
//			cpp_int t2 = pow(a, (size_t) r / 2) - 1;
//			cout << "t2: " << t2 << endl;
			auto f2 = boost::math::gcd(t2, N);
			cout << f1 << " " << f2 << endl;
		}
	}

	void factorize(cpp_int a, cpp_int N) {
		cpp_int r = Shor::find_period(a, N);
		factorize(a, r, N);
	}

	vector<size_t> create_table(size_t a, size_t N, size_t Q) {
		vector<size_t> f(Q);
		cpp_int t = 1;
		for (size_t x = 0; x < Q; ++x) {
			f[x] = (size_t) t;
			t *= (cpp_int) a;
			t = t % N;
		}
		return f;
	}

	vector<size_t> collect_z(const vector<size_t> f) {
		vector<size_t> z;
		for (size_t i = 1; i < f.size(); ++i) {
			z.push_back(f[i]);
			if (f[i] == 1) { break; }
		}
		return z;
	}

	vector<complex<double>> phases(const vector<size_t> f,
		const size_t z, const cpp_int N, size_t Q) {

		// number of qubits
		auto two_pi_i = QM::two_pi * QM::im / (double) Q;
		vector<complex<double>> phase(f.size());
		for (size_t y = 0; y < f.size(); ++y) {
//			cout << y << " ";
			complex<double> ph = 0;
			for (size_t x = 0; x < f.size(); ++x) {
				if (f[x] == z) {
					auto add = exp(two_pi_i * (double) (x * y));
					ph += add;
				}
			}
			phase[y] = ph / (double) Q;
		}
		return phase;
	}

	void predict_phases(const cpp_int a, const cpp_int N, const size_t q) {
		size_t Q = pow(2, 2 * q);
		auto f = Shor::create_table((size_t) a, (size_t) N, Q);
		for (size_t x = 0; x < f.size(); ++x) {
			cout << x << " - " << f[x] << endl;
		}
		auto z = Shor::collect_z(f);
		for (auto y: z) {
			cout << y << " ";
		}
		cout << endl;
		cout << "period: " << z.size();
//		for (size_t m = 0; m < z.size(); ++m) { /// all m will give the same distribution
		{
			size_t m = 0;
			auto phases = Shor::phases(f, z[m], N, Q);
			ofstream os(to_string((size_t) N) + "." + to_string((size_t) a) + "." + to_string((size_t) f[m]) + ".dat");
			for (size_t i = 0; i < f.size(); ++i) {
				os << i << " " << abs(phases[i]) << endl;
			}
		}
		getchar();
	}

	vector<size_t> check_criteria(vector<size_t> s, vector<size_t> d,
		size_t N, double theta, size_t Q) {

		vector<size_t> d_over;
		for (size_t i = 0; i < s.size(); ++i) {
			d_over.push_back(d[i] / s[i]);
		}

		vector<size_t> r;
		for (size_t i = 0; i < s.size(); ++i) {
			if ((s[i] < N) && (abs(theta - d_over[i]) < (1. / 2. * Q))) {
				r.push_back((s[i]));
			}
		}
		return r;
	}

	void gcdx(cpp_int& g, cpp_int& x, cpp_int& y, cpp_int a, cpp_int b) {
		cpp_int x0 = 0;
		cpp_int x1 = 1;
		cpp_int y0 = 1;
		cpp_int y1 = 0;
		while (a != 0) {
			cpp_int q = b / a;
			cpp_int tmp = a;
			a = b % a;
			b = tmp;
			tmp = y0;
			y0 = y1;
			y1 = tmp - q * y1;
			tmp = x0;
			x0 = x1;
			x1 = tmp - q * x1;
		}
		g = b;
		x = x0;
		y = y0;
	}

	long_integer LargeEuclidicInverse(const long_integer& a, const long_integer& N) {
		assert(a.size() <= 32);
		assert(N.size() <= 32);
		size_t aa = a.convert();
		size_t NN = N.convert();
		cpp_int g = 0;
		cpp_int x = 0;
		cpp_int y = 0;
		gcdx(g, x, y, aa, NN);
		cpp_int max = 0;
		max -= NN;
		if (x > max) { x += NN; }
		if (y > max) { y += NN; }
//	cout << "g = " << g << endl;
//	cout << "x_ = " << x << endl;
//	cout << "y = " << y << endl;
		if (g != 1) {
			cerr << "gcd(a, N) is not 1!\nPlease choose valid numbers.\n";
		}
		size_t z = (size_t) x;
		return long_integer(z, a.size());
	}
}