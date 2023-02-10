//
// Created by Roman Ellerbrock on 9/1/22.
//

#include "ShorSparsetensor.h"
#include <boost/math/special_functions/prime.hpp>

complex<double> Rphase(const vector<bool>& m) {
	double phi = 0;
	size_t j = m.size() + 1;
/*	cout << "measurements:\n";
	for (auto x : m) {
		cout << x << ", ";
	}
	cout << endl;
*/
//	cout << "phases:\n";
	for (int k = 2; k <= j; ++k) {
		auto add = m[j - k] / pow(2., k);
//		cout << k << " " << m[j - k] << ", phase: " << add << endl;
		phi += add;
	}
	return exp(-QM::two_pi * QM::im * phi);
}

Matrixcd X() {
	Matrixcd x(2, 2);
	x(1, 0) = 1;
	x(0, 1) = 1;
	return x;
}

void runShorSparsetensorVerbose(size_t a, size_t N, size_t n, mt19937& gen) {
	ShorSparsetensor state(a, N, n, gen);

	state.print();
	for (int L = 2 * n - 1; L >= 0; --L) {
		cout << "Start of cycle: " << L << endl;

		cout << "Apply H:\n";
		state.applyH();
		state.print();
		getchar();

		size_t aa = (size_t) Uaa(a, L, n, N);
		cout << L << " - " << aa << endl;
		cout << "Apply Ua(" << aa << ") for L = " << L << "/" << 2 * n << ":\n";
		getchar();
		state.applyUa(aa);
		state.print();
		getchar();

		cout << "Apply Phase + H:\n";
		state.applyPhaseAndH();
		state.print();
		getchar();

		cout << "measure & X:\n";
		state.measureAndX();
		state.print();
		getchar();
	}

	state.printResult();
}

cpp_int Uaa(size_t a, size_t L, size_t n, cpp_int N) {
	cpp_int aa = a % N;
	for (size_t l = 0; l < L; ++l) {
//			cout << l << " - " << aa << endl;
		aa = (aa * aa) % N;
	}
	return aa;
}

void runShorOnCoPrimes(size_t m, size_t n, mt19937& gen, size_t nsample) {
	cpp_int p = boost::math::prime(m);
	cpp_int q = boost::math::prime(n);
	cpp_int N = p * q;
	cout << "N = " << N;
	cout << " = " << p << " * " << q << " | (" << m << ", " << n << ")";
	/// search l by increasing 2^l until N fits into l bits.
	size_t l = 5;
	while ((N - pow((cpp_int) 2, l)) > 0) {
		l++;
		if (l > 100) {
			cerr << "Cannot find power of 2.\n";
			exit(1);
		}
	}
	cout << ", l = " << l;

	/// run Shor's
	cpp_int a = findA(N, gen);
	ShorSparsetensor state((size_t) a, (size_t) N, l, gen);
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (size_t i = 0; i < nsample; ++i) {
		cout << i;
		cout.flush();
//		while (!state.reinitA()) { }
		state.sample();
//		cout << ", a = " << state.a_ << "\t";
//		state.printResult();
//		runShorSparsetensor(state, (size_t)a, (size_t)N, l, gen);
		state.clear();
	}
	end = std::chrono::system_clock::now();
	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000000.;
	auto timePerBitstring = time / (double) nsample;
	auto bitstringsPerSecond = 1. / timePerBitstring;
	cout << "s(r) = " << state.s_r_ << ", s(factor(r)) = " << state.s_factor_ << ", s_total = " << state.s_tot_ << endl;
	cout << "nSample: " << nsample << ", Time: " << time << "s, time/bitstring: " << timePerBitstring
		 << "s, bitstrings/s: " << bitstringsPerSecond << endl;
}

void runShorOnCoPrimes(size_t n, mt19937& gen, size_t nsample) {
	ShorSparsetensor state(gen);
	size_t range = floor((double) n / 4.);
	uniform_int_distribution<size_t> dist(n - range, n + range);
	cout << "range(n) : [" << n - range << ", " << n + range << "]\n";
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (size_t M = 0; M < nsample; ++M) {
		cout << M << ", ";
		size_t n_trys = 0;
		size_t m = 0;
		n = 0;
		while (m == n) {
			m = dist(gen);
			n = dist(gen);
			n_trys++;
			if (n_trys > 1000) {
				cerr << "Finding prime indices failed. Only found squares in 1000 trys.\n";
				exit(1);
			}
		}

		cpp_int p = boost::math::prime(m);
		cpp_int q = boost::math::prime(n);
		cpp_int N = p * q;
		cout << "N = " << N;
		cout << " = " << p << " * " << q << " | (" << m << ", " << n << ")";
		/// search l by increasing 2^l until N fits into l bits.
		size_t l = 5;
		while ((N - pow((cpp_int) 2, l)) > 0) {
			l++;
			if (l > 100) {
				cerr << "Cannot find power of 2.\n";
				exit(1);
			}
		}
		cout << ", l = " << l;
		cout.flush();

		cpp_int a = findA(N, gen);
		cout << ", a = " << a << " ";
		cout.flush();

		state.initialize((size_t) N, (size_t) a, l);
		state.sample();
		state.printResult();
		state.clear();
	}
	end = std::chrono::system_clock::now();
	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count() / 1000000.;
	auto timePerBitstring = time / (double) nsample;
	auto bitstringsPerSecond = 1. / timePerBitstring;
	cout << "s(r) = " << state.s_r_ << ", s(factor(r)) = " << state.s_factor_ << ", s_total = " << state.s_tot_
		 << ", avg(n) = " << state.sumn_ / (double) nsample << endl;
	cout << "nSample: " << nsample << ", Time: " << time << "s, time/bitstring: " << timePerBitstring
		 << "s, bitstrings/s: " << bitstringsPerSecond << endl;
}

void ShorSparsetensor::sample() {
//	size_t M = (size_t) Ln(0.75, 1., n_);
//	cout << "dim = " << (size_t) sqrt((double) ys_.shape_.totalDimension()) << endl;
//	cout << endl;
	for (int L = 2 * n_ - 1; L >= 0; --L) {
		applyH();
		cpp_int aac = Uaa(a_, L, n_, N_);
		size_t aa = (size_t) aac;
		applyUa(aa);
//		cout << L << " ";
//		checkRank();
//		prune(M);
		applyPhaseAndH();
		measureAndX();
	}
//	cout << ", M = " << M << "\t";
	evaluateStatistics();
}

cpp_int findA(const cpp_int& N, mt19937& gen) {
	uniform_int_distribution<cpp_int> dist((cpp_int) 2, N - 1);
	cpp_int r = 0;
	cpp_int a = 0;
	while (r != 1) {
		a = dist(gen);
		r = gcd(a, N);
	}
	return a;
}

cpp_int convert(const vector<bool>& m) {
	cpp_int n = 0;
	cpp_int two = 2;

	for (size_t i = 0; i < m.size(); ++i) {
		n += ((size_t) m[i]) * pow(two, i);
	}
	return n;
}

cpp_int convertBack(const vector<bool>& m) {
	cpp_int n = 0;
	cpp_int two = 2;

	for (size_t i = 0; i < m.size(); ++i) {
		size_t p = m.size() - 1 - i;
		n += ((size_t) m[i]) * pow(two, p);
	}
	return n;
}

vector<int> continuousFraction(double r, size_t nmax) {
	vector<int> a;
	for (size_t i = 0; i < nmax; ++i) {
		int j = floor(r + 1e-10);
		r -= j;
		a.push_back(j);
		if (abs(r) < 1e-8) { break; }
		r = 1. / r;
	}
	return a;
}

pair<size_t, size_t> fractionApproximation(const vector<int>& a, size_t n) {
	size_t h1 = 1;
	size_t h2 = 0;
	size_t k1 = 0;
	size_t k2 = 1;
	for (size_t i = 0; i < n; ++i) {
		size_t hi = a[i] * h1 + h2;
		size_t ki = a[i] * k1 + k2;
		h2 = h1;
		k2 = k1;
		h1 = hi;
		k1 = ki;
	}
	return {h1, k1};
}

cpp_int Ln(double alpha, size_t c, size_t nin) {
	double n = log(pow(2., nin));
	return (cpp_int) exp(((double) c) * pow(n, alpha) * pow(log(1. * n), 1. - alpha));
}

vector<cpp_int> screen_a(const cpp_int& N, mt19937& gen, size_t M) {
	vector<cpp_int> r;
	for (size_t i = 0; i < M; ++i) {
		cpp_int a = findA(N, gen);
		r.push_back(Shor::find_period(a, N));
	}

	map<cpp_int, size_t> hist;
	for (const auto& c: r) {
		if (hist.find(c) != hist.end()) {
			hist[c]++;
		} else {
			hist[c] = 1;
		}
	}
	for (const auto& p: hist) {
		cout << p.first << "\t" << (double) p.second / (double) r.size() << endl;
	}
	return r;
}

void screen_a2(size_t n, mt19937& gen, size_t nsample) {

	vector<cpp_int> r;

	size_t range = floor((double) n / 4.);
	uniform_int_distribution<size_t> dist(n - range, n + range);
	cout << "range(n) : [" << n - range << ", " << n + range << "]\n";
	std::chrono::time_point<std::chrono::system_clock> start, end;
	start = std::chrono::system_clock::now();
	for (size_t M = 0; M < nsample; ++M) {
		size_t n_trys = 0;
		size_t m = 0;
		n = 0;
		while (m == n) {
			m = dist(gen);
			n = dist(gen);
			n_trys++;
			if (n_trys > 1000) {
				cerr << "Finding prime indices failed. Only found squares in 1000 trys.\n";
				exit(1);
			}
		}

		cpp_int p = boost::math::prime(m);
		cpp_int q = boost::math::prime(n);
		cpp_int N = p * q;

		cpp_int a = findA(N, gen);
		r.push_back(Shor::find_period(a, N));


	}
}
