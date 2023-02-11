//
// Created by Roman Ellerbrock on 9/1/22.
//

#ifndef SHORSPARSETENSOR_H
#define SHORSPARSETENSOR_H
#include "Shor.h"
#include <random>

Tensorcd createOne();
Tensorcd createZero();

complex<double> Rphase(const vector<bool>& m);
Matrixcd X();

cpp_int convert(const vector<bool>& m);
cpp_int convertBack(const vector<bool>& m);
cpp_int Uaa(size_t a, size_t L, size_t n, cpp_int N);
vector<int> continuousFraction(double r, size_t nmax);
pair<size_t, size_t> fractionApproximation(const vector<int>& a, size_t n);
cpp_int findA(const cpp_int& N, mt19937& gen);
/// note: n is the number of qubits, not the number N that is factorized
cpp_int Ln(double alpha, size_t c, size_t n);

class ShorSparsetensor {
public:

	ShorSparsetensor(mt19937& gen)
		: gen_(gen) {}

	void initialize(size_t a, size_t N, size_t n) {
		a_ = a;
		N_ = N;
		n_ = n;

		H_ = Matrixcd(2, 2);
		H_(0, 0) = 1. / sqrt(2.);
		H_(1, 0) = 1. / sqrt(2.);
		H_(0, 1) = 1. / sqrt(2.);
		H_(1, 1) = -1. / sqrt(2.);

		r_ = (size_t) Shor::find_period(a, N_);
//		cout << "\nPeriod: " << r_ << endl;
		table_ = Shor::create_table(a, N_, N_);

		ys_ = Tensorcd({2, r_});
		ys_(0, 0) = 1.;
		y1_ = Tensorcd(ys_.shape_);
	}

	ShorSparsetensor(size_t a, size_t N, size_t n, mt19937& gen)
		: a_(a), N_(N), n_(n), gen_(gen), size_(1) {

		initialize(a, N, n);
		cout << "\nPeriod: " << r_ << endl;
	}

	cpp_int reinitA() {
		a_ = (size_t) findA(N_, gen_);
		r_ = (size_t) Shor::find_period(a_, N_);
		table_ = Shor::create_table(a_, N_, N_);
		return r_;
	}

	void clear() {
		ys_.zero();
		ys_(0, 0) = 1.;
		y1_.zero();
		size_ = 1;
		measurements_.clear();
	}

	void applyH() {
		ys_ = matrixTensor(H_, ys_, 0);
	}

	void applyUa(size_t aa) {

		if (aa == 1) { return; }
		size_t shift = find_shift(aa);

		Tensorcd& y0 = ys_;
		y1_.zero();
		size_t y;
		for (size_t i = 0; i < r_; ++i) {
			y = (i + shift) % r_;
			y1_(1, y) = y0(1, i);
			y0(1, i) = 0;
		}
		y0 += y1_;
	}

	bool elementIsOccupied(const Tensorcd& y0, size_t j, double eps) const {
		double x = pow(abs(y0(0, j)), 2) + pow(abs(y0(1, j)), 2);
		if (x > eps) {
			return true;
		} else {
			return false;
		}
	}

	void prune(size_t maxsize) {

		if (maxsize == r_) { return; }
		if (maxsize <= 0) { maxsize = 1; }
		double eps = 1e-12;

		size_ = 0;
		for (size_t i = 0; i < r_; ++i) {
			if (elementIsOccupied(ys_, i, eps)) {
				size_++;
			}
		}

		uniform_int_distribution<size_t> dist(0, r_ - 1);
		while (size_ > maxsize) {
			while (true) {
				size_t x = dist(gen_);
				if (elementIsOccupied(ys_, x, eps)) {
					ys_(0, x) = 0.;
					ys_(1, x) = 0.;
					size_--;
					break;
				}
			}
		}
	}

	void applyPhaseAndH() {
		auto phase = Rphase(measurements_);
		for (size_t i = 0; i < r_; ++i) {
			ys_(1, i) *= phase;
		}
		ys_ = matrixTensor(H_, ys_, 0);
	}

	void measureAndX() {
		double p0 = 0.;
		for (size_t i = 0; i < r_; ++i) {
			p0 += pow(abs(ys_(0, i)), 2);
		}

//		cout << "p(|0>) = " << p0 << endl;
		double r = dist_(gen_);
		size_t act = 0;
		size_t off = 1;
		if (r > p0) {
			act = 1;
			off = 0;
		}
		measurements_.push_back(act);

		/// Apply X & re-normalize state
		double norm = 0.;
		for (size_t i = 0; i < r_; ++i) {
			ys_(0, i) = ys_(act, i);
			ys_(1, i) = 0.;
			norm += pow(abs(ys_(0, i)), 2);
		}

		norm = 1. / sqrt(norm);
		ys_ *= norm;
	}

	void applyX() {
		Matrixcd x = X();
		ys_ = matrixTensor(x, ys_, 0);
	}

	size_t find_shift(cpp_int aa) {
		size_t shift = 0;
		bool found = false;
		for (size_t i = 0; i < table_.size(); ++i) {
			if (table_[i] == aa) {
				shift = i;
				found = true;
				break;
			}
		}
		if (!found) {
			cerr << "Did not find number aa = " << aa << " in table with " << table_.size() << " elements.\n";
			cout << "r = " << r_ << endl;
			cout << "N_ = " << N_ << endl;
			exit(1);
		}
		return shift;
	}

	void print() const {
		cout << "State vector:\n";
		for (size_t i = 0; i < r_; ++i) {
			cout << table_[i] << "\t(" << ys_(0, i) << ", " << ys_(1, i) << ")\n";
		}
	}

	void evaluateStatistics() {
		auto num = convert(measurements_);
		auto a = continuousFraction((double) num / pow(2., 2 * n_), n_);
		for (size_t i = 0; i <= a.size(); ++i) {
			auto hk = fractionApproximation(a, i);
			auto s = hk.second;
			if (s <= 1) { continue; }
			if (s > N_) { break; }
			if ((r_ % s) == 0 && (double) r_ / (double) s < log(r_)) { s_r_++; }
			if (r_ % s == 0 && s < r_) { s_factor_++; };
			s_tot_++;
		}
		sumn_ += n_;
	}

	void printResult() const {
//		cout << "Results of simulation: ";
/*		for (const bool m : measurements_) {
			cout << m << " ";
		}
		cout << " | " <<; */
		auto num = convert(measurements_);
		cout << num << "\t, r = " << r_ << "\t | s = {";
		auto a = continuousFraction((double) num / pow(2., 2 * n_), n_);
		size_t nnn = 0;
		for (size_t i = 0; i <= a.size(); ++i) {
			auto hk = fractionApproximation(a, i);
			auto s = hk.second;
			if (s <= 1) { continue; }
			if (s > N_) { break; }
			if (nnn > 0) { cout << ", "; }
			if (r_ % s == 0) { cout << "\x1b[32m"; } //else { cout << "\x1b[43m"; }
			cout << s;
			if (r_ % s == 0) { cout << "\x1b[0m"; }
			nnn++;
		}
		cout << "}";
		cout << endl;
	}

	void checkRank() const {
		size_t dimtot = ys_.shape_.totalDimension();
		size_t fac = sqrt((double) dimtot);
		size_t dim1 = dimtot / fac;
		size_t dim2 = round(dimtot / (double)dim1);
//		cout << dim1*dim2 - dimtot << endl;
//		cout << dim1 << " " << dim2 << endl;
		auto tmp = ys_;
		tmp.reshape({dimtot});
		auto y2 = tmp.adjustDimensions({dim1*dim2});
		y2.reshape({dim1, dim2});
		Matrixcd y2mat = toMatrix(y2);
//		cout << y2mat.dim1() << " " << y2mat.dim2() << endl;
//		y2mat.print();
		auto x = svd(y2mat);
		auto sigma = get<2>(x);
		for (size_t i = 0; i < sigma.dim(); ++i) {
			sigma(i) = pow(sigma(i), 2);
		}
		double norm = 0.;
		for (size_t i = 0; i < sigma.dim(); ++i) {
			norm += abs(sigma(i));
		}
		sigma /= norm;
		double perp = 0.;
		for (size_t i = 0; i < sigma.dim(); ++i) {
			perp -= log(abs(sigma(i))+1e-10) * sigma(i);
		}
		perp = exp(perp);
		cout << perp << endl;
//		sigma.print();
	}

	void sample();

	/// 2^31 - 1
	/// phase, measurement bit
	vector<bool> measurements_{};
	/// Replace map<cpp_int, Tensorcd> by the following two vectors:
	vector<size_t> table_;
	Tensorcd ys_;
	Tensorcd y1_;
	size_t a_;
	size_t N_;
	size_t n_;
	size_t r_;
	size_t size_{1};
	size_t sumn_{0};

	size_t s_factor_{0};
	size_t s_r_{0};
	size_t s_tot_{0};

	Matrixcd H_;

	mt19937& gen_;
	uniform_real_distribution<double> dist_{0., 1.};
};

void runShorSparsetensorVerbose(size_t a, size_t N, size_t n, mt19937& gen);
cpp_int findA(const cpp_int& N, mt19937& gen);
void runShorOnCoPrimes(size_t n, size_t m, mt19937& gen, size_t nsample);
void runShorOnCoPrimes(size_t n, mt19937& gen, size_t nsample);
vector<cpp_int> screen_a(const cpp_int& N, mt19937& gen, size_t M);

#endif //SHORSPARSETENSOR_H