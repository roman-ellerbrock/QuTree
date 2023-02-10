//
// Created by Roman Ellerbrock on 9/11/20.
//

#include "QuadraticSieve.h"

namespace Sieve {

	vector<cpp_int> createSieve(cpp_int N, cpp_int num) {
		cpp_int sq = (cpp_int) sqrt((double) N);

		vector<cpp_int> Vx((size_t)num);
		for (size_t i = 0; i < num; ++i) {
			Vx[i] = pow(1 + i + sq, 2) - N;
		}

		return Vx;
	}

	vector<size_t> factorBase(cpp_int num) {
		if (num >= 10000) {
			cerr << "lookup table for primes too small.\n";
			exit(1);
		}
		vector<size_t> p((size_t) num);
		for (size_t i = 0; i < num; ++i) {
			p[i] = boost::math::prime(i);
		}
		return p;
	}

	void quadraticSieve(cpp_int N) {
		size_t num = 100;
		auto Vx = createSieve(N, num);
		for (size_t i = 0; i < Vx.size(); ++i) {
			cout << i << "\t" << Vx[i] << endl;
		}
	}

}