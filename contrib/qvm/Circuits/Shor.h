//
// Created by Roman Ellerbrock on 8/31/20.
//

#ifndef SHOR_H
#define SHOR_H
#include "Arithmetic.h"
#include "QuantumCircuits.h"
#include "../InTimeQuantumInstruction.h"
#include "../ConditionalQuantumInstruction.h"
#include "../OutputInstruction.h"

namespace Shor {

	pair<vector<size_t>, vector<size_t>> denominator(double x, size_t N);

	vector<size_t> check_criteria(vector<size_t> s, vector<size_t> d,
		size_t N, double theta, size_t Q);

	/* Iterative Function to calculate (x^y)%p in O(log y) */
	int power(int x, unsigned int y, int p);

	cpp_int pow_mod(cpp_int a, size_t r, cpp_int N);

	void check(size_t a, size_t N, size_t m, const Register& x);

	cpp_int find_period(cpp_int a, cpp_int N);

	void factorize(cpp_int a, cpp_int r, cpp_int N);

	void factorize(cpp_int a, cpp_int N);

	vector<size_t> create_table(size_t a, size_t N, size_t Q);

	vector<size_t> collect_z(const vector<size_t> f);

	vector<complex<double>> phases(const vector<size_t> f,
		const size_t z, const cpp_int N, size_t Q);

	void predict_phases(const cpp_int a, const cpp_int N, const size_t q);
}

#endif //SHOR_H
