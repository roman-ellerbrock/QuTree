//
// Created by Roman Ellerbrock on 9/22/22.
//

#ifndef NUMBERTHEORY_H
#define NUMBERTHEORY_H
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <boost/math/special_functions/prime.hpp>
#include <map>

namespace number_theory {
	using namespace boost::multiprecision;
	using namespace std;

	using PrimeFactorization = map<cpp_int, size_t>; /// which factor & how often

	PrimeFactorization factorize(cpp_int N);
	bool is_square_free(const PrimeFactorization& f);

	vector<cpp_int> coprimes(const cpp_int& N);
	PrimeFactorization lcm(const vector<PrimeFactorization>& fs);
	PrimeFactorization gcd(const vector<PrimeFactorization>& fs);
	PrimeFactorization lambda(const PrimeFactorization& f);
	cpp_int product(const PrimeFactorization& f);

	cpp_int number_of_divisors(const PrimeFactorization& f);
	PrimeFactorization divisor(const PrimeFactorization& f, cpp_int I);
	vector<PrimeFactorization> factorizePeriods(const PrimeFactorization& f);

	cpp_int number_of_elements_with_order(const vector<PrimeFactorization>& G, const PrimeFactorization& m);
}

std::ostream& operator<<(std::ostream& os, const number_theory::PrimeFactorization & f);

#endif //NUMBERTHEORY_H
