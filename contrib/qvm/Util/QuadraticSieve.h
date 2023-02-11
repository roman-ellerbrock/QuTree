//
// Created by Roman Ellerbrock on 9/11/20.
//

#ifndef QUADRATICSIEVE_H
#define QUADRATICSIEVE_H
#include "Circuits/Arithmetic.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
#include <boost/math/special_functions/prime.hpp>
using namespace boost::multiprecision;

namespace Sieve {

	vector<cpp_int> createSieve(cpp_int N, cpp_int num);

	vector<size_t> factorBase(cpp_int num);

	void quadraticSieve(cpp_int N);

}

#endif //QUADRATICSIEVE_H
