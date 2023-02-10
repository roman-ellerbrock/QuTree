#ifndef GATE_ARITHMETIC_H
#define GATE_ARITHMETIC_H

#include "TreeOperators/SOPVector.h"
#include "GateOperators.h"
#include "../Util/long_integer.h"
#include "Register.h"
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
using namespace boost::multiprecision;


namespace Circuits {

/*!
 * \namespace GateArithmetic
 * \brief Here, fundamental arithmetic expressions for quantum computing are defined.
 *
 * The namespace features basic arithmetic operations like addition, multiplication, exponentiation.
 * For the circuits it is important to stick to a consistens endian convention for binary numbers.
 * This namespace uses little-endian convention and relies on a little-endian QFT.
 *
 * [1] "Massively Parallel ...Eleven Years Later", Computer Physics Communication (2019)
 * */
 	template <class Z>
	SOPVectorcd const_add(const Register& b, const Z& a,
		bool adjungate, size_t approx = 0);

	template <class Z>
	SOPVectorcd c_const_add(const Register& b, const Z& a,
		size_t control, bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd cc_const_add(const Register& b, const Z& a,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd const_subst(const Register& b, const Z& a,
		bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd cc_add_mod(const Register& b, size_t zero,
		const Z& a, const Z& N,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd c_mult_mod(const Register& b, const Register& x,
		Z a, const Z& N,
		size_t mode_0, size_t control, bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd Ua(const Register& x, const Register& b, size_t mode_0,
		size_t control, const Z& a, const Z& N,
		bool adjungate, size_t approx);

	template <class Z>
	SOPVectorcd cUa(Z a, const Register& control, const Register& x,
		const Register& b, const Register& zero, const Z& N,
		size_t l, bool adjungate, size_t approx);

	SOPVectorcd Rprime(const Register& t, const vector<size_t>& measurements, size_t j, bool adjoint);

	template <class Z>
	SOPVectorcd Shors(const Register& M, const Register& x, const Register& b,
		const Register& zero, Z a, const Z& N,
		bool adjungate, size_t approx);

	SOPVectorcd Set_Number(const Register& x, const long_integer& a);

	SOPVectorcd c_swap(size_t control, size_t target1, size_t target2);

	SOPVectorcd Random_Number(size_t mode_a, size_t n_bit);

	SOPVectorcd Set_Number(const long_integer& a, size_t mode_a, size_t n_bit);

	SOPVectorcd add_qft(const Register& a, const Register& b,
		bool adjungate, size_t approx = 0);

	/// ==================================================================================

	/**
	 *
	 */
	vector<size_t> binary_rep(size_t a, size_t n_bit);

	vector<size_t> padded(vector<size_t> a, size_t n_bit);

	vector<size_t> bitshift(const vector<size_t>& bit_a, size_t x);

	long double compute_x(const long_integer& a, size_t j);
}

#endif //GATE_ARITHMETIC_H
