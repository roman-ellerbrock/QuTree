#ifndef GATE_ARITHMETIC_H
#define GATE_ARITHMETIC_H

#include "SOPVector.h"
#include "GateOperators.h"
#include "long_integer.h"
#include "Register.h"

namespace GateArithmetic {

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

	SOPVector const_add(const Register& b, const long_integer& a,
		bool adjungate, size_t approx = 0);

	SOPVector c_const_add(const Register& b, const long_integer& a,
		size_t control, bool adjungate, size_t approx);

	SOPVector cc_const_add(const Register& b, const long_integer& a,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	SOPVector const_subst(const Register& b, const long_integer& a,
		bool adjungate, size_t approx);

	SOPVector cc_add_mod(const Register& b, size_t zero,
		const long_integer& a, const long_integer& N,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	SOPVector c_mult_mod(const Register& b, const Register& x,
		const long_integer& a, const long_integer& N,
		size_t mode_0, size_t control, bool adjungate, size_t approx);

	SOPVector Ua(const Register& x, const Register& b,
		size_t mode_0, size_t control,
		const long_integer& a, const long_integer& N,
		size_t adjungate, size_t approx);

	SOPVector Shors(const Register& M, const Register& x, const Register& b,
		size_t mode_0, long_integer a, const long_integer& N,
		size_t adjungate, size_t approx);

	SOPVector c_phi_mac(const Register& b, const Register& x,
		const long_integer& a,
		bool adjungate, size_t approx, size_t control);

	SOPVector add_qft(const Register& a, const Register& b,
		bool adjungate, size_t approx = 0);

	SOPVector Set_Number(const Register& x, const long_integer& a);

	/// ==================================================================================

	SOPVector Random_Number(size_t mode_a, size_t n_bit);

	SOPVector Set_Number(const long_integer& a, size_t mode_a, size_t n_bit);

	SOPVector add_qft(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx = 0);

	SOPVector const_add(const long_integer& a, size_t n_bit, size_t mode_b, bool adjungate, size_t approx = 0);

	SOPVector const_subst(const long_integer& bin_a, size_t n_bit, size_t mode_b, bool adjungate, size_t approx);

	SOPVector cc_const_add(const long_integer& a, size_t n_bit, size_t mode_b,
		bool adjungate, size_t approx, size_t control1, size_t control2);

	SOPVector c_phi_mac(const long_integer& a, size_t n_bit, size_t mode_x,
		size_t mode_b, bool adjungate, size_t approx, size_t control);

	SOPVector cc_add_mod(const long_integer& bin_a, const long_integer& bin_N, size_t n_bit,
		size_t mode_b, size_t mode_0, bool adjungate, size_t approx,
		size_t control1, size_t control2);

	SOPVector c_mult_mod(const long_integer& a, const long_integer& N,
		size_t n_bit, size_t mode_x, size_t mode_b, size_t mode_0,
		bool adjungate, size_t approx, size_t control);

	SOPVector Ua(const long_integer& a, const long_integer& a_inv, const long_integer& N,
		size_t n_bit, size_t mode_x, size_t mode_b, size_t mode_0,
		bool adjungate, size_t approx, size_t control);

	SOPVector c_swap(size_t control, size_t target1, size_t target2);

	SOPVector exp_mod(const long_integer& a, const long_integer& N,
		const mctdhBasis& basis, bool adjungate, size_t approx);

	/**
	 *
	 */
	vector<size_t> binary_rep(size_t a, size_t n_bit);

	vector<size_t> padded(vector<size_t> a, size_t n_bit);

	vector<size_t> bitshift(const vector<size_t>& bit_a, size_t x);

	long double compute_x(const long_integer& a, size_t j);
}

#endif //GATE_ARITHMETIC_H
