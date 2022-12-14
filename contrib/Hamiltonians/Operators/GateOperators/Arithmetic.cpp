#include "Arithmetic.h"
#include "GateOperators.h"
#include "QFT.h"

namespace GateArithmetic {

	SOPVector add_qft(const Register& a, const Register& b, bool adjungate, size_t approx) {
		using namespace GateQFT;
		SOPVector stack = QFT(b, adjungate, approx);
		stack.append(ControlledRotation(a, b, adjungate, approx));
		stack.append(iQFT(b, adjungate, approx));
		return stack;
	}

	SOPVector const_add(const Register& b, const long_integer& a, bool adjungate, size_t approx) {
		using namespace GateOperators;
		using namespace GateQFT;
		MPO M;
		size_t n_bit = b.Size();
		for (int j = 0; j < n_bit; j++) {
			long double x = compute_x(a, j);
			auto R = make_shared<Rx>(x, adjungate);
			M.push_back(R, b.Last() - j);
		}
		SOP S(M, 1.0);
		SOPVector stack;
		stack.append(S);
		return stack;
	}

	SOPVector c_const_add(const Register& b, const long_integer& a, size_t control,
		bool adjungate, size_t approx) {
		SOPVector stack = const_add(b, a, adjungate, approx);
		return GateOperators::DistribSingleControl(stack, control);
	}

	SOPVector cc_const_add(const Register& b, const long_integer& a,
		size_t control1, size_t control2, bool adjungate, size_t approx) {
		SOPVector stack = const_add(b, a, adjungate, approx);
		return GateOperators::DistribDoubleControl(stack, control1, control2);
	}

	SOPVector const_subst(const Register& b, const long_integer& a, bool adjungate, size_t approx) {
		return const_add(b, a, !adjungate, approx);
	}

	SOPVector c_phi_mac(const Register& b, const Register& x, const long_integer& a,
		bool adjungate, size_t approx, size_t control) {

		SOPVector stack;
		for (size_t xi = 0; xi < x.Size(); ++xi) {
			size_t controlx = x.Begin() + xi;
			long_integer tmp(a);
			tmp.bitshift(xi);
			SOPVector ccadd = cc_const_add(b, tmp, control, controlx, adjungate, approx);
			stack.append(ccadd);
		}

		return stack;
	}

	SOPVector cc_add_mod(const Register& b, size_t zero,
		const long_integer& a, const long_integer& N,
		size_t control1, size_t control2, bool adjungate, size_t approx) {
		using namespace GateOperators;
		using namespace GateQFT;
		/// cc-Modular adder (b) -> ((b + a) % N) where a, N are classical numbers
		/// This design follows Ref. [1] Fig. 5. See paper for an explanation.
		/// Circuit expects a quint b in Fourier representation

		size_t largest_bit = b.LargestBit();
		SOPVector stack;
		stack.append(cc_const_add(b, a, control1, control2, adjungate, approx));

		stack.append(const_subst(b, N, adjungate, approx));

		stack.append(iQFT(b, adjungate, approx));
		stack.append(CNot(largest_bit, zero));
		stack.append(QFT(b, adjungate, approx));

		stack.append(c_const_add(b, N, zero, adjungate, approx));
		stack.append(cc_const_add(b, a, control1, control2, !adjungate, approx));
		stack.append(iQFT(b, adjungate, approx));
		stack.append(MPO(X, largest_bit));
		stack.append(CNot(largest_bit, zero));
		stack.append(MPO(X, largest_bit));
		stack.append(QFT(b, adjungate, approx));
		stack.append(cc_const_add(b, a, control1, control2, adjungate, approx));

		return stack;
	}

	SOPVector c_mult_mod(const Register& b, const Register& x,
		const long_integer& a, const long_integer& N,
		size_t mode_0, size_t c, bool adjungate, size_t approx) {
		using namespace GateQFT;
		/// cc-Modular adder (b) -> ((b + a*x) % N) where a, N are classical numbers
		/// This design follows Ref. [1] Fig. 6. See paper for an explanation.
		/// Circuit expects a quint b computational (not Fourier) representation

//		assert(x.Size() + 1 == b.Size());
		SOPVector stack;
		/// Transform to Fourier space
		stack.append(QFT(b, adjungate, approx));

		/// Perform sequence of modular adders for qubit in x
		for (int i = 0; i < x.Size(); ++i) {
			long_integer ashift(a);
			ashift.bitshift_mod(i, N);
			size_t c2 = x.SmallestBit() - i;
			stack.append(cc_add_mod(b, mode_0, ashift, N, c, c2, adjungate, approx));
		}

		/// Transform back to original space
		stack.append(iQFT(b, adjungate, approx));
		return stack;
	}

	SOPVector Ua(const Register& x, const Register& b,
		size_t mode_0, size_t control,
		const long_integer& a, const long_integer& N,
		size_t adjungate, size_t approx) {
		/// cc-Modular adder (x, 0) -> ((a*x) % N, 0) where a, N are classical numbers.
		/// This design follows Ref. [1] Fig. 7. See paper for an explanation.
		/// Circuit expects a quint x computational (not Fourier) representation.
		/// This is the main building block for Shor's algorithm.

		SOPVector stack;
		/// Do a regular c_mult_mod (b -> b + x*a % N)
		stack.append(c_mult_mod(b, x, a, N, mode_0, control, adjungate, approx));

		/// Do a n-qubit c_swap from one-qubit c_swap
		for (size_t i = 0; i < x.Size(); ++i) {
			size_t idx_b = b.Last() - i;
			size_t idx_x = x.Last() - i;
			stack.append(c_swap(control, idx_x, idx_b));
		}

		/// Run the adjoint of c_mult_mod (therefore reverse and adjungate)
		auto a_inv = EuclidicInverse(a, N);
		cout <<  "a_inv\n";
		a_inv.print();
		auto tmp = c_mult_mod(b, x, a_inv, N, mode_0, control, !adjungate, approx);
		reverse(tmp.begin(), tmp.end());
		stack.append(tmp);

		return stack;
	}

	SOPVector Shors(const Register& M, const Register& x, const Register& b,
		size_t mode_0, long_integer a, const long_integer& N,
		size_t adjungate, size_t approx) {
		/**
		 * Shor's algorithm
		 * @param M measurement register
		 * @param x measurement register
		 * @param b measurement register
		 * @param mode_0 measurement register
		 * @param a guess integer of Shor's
		 * @param N Integer that is to be factorized
		 *
		 * \brief This function performs the main part of Shor's algorithm a^M mod N
		 *
		 * The stack generated by this function contains the modular exponentiation,
		 * note that for the full circuit Hadamard gates have to be applied on every
		 * qubit in M to prepare the initial state and an inverse QFT has to be
		 * applied after performing the present circuit.
		 */

		SOPVector stack;

		for (size_t i = 0; i < M.Size(); ++i) {
			int control = M.Last() - i;
//			int control = M.Last();
			cout << control <<  endl;
			a.print();
			a.print_4byte();
			stack.append(Ua(x, b, mode_0, control, a, N, adjungate, approx));
			a = mult_mod(a, a, N);
//			if (i == 3) { break; }
		}

		return stack;
	}

	SOPVector Set_Number(const Register& x, const long_integer& a) {
		assert(x.Size() == a.size());
		MPO M;
		for (size_t k = 0; k < x.Size(); ++k) {
			if (a[k] == 1) {
				M.push_back(GateOperators::X, x.Last() - k);
			}
		}
		SOP S(M);
		SOPVector stack;
		stack.append(S);
		return stack;
	}

	/// ==================================================================================

	SOPVector Random_Number(size_t mode_a, size_t n_bit) {
//		srand(time(nullptr));
//		int seed = rand();
		int seed = 4661242;
		uniform_real_distribution<double> dist(0., 1.);
		mt19937 rand_gen(seed);
		MPO M;
		for (size_t k = 0; k < n_bit; ++k) {
			if (dist(rand_gen) > 0.5) {
				M.push_back(GateOperators::X, k + mode_a);
			}
		}
		SOP S(M);
		SOPVector stack;
		stack.append(S);
		return stack;
	}

	SOPVector Set_Number(const long_integer& bin_a, size_t mode_a, size_t n_bit) {

		MPO M;
		for (size_t k = 0; k < bin_a.size(); ++k) {
			if (bin_a[k] == 1) {
				M.push_back(GateOperators::X, mode_a + n_bit - 1 - k);
			}
		}
		SOP S(M);
		SOPVector stack;
		stack.append(S);
		return stack;
	}

	// Return x for single particle operator A_j
	// Converts a from big endian into little endian PLZ FIX
	long double compute_x(const long_integer& a, size_t j) {
		// Compute x
		long double x = 0.0;
		for (int k = 0; k <= j; k++) {
			x += a[j - k] / pow(2, k + 1);
//			cout << a[j - k] / pow(2, k + 1) << " ";
		}
//		cout << endl;
		return x;
	}

	SOPVector add_qft(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx) {
		using namespace GateQFT;

		cout << "Initializing QFT adder circuit ...\n";
		cout.flush();
		SOPVector stack = QFT(mode_b, n_bit, adjungate, approx);

		SOPVector rot = ControlledRotation(n_bit, mode_a, mode_b, adjungate, approx);
		stack.insert(stack.end(), rot.begin(), rot.end());

		SOPVector iqft = iQFT(mode_b, n_bit, adjungate, approx);
		stack.insert(stack.end(), iqft.begin(), iqft.end());
		cout << "Initialized QFT adder circuit with a stacksize of " << stack.size() << " operators.\n";
		return stack;
	}

	SOPVector const_add(const long_integer& bin_a, size_t n_bit, size_t mode_b, bool adjungate, size_t approx) {
		using namespace GateOperators;

		// Build A_j
		MPO M;
		for (int j = 0; j < n_bit; j++) {
			long double x = compute_x(bin_a, j);
			if (x > pow(2, -int(approx))) {
				auto R = make_shared<Rx>(x, adjungate);
				M.push_back(R, mode_b + n_bit - 1 - j);
			}
		}
		SOP S(M, 1.0);
		SOPVector stack;
		stack.append(S);
		return stack;
	}

	SOPVector const_subst(const long_integer& bin_a, size_t n_bit, size_t mode_b,
		bool adjungate, size_t approx) {
		return const_add(bin_a, n_bit, mode_b, !adjungate, approx);
	}

	SOPVector c_const_add(const long_integer& bin_a, size_t n_bit, size_t mode_b, bool adjungate, size_t approx,
		size_t control) {
		SOPVector stack = const_add(bin_a, n_bit, mode_b, adjungate, approx);
		return GateOperators::DistribSingleControl(stack, control);
	}

	SOPVector cc_const_add(const long_integer& a, size_t n_bit, size_t mode_b,
		bool adjungate, size_t approx,
		size_t control1, size_t control2) {
		SOPVector stack = const_add(a, n_bit, mode_b, adjungate, approx);
		return GateOperators::DistribDoubleControl(stack, control1, control2);
	}

	// Phi MAC from Fig 8 of Pavlidias & Gizopoulus
	// If control on, multiply qubit int x by const int a and accum with qubit int b
	// Qubit b must have 2n bits, qubit x has n bits
	SOPVector c_phi_mac(const long_integer& a, size_t n_bit, size_t mode_x, size_t mode_b,
		bool adjungate, size_t approx, size_t control) {

		SOPVector stack;
		for (size_t xi = 0; xi < n_bit; ++xi) {
			size_t controlx = mode_x + xi;
			long_integer tmp(a);
			tmp.bitshift(xi);
			SOPVector ccadd = cc_const_add(tmp, 2 * n_bit, mode_b, adjungate, approx, control, controlx);
			stack.append(ccadd);
		}

		return stack;
	}

	// Double controlled phi add a mod N from Figure 5 from Beauregard 2003
	SOPVector cc_add_mod(const long_integer& a, const long_integer& N, size_t n_bit,
		size_t mode_b, size_t mode_0, bool adjungate, size_t approx,
		size_t control1, size_t control2) {
		using namespace GateQFT;
		using namespace GateOperators;

		size_t largest_bit = mode_b;

		SOPVector stack;
		stack.append(cc_const_add(a, n_bit, mode_b, adjungate, approx, control1, control2));

		stack.append(const_subst(N, n_bit, mode_b, adjungate, approx));

		stack.append(iQFT(mode_b, n_bit, adjungate, approx));
		stack.append(CNot(largest_bit, mode_0));
		stack.append(QFT(mode_b, n_bit, adjungate, approx));

		stack.append(c_const_add(N, n_bit, mode_b, adjungate, approx, mode_0));
		stack.append(cc_const_add(a, n_bit, mode_b, !adjungate, approx, control1, control2));
		stack.append(iQFT(mode_b, n_bit, adjungate, approx));
		stack.append(MPO(X, largest_bit));
		stack.append(CNot(largest_bit, mode_0));
		stack.append(MPO(X, largest_bit));
		stack.append(QFT(mode_b, n_bit, adjungate, approx));
		stack.append(cc_const_add(a, n_bit, mode_b, adjungate, approx, control1, control2));

		return stack;
	}

	// Controlled CMULT(a)MOD(N) from Figure 6 from Beauregard 2003
	SOPVector c_mult_mod(const long_integer& a, const long_integer& N,
		size_t n_bit, size_t mode_x, size_t mode_b, size_t mode_0,
		bool adjungate, size_t approx, size_t c) {
		using namespace GateQFT;

		SOPVector stack;
		stack.append(QFT(mode_b, n_bit, adjungate, approx));

		for (int i = 0; i < n_bit - 1; ++i) {
			long_integer ashift(a);
			ashift.bitshift_mod(i, N);
			size_t c2 = mode_x + n_bit - 1 - i;

			cout << "For i = " << i << " add the number: \n";
			ashift.print_4byte();
			cout << "Control mode for x:" << c2 << endl;

			SOPVector add_mod = cc_add_mod(ashift, N, n_bit, mode_b, mode_0,
				adjungate, approx, c, c2);
			stack.append(add_mod);
		}

		stack.append(iQFT(mode_b, n_bit, adjungate, approx));

		return stack;
	}

	SOPVector Ua(const long_integer& a, const long_integer& a_inv, const long_integer& N,
		size_t n_bit, size_t mode_x, size_t mode_b, size_t mode_0,
		bool adjungate, size_t approx, size_t control) {

		SOPVector stack;
		stack.append(c_mult_mod(a, N, n_bit, mode_x, mode_b, mode_0, adjungate, approx, control));

		/// Do a n-qubit c_swap from one-qubit c_swap
		for (size_t i = 0; i < n_bit; ++i) {
			size_t idx_b = mode_b + i;
			size_t idx_x = mode_x + i;
			stack.append(c_swap(control, idx_x, idx_b));
		}

		/// Run the adjoint of c_mult_mod (therefore reverse and adjungate)
		auto tmp = c_mult_mod(a_inv, N, n_bit, mode_x, mode_b, mode_0, !adjungate, approx, control);
		reverse(tmp.begin(), tmp.end());
		stack.append(tmp);

		return stack;
	}

	//@TODO: Add extended euclid's

	SOPVector c_swap(size_t control, size_t target1, size_t target2) {
		using namespace GateOperators;
		SOPVector stack;
		stack.append(CNot(target2, target1));
		stack.append(Toffoli(control, target1, target2));
		stack.append(CNot(target2, target1));
		return stack;
	}

	SOPVector exp_mod(const long_integer& a, const long_integer& N,
		const mctdhBasis& basis, bool adjungate, size_t approx) {
		/// f has to be 4n+2 qubits
		size_t f = basis.nLeaves();
		size_t n_bit = (f - 2) / 4;
		assert(f == 4 * n_bit + 2);
		size_t mode_measure_space = 0;
		size_t mode_x = 2 * n_bit;
		size_t mode_b = mode_x + n_bit;
		size_t n_bit_b = n_bit + 1;
		size_t mode_0 = mode_b + n_bit_b;

		SOPVector stack;
		assert(0);
		for (size_t i = 0; i < 2 * n_bit; ++i) {
			long_integer ashift = a;
			size_t exponent = pow(2, i);
			/// @TODO: get a modular exponential operation
//			ashift.exp_mod(exponent, N);

			/// @TODO: get a_invs
			long_integer a_shift_inv;

			size_t c = mode_measure_space + 2 * n_bit - 1 - i;
			stack.append(Ua(ashift, a_shift_inv, N, n_bit,
				mode_x, mode_b, mode_0, adjungate, approx, c));
		}
		return stack;
	}
} // end namespace GateArithmetic
