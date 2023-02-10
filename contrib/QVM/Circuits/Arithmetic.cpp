#include "Arithmetic.h"
#include "GateOperators.h"
#include "QFT.h"

namespace Circuits {

	SOPVectorcd add_qft(const Register& a, const Register& b, bool adjungate, size_t approx) {
		SOPVectorcd stack = QFT(b, adjungate, approx);
		stack.append(controlledRotation(a, b, adjungate, approx));
		stack.append(iQFT(b, adjungate, approx));
		return stack;
	}

	template <class Z>
	SOPVectorcd const_add(const Register& b, const Z& a, bool adjungate, size_t approx) {
		MLOcd M;
		size_t n_bit = b.size();
		for (int j = 0; j < n_bit; j++) {
			long double x = compute_x(long_integer((size_t) a, 64), j);
			auto R = make_shared<Rx>(x, adjungate);
			M.push_back(R, b.back() - j);
		}
		SOPcd S(M, 1.0);
		SOPVectorcd stack;
		stack.append(S);
		return stack;
	}

	template <class Z>
	SOPVectorcd c_const_add(const Register& b, const Z& a, size_t control,
		bool adjungate, size_t approx) {
		SOPVectorcd stack = const_add<Z>(b, a, adjungate, approx);
		return distribSingleControl(stack, control);
	}

	template <class Z>
	SOPVectorcd cc_const_add(const Register& b, const Z& a,
		size_t control1, size_t control2, bool adjungate, size_t approx) {
		SOPVectorcd stack = const_add<Z>(b, a, adjungate, approx);
		return distribDoubleControl(stack, control1, control2);
	}

	template <class Z>
	SOPVectorcd const_subst(const Register& b, const Z& a, bool adjungate, size_t approx) {
		return const_add<Z>(b, a, !adjungate, approx);
	}

	template <class Z>
	SOPVectorcd cc_add_mod(const Register& b, size_t zero,
		const Z& a, const Z& N,
		size_t control1, size_t control2, bool adjungate, size_t approx) {
		/// cc-Modular adder (b) -> ((b + a) % N) where a, N are classical numbers
		/// This design follows Ref. [1] Fig. 5. See paper for an explanation.
		/// Circuit expects a quint b in Fourier representation

		size_t largest_bit = b.largestBit();
		SOPVectorcd stack;
		stack.append(cc_const_add<Z>(b, a, control1, control2, adjungate, approx));

		stack.append(const_subst<Z>(b, N, adjungate, approx));

		stack.append(iQFT(b, adjungate, approx));
		stack.append(CNot(largest_bit, zero));
		stack.append(QFT(b, adjungate, approx));

		stack.append(c_const_add<Z>(b, N, zero, adjungate, approx));
		stack.append(cc_const_add<Z>(b, a, control1, control2, !adjungate, approx));
		stack.append(iQFT(b, adjungate, approx));
		stack.append(MLOcd(X, largest_bit));
		stack.append(CNot(largest_bit, zero));
		stack.append(MLOcd(X, largest_bit));
		stack.append(QFT(b, adjungate, approx));
		stack.append(cc_const_add<Z>(b, a, control1, control2, adjungate, approx));

		return stack;
	}

	template <class Z>
	SOPVectorcd c_mult_mod(const Register& b, const Register& x,
		Z a, const Z& N,
		size_t mode_0, size_t c, bool adjungate, size_t approx) {
		/// cc-Modular adder (b) -> ((b + a*x) % N) where a, N are classical numbers
		/// This design follows Ref. [1] Fig. 6. See paper for an explanation.
		/// Circuit expects a quint b computational (not Fourier) representation

//		assert(x.Size() + 1 == b.Size());
		SOPVectorcd stack;
		/// Transform to Fourier space
		stack.append(QFT(b, adjungate, approx));

		/// Perform sequence of modular adders for qubit in x
		for (int i = 0; i < x.size(); ++i) {
//			long_integer ashift(a);
//			ashift.bitshift_mod(i, N);
			Z ashift = a;
			a = (2 * a) % N;
			size_t c2 = x.smallestBit() - i;
//			cout << "a=" << ashift << endl;
			stack.append(cc_add_mod<Z>(b, mode_0, ashift, N, c, c2, adjungate, approx));
		}

		/// Transform back to original space
		stack.append(iQFT(b, adjungate, approx));
		return stack;
	}

	template <class Z>
	SOPVectorcd Ua(const Register& x, const Register& b,
		size_t mode_0, size_t control, const Z& a, const Z& N,
		bool adjungate, size_t approx) {
		/// cc-Modular adder (x, 0) -> ((a*x) % N, 0) where a, N are classical numbers.
		/// This design follows Ref. [1] Fig. 7. See paper for an explanation.
		/// Circuit expects a quint x computational (not Fourier) representation.
		/// This is the main building block for Shor's algorithm.

		SOPVectorcd stack;
		/// Do a regular c_mult_mod (b -> b + x*a % N)
		stack.append(c_mult_mod<Z>(b, x, a, N, mode_0, control, adjungate, approx));

		/// Do a n-qubit c_swap from one-qubit c_swap
		for (size_t i = 0; i < x.size(); ++i) {
			size_t idx_b = b.back() - i;
			size_t idx_x = x.back() - i;
			stack.append(c_swap(control, idx_x, idx_b));
		}

		/// Run the adjoint of c_mult_mod (therefore reverse and adjungate)
		auto a_inv = boost::integer::mod_inverse(a, N);
//		auto a_inv = LargeEuclidicInverse(a, N);
//		auto a_inv = EuclidicInverse(a, N);
//		cout <<  "a_inv\n";
//		a_inv.print();
		auto tmp = c_mult_mod<Z>(b, x, a_inv, N, mode_0, control, !adjungate, approx);
		reverse(tmp.begin(), tmp.end());
		stack.append(tmp);

		return stack;
	}

	template <class Z>
	SOPVectorcd Shors(const Register& M, const Register& x, const Register& b,
		const Register& zero, Z a, const Z& N,
		bool adjungate, size_t approx) {
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

		SOPVectorcd stack;
		cout << "Shors Algorithm with 4n+2 qubits.\n";
		cout << "N = " << N << endl;
		cout << "a = " << a << endl;
		cout << "size(M) = " << M.size() << endl;
		cout << "size(x) = " << x.size() << endl;
		cout << "size(b) = " << b.size() << endl;
		cout << "size(zero) = " << zero.size() << endl;

		for (size_t i = 0; i < M.size(); ++i) {
			int control = M.back() - i;
			Register c(control, 1);
//			int control = M.Last();
			cout << control <<  endl;
//			a.print();
//			a.print_4byte();
			stack.append(Ua<Z>(x, b, zero.smallestBit(), control, a, N, adjungate, approx));
//			stack.append(cUa<Z>(a, c, x, b, zero, N, i, adjungate, approx));
//			a = mult_mod(a, a, N);
			a = (a * a) % N;
//			if (i == 3) { break; }
		}

		return stack;
	}

	SOPVectorcd Rprime(const Register& t, const vector<size_t>& measurements, size_t j, bool adjoint) {
		/**
		 * \brief combinations of rotations Rj for semiclassical implementation of Shor's.
		 *
		 *
		 */
		/// calculate angle
		double x = 0.;
		for (size_t k = 2; k < j; ++k) {
			x += measurements[j - k];
		}
		auto R = make_shared<Rx>(-x, adjoint);
		MLOcd Rvec;
		Rvec.push_back(R, t.front());
		SOPVectorcd stack;
		stack.append(Rvec);
		return stack;
	}

	template <class Z>
	SOPVectorcd cUa(Z a, const Register& control, const Register& x,
		const Register& b, const Register& zero, const Z& N,
		size_t l, bool adjungate, size_t approx) {
		SOPVectorcd stack;
		/**
		 * \brief Single controlled Ua_a^{2^l}
		 *
		 * It is a single controlled Ua operation for a given power of l
		 * (see Fig. 2 in https://arxiv.org/pdf/quant-ph/0001066.pdf)
		 *
		 * @param a integer modified for this layer
		 *
		 */

//		cpp_int z = a.convert_cpp();
//		cpp_int NN = N.convert_cpp();
		for (size_t i = 0; i < l; ++i) {
//			a = mult_mod(a, a, N);
			cout << i << " -- " << a << endl;
			a = (a * a) % N;
		}
		cout << "a = " << a << endl;
//		long_integer tmp((size_t) z, a.size());
//		a = tmp;
		stack.append(Ua<Z>(x, b, zero.front(), control.front(), a, N, adjungate, approx));
		return stack;
	}

	SOPVectorcd Set_Number(const Register& x, const long_integer& a) {
		assert(x.size() == a.size());
		MLOcd M;
		for (size_t k = 0; k < x.size(); ++k) {
			if (a[k] == 1) {
				M.push_back(Circuits::X, x.back() - k);
			}
		}
		SOPcd S(M);
		SOPVectorcd stack;
		stack.append(S);
		return stack;
	}

	/// ==================================================================================

	SOPVectorcd Random_Number(size_t mode_a, size_t n_bit) {
//		srand(time(nullptr));
//		int seed = rand();
		int seed = 4661242;
		uniform_real_distribution<double> dist(0., 1.);
		mt19937 rand_gen(seed);
		MLOcd M;
		for (size_t k = 0; k < n_bit; ++k) {
			if (dist(rand_gen) > 0.5) {
				M.push_back(Circuits::X, k + mode_a);
			}
		}
		SOPcd S(M);
		SOPVectorcd stack;
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

	SOPVectorcd add_qft(size_t n_bit, size_t mode_a, size_t mode_b, bool adjungate, size_t approx) {

		cout << "Initializing QFT adder circuit ...\n";
		cout.flush();
		SOPVectorcd stack = QFT(mode_b, n_bit, adjungate, approx);

		SOPVectorcd rot = controlledRotation(n_bit, mode_a, mode_b, adjungate, approx);
		stack.insert(stack.end(), rot.begin(), rot.end());

		SOPVectorcd iqft = iQFT(mode_b, n_bit, adjungate, approx);
		stack.insert(stack.end(), iqft.begin(), iqft.end());
		cout << "Initialized QFT adder circuit with a stacksize of " << stack.size() << " operators.\n";
		return stack;
	}

	SOPVectorcd c_swap(size_t control, size_t target1, size_t target2) {
		SOPVectorcd stack;
		stack.append(CNot(target2, target1));
		stack.append(Toffoli(control, target1, target2));
		stack.append(CNot(target2, target1));
		return stack;
	}

	template SOPVectorcd const_add(const Register& b, const cpp_int& a,
		bool adjungate, size_t approx = 0);

	template SOPVectorcd c_const_add(const Register& b, const cpp_int& a,
		size_t control, bool adjungate, size_t approx);

	template SOPVectorcd cc_const_add(const Register& b, const cpp_int& a,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	template SOPVectorcd const_subst(const Register& b, const cpp_int& a,
		bool adjungate, size_t approx);

	template SOPVectorcd cc_add_mod(const Register& b, size_t zero,
		const cpp_int& a, const cpp_int& N,
		size_t control1, size_t control2,
		bool adjungate, size_t approx);

	template SOPVectorcd c_mult_mod(const Register& b, const Register& x,
		cpp_int a, const cpp_int& N,
		size_t mode_0, size_t control, bool adjungate, size_t approx);

	template SOPVectorcd Ua(const Register& x, const Register& b, size_t mode_0,
		size_t control, const cpp_int& a, const cpp_int& N,
		bool adjungate, size_t approx);

	template SOPVectorcd cUa(cpp_int a, const Register& control, const Register& x,
		const Register& b, const Register& zero, const cpp_int& N,
		size_t l, bool adjungate, size_t approx);

	template SOPVectorcd Shors(const Register& M, const Register& x, const Register& b,
		const Register& zero, cpp_int a, const cpp_int& N, bool adjungate, size_t approx);

} // end namespace GateArithmetic
