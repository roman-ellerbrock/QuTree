//
// Created by Roman Ellerbrock on 4/27/20.
//
#include "Arithmetic.h"
#include "QuantumCircuits.h"
#include "../InTimeQuantumInstruction.h"
#include "../ConditionalQuantumInstruction.h"
#include "../OutputInstruction.h"
#include "Shor.h"
#include "QFT.h"
#include "../RandomProjectorInstruction.h"

#include <boost/multiprecision/cpp_int.hpp>
#include <boost/integer/mod_inverse.hpp>
using namespace boost::multiprecision;


namespace Circuits {

	YAML::Node Rspecial() {
		string str = "instr: intimequantum\nclassical_reg: measurement\ntype: Rspecial\nreg: [M]";
		stringstream ss(str);
		YAML::Node node = YAML::Load(ss);
		return node;
	}

    QuantumCircuit ShorFull(const Register& reg,
                        const cpp_int& a, const cpp_int& N) {

        /// prepare x, b and zero
        const Register& M = reg["M"];
        const Register& x = reg["x"];
        const Register& b = reg["b"];
        const Register& zero = reg["zero"];
        //size_t approx = 48;
        size_t approx = 0;

        cout << "Running Shors Algorithm with 4n+2 qubits.\n";
        cout << "N = " << N << endl;
        cout << "a = " << a << endl;
        cout << "size(M) = " << M.size() << endl;
        cout << "size(x) = " << x.size() << endl;
        cout << "size(b) = " << b.size() << endl;
        cout << "size(zero) = " << zero.size() << endl;

        // check phases
        //Shor::predict_phases(a, N, x.size());

        QuantumCircuit circ;

        /// apply Hadamard
        SOPVectorcd Hc(HadamardChain(M));
        circ.emplace_back(make_shared<QuantumInstruction>(Hc, Hc));

        /// Prepare x=1
        SOPVectorcd Xc(Set_Number(x, long_integer(1, x.size())));
        circ.emplace_back(make_shared<QuantumInstruction>(Xc, Xc));

        /// Call Ua gate cycles
        for (int l = M.size() - 1; l >= 0; --l) {
            int control = M.front() + l;
            Register c(control, 1);

            auto uca = Circuits::cUa<cpp_int>(a, c, x, b, zero, N, l, false, approx);
            auto uca_adj = Circuits::cUa<cpp_int>(a, c, x, b, zero, N, l, true, approx);
            circ.emplace_back(make_shared<QuantumInstruction>(uca, uca_adj));
        }
        // The version below doesn't work becasue this Shors fn already adds the adj
        //SOPVectorcd Uas(Shors(M, x, b, zero, a, N, false, approx));
        //SOPVectorcd Uas_adj(Shors(M, x, b, zero, a, N, true, approx));
        //circ.emplace_back(make_shared<QuantumInstruction>(Uas, Uas_adj));

        /// iQFT
        auto iqft = Circuits::iQFT(M, false);
        auto iqft_adj = Circuits::iQFT(M, true);
        circ.emplace_back(make_shared<QuantumInstruction>(iqft, iqft_adj));

        return circ;
    }

	QuantumCircuit Shor(const Register& reg,
		const cpp_int& a, const cpp_int& N) {

		/// prepare x, b and zero
		const Register& M = reg["M"];
		const Register& x = reg["x"];
		const Register& b = reg["b"];
		const Register& zero = reg["zero"];
		size_t approx = 48;

/*
		Shor::predict_phases(a, N, x.Size());
		Sieve::quadraticSieve(N);
		cout << "done.\n";
		getchar();
		*/

		/// apply Hadamard
		SOPVectorcd Hc(HadamardChain(M));

		QuantumCircuit circ;
		/// Prepare x=1
		SOPVectorcd Xc(Set_Number(x, long_integer(1, x.size())));
		circ.emplace_back(make_shared<QuantumInstruction>(Xc, Xc));
		int L = 2 * x.size();
//		L=2;
		for (int l = L - 1; l >= 0; --l) {
			circ.emplace_back(make_shared<QuantumInstruction>(Hc, Hc));

			cout << "cUa:\n";
			/// Ua^{2^l}
			auto t = a;
			for (int i = 0; i < l; ++i) {
				t = (t * t) % N;
			}
			if (t != 1) {
				auto uca = Circuits::cUa<cpp_int>(a, M, x, b, zero, N, l, false, approx);
				auto uca_adj = Circuits::cUa<cpp_int>(a, M, x, b, zero, N, l, true, approx);
				circ.emplace_back(make_shared<QuantumInstruction>(uca, uca_adj));
			}

			size_t j = L - l;
			if (j > 1) { /// (j >= 2) ? R'
				auto node = Rspecial();
				circ.emplace_back(make_shared<InTimeQuantumInstruction>(node, reg));
			}

			/// Hadamard
			circ.emplace_back(make_shared<QuantumInstruction>(Hc, Hc));

			/// Measure (append to measurement list)
			circ.emplace_back(make_shared<MeasurementInstruction>(M, "measurement", true));
			circ.emplace_back(make_shared<MeasurementInstruction>(M, "m0", false));

			/// Restore measurement qubit |0> state for next cycle
			vector<Register> tars = {M};
			circ.emplace_back(make_shared<ConditionalQuantumInstruction>("m0", "1", "X", tars));

			circ.emplace_back(make_shared<OutputInstruction>());
		}

//		auto qft = Circuits::QFT(x, false);
//		auto qft_adj = Circuits::QFT(x, true);
//		circ.emplace_back(make_shared<QuantumInstruction>(qft, qft_adj));
//		circ.emplace_back(make_shared<OutputInstruction>());

		return circ;
	}

	SOPVectorcd cUa(const YAML::Node& node,
		const Register& control, const Register& x, const Register& b,
		const Register& zero, bool adjoint, size_t approx) {
		/**
		 *
		 */

		auto a = (cpp_int) eval<size_t>(node, "a");
		auto N = (cpp_int) eval<size_t>(node, "N");
//		auto as = eval<size_t>(node, "a");
//		auto Ns = eval<size_t>(node, "N");
		auto l = eval<size_t>(node, "l");
//		long_integer a(as, b.Size());
//		long_integer N(Ns, b.Size());
		return Circuits::cUa<cpp_int>(a, control, x, b, zero, N, l, adjoint, approx);
	}

	SOPVectorcd cMULTmod(const YAML::Node& node,
		const Register& control, const Register& x, const Register& b,
		const Register& zero, bool adjoint) {

/*		auto as = eval<size_t>(node, "a");
		auto Ns = eval<size_t>(node, "N");
		long_integer a(as, b.Size());
		long_integer N(Ns, b.Size());
*/
		auto a = (cpp_int) eval<size_t>(node, "a");
		auto N = (cpp_int) eval<size_t>(node, "N");

		size_t approx = 48;
		return c_mult_mod<cpp_int>(b, x, a, N, zero.front(), control.front(), adjoint, approx);
	}

	SOPVectorcd ccADDmod(const YAML::Node& node,
		const Register& c1, const Register& c2, const Register& b,
		const Register& zero, bool adjoint) {

/*		auto as = eval<size_t>(node, "a");
		auto Ns = eval<size_t>(node, "N");
		long_integer a(as, b.Size());
		long_integer N(Ns, b.Size());*/
		auto a = (cpp_int) eval<size_t>(node, "a");
		auto N = (cpp_int) eval<size_t>(node, "N");

		size_t approx = 64;
		return cc_add_mod<cpp_int>(b, zero.front(), a, N, c1.front(),
			c2.front(), adjoint, approx);
	}

	SOPVectorcd Rspecial(const YAML::Node& node,
		const Register& tar, const QVMState& state,
		bool adjoint) {

		if (node["classical_reg"]) {
			auto par = eval<string>(node, "classical_reg");
			const Measurements::Measurement& meas = state.classical_reg_.at(par);
			assert(tar.size() == 1);
			size_t target = tar.front();
//			cout << "target: " << target << endl;
			MLOcd M(Rspecial(meas, adjoint), target);
			SOPcd H(M, 1.);
//			cout << "Hsize: " << H.size() << endl;
			return SOPVectorcd(H);
		} else {
			cerr << "Missing classical_reg flag" << endl;
			exit(1);
		}
	}

	shared_ptr<LeafOperatorcd> Rspecial(
		const Measurements::Measurement& m, bool adjoint) {

		double phi = 0;
		size_t j = m.size();
//		cout << "Rspecial:\n";
//		cout << "j = " << j << endl;
//		for (size_t i = 0; i < j; ++i) {
//			cout << i << " - " << m[i] << endl;
//		}
		for (int k = 2; k <= j; ++k) {
			auto add = m[j - k] / pow(2., k);
			phi += add;
//			cout << j - k << " " << add << " " << phi << endl;
		}
		return make_shared<Rx>(-phi, adjoint);
	}

	SOPVectorcd transversalIsing(size_t steps) {
		SOPVectorcd I;

		MLOcd Hs;
		for (size_t q = 0; q < 10; ++q) {
			Hs.push_back(Hadamard, q);
		}
		I.append(Hs);
		I.append(ising(0, 1, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(0, 1, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(2, 3, steps));
		I.append(ising(1, 2, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(4, 5, steps));
		I.append(ising(3, 4, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));
		I.append(ising(6, 7, steps));
		I.append(ising(5, 6, steps));
		I.append(ising(8, 9, steps));
		I.append(ising(7, 8, steps));

		return I;

/*		SOPVectorcd J;
		mt19937 gen(0101);
		uniform_real_distribution<double> dist(0., 1.);
		double p = 0.1;
		cout << "p=?\n";
		cin >> p;
		for (size_t i = 0; i < I.size(); ++i) {
			J.push_back(I[i]);
			MLOcd M;
			for (size_t q = 0; q < 10; ++q) {
				if (!I[i].isActive(q)) { continue; }
				cout << "i = " << i << " | " << q << endl;
				if (dist(gen) < p) {
					M.push_back(randomProject(gen), q);
				}
			}
			if (M.size() > 0) {
				J.append(M);
			}
		}
		cout << "size(J) = " << J.size() << " | size(I) = " << I.size() << endl;
		return J;
		*/
	}

	SOPVectorcd isingStep(const Register& reg, size_t steps) {
		SOPVectorcd circ;
		for (size_t i = reg.front(); i <= reg.back(); i += 2) {
			circ.append(ising(i, i+1, steps));
		}
		for (size_t i = reg.front()+1; i < reg.back(); i += 2) {
			circ.append(ising(i, i+1, steps));
		}
		return circ;
	}

	QuantumCircuit isingIntegration(const Register& reg, size_t time, double p, size_t steps) {
		QuantumCircuit circ;
//		circ.push_back(make_shared<OutputInstruction>(*os));
		cout << "time, steps: " << time * steps << endl;
		for (size_t d = 0; d < time * steps; ++d) {
			SOPVectorcd sops(isingStep(reg, steps));
			circ.push_back(make_shared<QuantumInstruction>(sops, sops));
//			circ.push_back(make_shared<OutputInstruction>(*os));
		}
		return circ;
	}

	QuantumCircuit isingModel(double p, size_t steps) {
		QuantumCircuit circ;
		auto I = transversalIsing(steps);

		mt19937 gen(time(nullptr));
		uniform_real_distribution<double>dist(0., 1.);
		for (const auto& i : I) {
			SOPVectorcd isop(i);
			circ.push_back(make_shared<QuantumInstruction>(isop, isop));

			for (size_t q = 0; q < 10; ++q) {
				if (!i.isActive(q)) { continue; }
				if (dist(gen) < p) {
					circ.push_back(make_shared<RandomProjectorInstruction>(q));
				}
			}
		}
		return circ;
	}

}

