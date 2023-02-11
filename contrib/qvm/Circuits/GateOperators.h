//
// Created by Roman Ellerbrock on 2019-04-19.
//

#ifndef MCTDH_QUBIT_OPERATORS_H
#define MCTDH_QUBIT_OPERATORS_H

#include "Util/QMConstants.h"
#include "Register.h"
#include "TreeShape/Tree.h"
#include "TreeOperators/SOPVector.h"
#include "TreeOperators/LeafFunction.h"

namespace Circuits {

	/*!
	 * \namespace GateOperators
	 * \brief This is a library of Gate operators to build quantum circuits.
	 *
	 * The nomenclature of the operations is related to the standard notation in the
	 * community. A detailed list of many relevant gates can be found in [1, 2].
	 *
	 * [1] "Massively Parallel ...Eleven Years Later", Computer Physics Communication (2019)
	 * [2] "..." https://arxiv.org/pdf/quant-ph/0511250.pdf
	 *
	 */

	/// =========================
	/// Single-Qubit operators
	/// =========================

	/// These are fundamental Gate Operations
	void Hadamard(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void X(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Y(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Z(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void S(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Sdagger(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void set0(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void set1(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void identity(const LeafInterface& grid, Tensorcd& xPhi, const Tensorcd& phi);
	Matrixcd T();
	Matrixcd sqrtX();
	Matrixcd sqrtY();
	Matrixcd sigma_x();
	Matrixcd sigma_y();
	Matrixcd sigma_z();
	Matrixcd rX(double theta);
	Matrixcd rY(double theta);
	Matrixcd rZ(double theta);
	Matrixcd U(double theta, double phi, double lambda);

	SOPVectorcd ising(size_t q0, size_t q1, size_t steps);
    SOPcd givens(double theta, size_t c, size_t t);
    SOPcd cGivens(double theta, size_t control, size_t q1, size_t q2);
	LeafFunctioncd randomProject(mt19937& gen);

	class PrimitiveProjector: public LeafOperatorcd {
	public:
		explicit PrimitiveProjector(size_t idx)	: idx_(idx) {}
		~PrimitiveProjector() = default;
		void apply(const LeafInterface& grid, Tensorcd& PPhi, const Tensorcd& Phi) const override {
			const TensorShape& tdim = Phi.shape();
			PPhi.zero();
			assert(idx_ < tdim.lastBefore());
			for (size_t n = 0; n < tdim.lastDimension(); ++n) {
				PPhi(idx_, n) = Phi(idx_, n);
			}
		}
	private:
		size_t idx_;
	};

	class Rk: public LeafOperatorcd {
		/// This class represents the R(k) Phase shift operation
	public:
		/// The phase will be converted to a factor: exp(2*pi*im*phase_)
		explicit Rk(size_t n, bool adjungate = false);
		~Rk() = default;
		/// Efficient applying routine
		void apply(const LeafInterface& grid, Tensorcd& RAcoeff, const Tensorcd& Acoeff) const override;
	private:
		complex<double> factor;
	};

	class Rx: public LeafOperatorcd {
		/// This class represents the R(k) Phase shift operation
	public:
		/// The phase will be converted to a factor: exp(2*pi*im*phase_)
		explicit Rx(long double n, bool adjungate = false);
		~Rx() = default;
		/// Efficient applying routine
		void apply(const LeafInterface& grid, Tensorcd& RAcoeff, const Tensorcd& Acoeff) const override;
	private:
		complex<double> factor;
	};

	class SetBlack: public LeafOperatorcd {
	public:
		SetBlack(double epsilon_, bool adjungate_)
			: epsilon(epsilon_), adjungate(adjungate_) {}
		~SetBlack() = default;
		void apply(const LeafInterface& grid, Tensorcd& xPhi,
			const Tensorcd& phi) const override {
			if (adjungate) {
				for (size_t n = 0; n < 2; ++n) {
					/// white
					xPhi(0, n) = epsilon * phi(0, n) + sqrt(1. - pow(epsilon, 2)) * phi(1, n);
					/// black
					xPhi(1, n) = 0.;
				}
			} else {
				for (size_t n = 0; n < 2; ++n) {
					/// white
					xPhi(0, n) = epsilon * phi(0, n);
					/// black
					xPhi(1, n) = sqrt(1. - pow(epsilon, 2)) * phi(0, n);
				}
			}
		}

	private:
		bool adjungate;
		double epsilon;
	};

	/// =========================
	/// Two-Qubit operators
	/// =========================

	/// Build single-control gates (e.g. cX == CNOT gate)
	SOPcd makeCGate(const SOPcd& sop, size_t control);

	/// Overloads purely for typecasting convenience
	inline SOPcd makeCGate(const MLOcd& mpo, size_t control) { return makeCGate(SOPcd(mpo), control); }

	inline SOPcd makeCGate(const LeafFunctioncd& spo, size_t control, size_t target) {
		return makeCGate(MLOcd(spo, target), control);
	}

	inline SOPcd makeCGate(const LeafFuncd& spo, size_t control, size_t target) {
		return makeCGate(MLOcd(spo, target), control);
	}

	inline SOPcd makeCGate(shared_ptr<LeafOperatorcd> spo, size_t control, size_t target) {
		MLOcd M;
		M.push_back(spo, target);
		return makeCGate(M, control);
	}

	/// Build double-control gates (e.g. ccX == CCNOT == Toffoli gate)
	SOPcd makeCCGate(const SOPcd& sop, size_t control1, size_t control2);

	/// Overloads purely for typecasting convenience
	inline SOPcd makeCCGate(const MLOcd& mpo, size_t control1, size_t control2) {
		return makeCCGate(SOPcd(mpo),
			control1,
			control2);
	}

	inline SOPcd makeCCGate(const LeafFunctioncd& spo, size_t control1, size_t control2, size_t target) {
		return makeCCGate(MLOcd(spo, target), control1, control2);
	}

	inline SOPcd makeCCGate(const LeafFuncd& spo, size_t control1, size_t control2, size_t target) {
		return makeCCGate(MLOcd(spo, target), control1, control2);
	}

	inline SOPcd makeCCGate(shared_ptr<LeafOperatorcd> spo, size_t control1, size_t control2, size_t target) {
		MLOcd M;
		M.push_back(spo, target);
		return makeCCGate(M, control1, control2);
	}

	/// Convenience functions
	SOPcd CNot(const Register& r1, const Register& r2);
	SOPcd CNot(size_t control, size_t target);
	SOPcd CZ(const Register& r1, const Register& r2);
	SOPcd CZ(size_t control, size_t target);
	SOPcd Uk(const Register& c, const Register& t, size_t k, bool conj_phase = false);
	SOPcd Uk(size_t control, size_t target, size_t k, bool conj_phase = false);
	SOPcd Ux(const Register& c, const Register& t, double x, bool conj_phase = false);
	SOPcd Ux(size_t c, size_t t, double x, bool conj_phase = false);
	SOPcd Toffoli(const Register& c1, const Register& c2, const Register& t);
	SOPcd Toffoli(size_t control1, size_t control2, size_t target);

	SOPVectorcd SWAP(const size_t c, const size_t t);
	SOPVectorcd iSWAP(const size_t c, const size_t t, double theta, bool adj);
	SOPcd iSWAPSOP(const size_t c, const size_t t, double theta, bool adj);

	SOPVectorcd distribute(const LeafFunPaircd& h, const Register& reg, bool adjoint);
	SOPVectorcd distribute(const LeafFuncd& h, const Register& reg, bool adjoint);
	SOPVectorcd distibute(const Matrixcd& h, const Register& reg, bool adjoint);
	SOPVectorcd distribSingleControl(const SOPVectorcd& stack, size_t control);
	SOPVectorcd distribDoubleControl(const SOPVectorcd& stack, size_t control1, size_t control2);
	SOPVectorcd distribute(const Register& reg, function<void(const LeafInterface&, Tensorcd&, const Tensorcd&)> f);

	/// =========================
	/// High-level operations
	/// =========================

	/// Building blocks for logic circuits
	SOPVectorcd Carry(const Register& c1, const Register& a, const Register& b, const Register& c2);
	SOPVectorcd Carry(size_t carry1, size_t a, size_t b, size_t carry2);
	SOPVectorcd RCarry(const Register& c1, const Register& a, const Register& b, const Register& c2);
	SOPVectorcd RCarry(size_t carry1, size_t a, size_t b, size_t carry2);
	SOPVectorcd SUM(const Register& a, const Register& b, const Register& sum);
	SOPVectorcd SUM(size_t a, size_t b, size_t sum);

	SOPVectorcd Adder_circuit(const Tree& tree, size_t n_bit, size_t mode_a, size_t mode_b,
		size_t mode_carry);

	SOPVectorcd RandomSuperposition(const Tree& tree, size_t start, size_t end, size_t n);
	SOPVectorcd HadamardChain(const Register& reg);

	// This is the CNot-Chain-benchmark
	SOPVectorcd CNotChain(const Register& reg);
	SOPVectorcd CNotBenchmark(const Register& reg);
	SOPVectorcd CNotGridRow(const Tree& tree, size_t n_row, size_t n_col, size_t offset = 0, size_t shift = 0);
	SOPVectorcd CNotGridCol(const Tree& tree, size_t n_row, size_t n_col, size_t offset = 0, size_t shift = 0);
	SOPVectorcd CNotGridPattern(const Tree& tree, size_t n_row, size_t n_col, size_t which);
	SOPVectorcd CNotGridSwipe(const Tree& tree, size_t n_row, size_t n_col);

	// This is a light version of the google-benchmark resulting in a thomas-porter distribution
	SOPVectorcd GoogleIterations(const Tree& tree, size_t n_row, size_t n_col, size_t n_swipes);

	// Convert indices
	size_t ConvertIdx2D(size_t idx, size_t nrow, size_t ncol);
}

#endif //MCTDH_QUBIT_OPERATORS_H
