//
// Created by Roman Ellerbrock on 2019-04-19.
//

#ifndef MCTDH_QUBIT_OPERATORS_H
#define MCTDH_QUBIT_OPERATORS_H

#include "Hamiltonian.h"
#include "QMConstants.h"
#include "Register.h"

namespace GateOperators {

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
	void Hadamard(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void X(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Y(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Z(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void S(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void T(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Sdagger(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Set0(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Set1(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void Identity(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
//	void SetBlack(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void sqrtX(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);
	void sqrtY(const PrimitiveBasis& grid, Tensorcd& xPhi, const Tensorcd& phi);


	class PrimitiveProjector: public BottomLayerSPO {
	public:
		explicit PrimitiveProjector(size_t idx)	: idx_(idx) {}
		~PrimitiveProjector() = default;
		void Apply(Tensorcd& PPhi, const PrimitiveBasis& grid, const Tensorcd& Phi) const override {
			const TensorDim& tdim = Phi.Dim();
			PPhi.Zero();
			for (size_t n = 0; n < tdim.getntensor(); ++n) {
				PPhi(idx_, n) = Phi(idx_, n);
			}
		}
	private:
		size_t idx_;
	};

	class Rk: public BottomLayerSPO {
		/// This class represents the R(k) Phase shift operation
	public:
		/// The phase will be converted to a factor: exp(2*pi*im*phase_)
		explicit Rk(size_t n, bool adjungate = false);
		~Rk() = default;
		/// Efficient applying routine
		void Apply(Tensorcd& RAcoeff, const PrimitiveBasis& grid, const Tensorcd& Acoeff) const override;
	private:
		complex<double> factor;
	};

	class Rx: public BottomLayerSPO {
		/// This class represents the R(k) Phase shift operation
	public:
		/// The phase will be converted to a factor: exp(2*pi*im*phase_)
		explicit Rx(long double n, bool adjungate = false);
		~Rx() = default;
		/// Efficient applying routine
		void Apply(Tensorcd& RAcoeff, const PrimitiveBasis& grid, const Tensorcd& Acoeff) const override;
	private:
		complex<double> factor;
	};

	class SetBlack: public BottomLayerSPO {
	public:
		SetBlack(double epsilon_, bool adjungate_)
			: epsilon(epsilon_), adjungate(adjungate_) {}
		~SetBlack() = default;
		void Apply(Tensorcd& xPhi, const PrimitiveBasis& grid,
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
	SOP MakeCGate(const SOP& sop, size_t control);

	/// Overloads purely for typecasting convenience
	inline SOP MakeCGate(const MPO& mpo, size_t control) { return MakeCGate(SOP(mpo), control); }

	inline SOP MakeCGate(const RefSPO& spo, size_t control, size_t target) {
		return MakeCGate(MPO(spo, target),
			control);
	}

	inline SOP MakeCGate(shared_ptr<BottomLayerSPO> spo, size_t control, size_t target) {
		MPO M;
		M.push_back(spo, target);
		return MakeCGate(M, control);
	}

	/// Build double-control gates (e.g. ccX == CCNOT == Toffoli gate)
	SOP MakeCCGate(const SOP& sop, size_t control1, size_t control2);

	/// Overloads purely for typecasting convenience
	inline SOP MakeCCGate(const MPO& mpo, size_t control1, size_t control2) {
		return MakeCCGate(SOP(mpo),
			control1,
			control2);
	}

	inline SOP MakeCCGate(const RefSPO& spo, size_t control1, size_t control2, size_t target) {
		return MakeCCGate(MPO(spo,
			target), control1, control2);
	}

	inline SOP MakeCCGate(shared_ptr<BottomLayerSPO> spo, size_t control1, size_t control2, size_t target) {
		MPO M;
		M.push_back(spo, target);
		return MakeCCGate(M, control1, control2);
	}

	/// Convenience functions
	SOP CNot(const Register& r1, const Register& r2);
	SOP CNot(size_t control, size_t target);
	SOP Uk(const Register& c, const Register& t, size_t k, bool conj_phase = false);
	SOP Uk(size_t control, size_t target, size_t k, bool conj_phase = false);
	SOP Toffoli(const Register& c1, const Register& c2, const Register& t);
	SOP Toffoli(size_t control1, size_t control2, size_t target);

	SOPVector DistribSingleControl(const SOPVector& stack, size_t control);
	SOPVector DistribDoubleControl(const SOPVector& stack, size_t control1, size_t control2);

	/// =========================
	/// High-level operations
	/// =========================

	/// Building blocks for logic circuits
	SOPVector Carry(const Register& c1, const Register& a, const Register& b, const Register& c2);
	SOPVector Carry(size_t carry1, size_t a, size_t b, size_t carry2);
	SOPVector RCarry(const Register& c1, const Register& a, const Register& b, const Register& c2);
	SOPVector RCarry(size_t carry1, size_t a, size_t b, size_t carry2);
	SOPVector SUM(const Register& a, const Register& b, const Register& sum);
	SOPVector SUM(size_t a, size_t b, size_t sum);

	SOPVector Adder_circuit(const mctdhBasis& basis, size_t n_bit, size_t mode_a, size_t mode_b,
		size_t mode_carry);

	SOPVector RandomSuperposition(const mctdhBasis& basis, size_t start, size_t end, size_t n);
	SOPVector ControlledSuperposition(const Register& reg, size_t n);
	SOPVector HadamardChain(const Register& reg);
	SOPVector Distribute(const Register& reg, function<void(const PrimitiveBasis&, Tensorcd&, const Tensorcd&)> f);

	// This is the CNot-Chain-benchmark
	SOPVector CNotChain(const Register& reg);
	SOPVector CNotBenchmark(const Register& reg);
	SOPVector CNotGridRow(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t offset = 0, size_t shift = 0);
	SOPVector CNotGridCol(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t offset = 0, size_t shift = 0);
	SOPVector CNotGridPattern(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t which);
	SOPVector CNotGridSwipe(const mctdhBasis& basis, size_t n_row, size_t n_col);

	// This is a light version of the google-benchmark resulting in a thomas-porter distribution
	SOPVector GoogleIterations(const mctdhBasis& basis, size_t n_row, size_t n_col, size_t n_swipes);

	// Convert indices
	size_t ConvertIdx2D(size_t idx, size_t nrow, size_t ncol);
}

#endif //MCTDH_QUBIT_OPERATORS_H
