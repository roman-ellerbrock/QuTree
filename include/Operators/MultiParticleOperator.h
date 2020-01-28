#pragma once
#include "Core/Tensor.h"
#include "TensorTreeBasis.h"
#include "TensorTree.h"
#include "SingleParticleOperator.h"
#include "SingleParticleOperatorFunction.h"
#include "SingleParticleOperatorMatrix.h"
#include "PotentialOperator.h"

template <typename T>
class MultiParticleOperator
	/*!
	 * \brief A MultiParticleOperator (MPO) is a product of general single-particle operators
	 *
	 * A MultiParticleOperator is one summand in a SumofProducts
	 * operator. It can be applied to a mctdhWavefunction resulting in
	 * a wavefunction with a different SPF-basis. The MPO is central for building hamiltonians.
	 *
	 * Please note:
	 * - The recommended operator to build MPO is RefSPO (due to performance)
	 * - SPO initialization might be removed in the future
	 * - mctdhWavfunctions with different SPF-basis
	 * functions cannot be added without loss of information.
	 */
{
public:
	/// Constructor without memory allocation
	MultiParticleOperator();

	/// Default destructor
	~MultiParticleOperator() = default;

	/// Construct a MPO from a single SPO
	MultiParticleOperator(const SPOM<T>& h, int mode_x)
		: MultiParticleOperator() {
		push_back(h, mode_x);
	}

	/// Construct a MPO from a single RefSPO
	MultiParticleOperator(const SPOf<T>& h, int mode_x)
		: MultiParticleOperator() {
		push_back(h, mode_x);
	}

	/// This routine manages how to apply a MPO
	Tensor<T> ApplyBottomLayer(Tensor<T> Acoeffs,
		const Leaf& phys) const;

	/// This routine manages how to apply a MPO
	Tensor<T> ApplyBottomLayer(Tensor<T> Acoeffs,
		const vector<int>& list, const PrimitiveBasis& grid) const;

	/// Push back a SPO to the MPO
	void push_back(shared_ptr<SPO<T>> h, int mode_x) {
		SingParOp.push_back(h);
		mode_.push_back(mode_x);
	}

	void push_back(const SPOM<T>& h, size_t mode) {
		auto* h_ptr = new SPOM<T>(h);
		h_ptr->Mode() = 0;
		SingParOp.push_back(shared_ptr<SPO<T>>(h_ptr));
		mode_.emplace_back(mode);
	}

	/// Push back a RefSPO to the MPO
	void push_back(const SPOf<T>& h, int mode_x) {
		auto *spo = new SPOf<T>(h);
		SingParOp.push_back(shared_ptr<SPO<T>>(spo));
		mode_.push_back(mode_x);
	}

	/// Access the i-th SPO in the MPO
	shared_ptr<SPO<T>> operator()(size_t i) {
		assert(i >= 0);
		assert(i < SingParOp.size());
		return SingParOp[i];
	}

	/// The number of SPOs in the MPO
	size_t size() const { return mode_.size(); }

	/// Access the i-th SPO in the MPO
	shared_ptr<SPO<T>> operator[](size_t i) const {
		assert(i < mode_.size());
		assert(i >= 0);
		return SingParOp[i];
	}

	/// Check whether a SPO acts on mode "mode_x"
	bool ModeIsActive(size_t mode_x) const {
		for (size_t i = 0; i < mode_.size(); i++) {
			if (mode_[i] == mode_x) { return true; }
		}
		return false;
	}

	/// Check whether the "part"-th SPO acts on mode "mode_x"
	bool ModeIsActive(size_t part, size_t mode_x) const {
		return (mode_[part] == mode_x);
	}

	/// Apply a MPO to a wavefunction. This routine is not optimized for performance.
	TensorTree<T> Apply(TensorTree<T> Psi, const TTBasis& basis) const;

	/// On which mode does the "part"-th SPO act?
	size_t Mode(size_t part) const {
		assert(part < size());
		return mode_[part];
	}

	/// Multiply two MPOs.
	friend MultiParticleOperator operator*(const MultiParticleOperator& A,
		const MultiParticleOperator& B) {
		MultiParticleOperator M = B;

		for (size_t i = 0; i < A.size(); i++) {
			M.push_back(A[i], A.Mode(i));
		}

		return M;
	}

	/// Check whether a Potential operator is included in the MPO. Important for CDVR.
	bool HasV() const { return hasV; }

	/// Get a reference to the PotentialOperator
	PotentialOperator& V() {
		return v;
	}

	/// Get a reference to the PotentialOperator
	const PotentialOperator& V() const {
		return v;
	}

	/// Set the PotentialOperator
	void SetV(const PotentialOperator& V);

	/// Return vector of all active modes in this operator
	const vector<int>& Modes()const { return mode_; }

protected:
	/// These are the SPOs
	vector<shared_ptr<SPO<T>>> SingParOp;
	/// These are the modes the SPOs act on
	vector<int> mode_;
	/// The potential operator
	PotentialOperator v;
	/// Is there a PotentialOperator?
	bool hasV;
};

template <typename T>
using MPO = MultiParticleOperator<T>;

typedef MPO<complex<double>> MPOcd;

