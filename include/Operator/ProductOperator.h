#pragma once
#include "Tensor/Tensor"
#include "Tree/Tree.h"
#include "TensorNetwork/TensorTree.h"

#include "LeafOperator.h"
#include "LeafFunction.h"
#include "LeafMatrix.h"
#include "PotentialOperator.h"

template<typename T>
class ProductOperator
	/*!
	 * \class MultiLeafOperator
	 * \ingroup Operators
	 * \brief A MultiLeafOperator (MLO) is a product of general single-particle operators
	 *
	 * A MultiLeafOperator is one summand in a SumofProducts
	 * operator. It can be applied to a TensorTree resulting in
	 * a wavefunction with a different bottomlayer-Tensors. The MLO is
	 * central for building hamiltonians.
	 *
	 * Please note:
	 * - The recommended operator to build MLO is SPOf or SPOM (due to performance)
	 * - SPO initialization might be removed in the future
	 * - TensorTrees with different lower-node Tensors
	 *   cannot (straightforwardly) be added without loss of information.
	 */
{
public:
	/// Constructor without memory allocation
	ProductOperator() : hasV_(false) {}

	/// Default destructor
	~ProductOperator() = default;

	ProductOperator(const shared_ptr<LeafOperator<T>>& h, size_t leaf_id)
		: hasV_(false) {
		push_back(h, leaf_id);
	}

	/// Construct a MLO from a single SPO
	ProductOperator(const LeafMatrix<T>& h, size_t leaf_id)
		: ProductOperator() {
		push_back(h, leaf_id);
	}

	ProductOperator(const Matrix<T> h, size_t leaf_id, size_t adjoint = false)
		: ProductOperator() {
		push_back(h, leaf_id, adjoint);
	}

	/// Construct a MLO from a single RefSPO
	ProductOperator(const LeafFunction<T>& h, size_t leaf_id)
		: ProductOperator() {
		push_back(h, leaf_id);
	}

	/// Construct a MLO from a single RefSPO
	ProductOperator(const LeafFun<T>& h, size_t leaf_id)
		: ProductOperator() {
		push_back(h, leaf_id);
	}

	/// Push back a SPO to the MLO
	void push_back(shared_ptr<LeafOperator<T>> h, size_t mode_x) {
		leafOperators_.push_back(h);
		targetLeaves_.push_back(mode_x);
	}

	void push_back(const LeafMatrix<T>& h, size_t mode) {
		push_back(make_shared<LeafMatrix<T>>(h), mode);
	}

	void push_back(const Matrix<T>& h, size_t mode, bool adj = false) {
		push_back(make_shared<LeafMatrix<T>>(h, adj), mode);
	}

	/// Push back a RefSPO to the MLO
	void push_back(const LeafFunction<T>& h, size_t mode) {
		push_back(make_shared<LeafFunction<T>>(h), mode);
	}

	void push_back(const LeafFun<T>& h, size_t mode) {
		push_back(make_shared<LeafFunction<T>>(h), mode);
	}

	void push_back(const LeafFunPair<T>& hs, size_t mode, bool adjoint = false) {
		if (adjoint) {
			push_back(hs.second, mode);
		} else {
			push_back(hs.first, mode);
		}
	}

	/// The number of SPOs in the MLO
	[[nodiscard]] size_t size() const { return targetLeaves_.size(); }

	/// Access the i-th SPO in the MLO
	shared_ptr<LeafOperator<T>> operator[](size_t i) const {
		assert(i < targetLeaves_.size());
		return leafOperators_[i];
	}

	/// apply a MLO to a wavefunction.
	void applyReference(TensorTree<T>& Psi) const;

	/// apply a MLO to a wavefunction. Call-by-value-version. This routine is not optimized for performance.
	TensorTree<T> apply(TensorTree<T> Psi) const;

	/// On which mode does the "part"-th SPO act?
	[[nodiscard]] size_t mode(size_t part) const {
		assert(part < size());
		return targetLeaves_[part];
	}

	/// Multiply two MPOs.
	ProductOperator operator*(const ProductOperator& B) const {
		ProductOperator<T> M = B;

		for (size_t i = 0; i < leafOperators_.size(); i++) {
			M.push_back(leafOperators_[i], targetLeaves_[i]);
		}

		return M;
	}

	/// This routine manages how to apply a MLO
	//	Tensor<T> apply(Tensor<T> A,
	//		const Leaf& leaf) const;

	/*	/// Check whether a SPO acts on mode "mode_x"
		bool isActive(size_t mode_x) const {
			for (size_t i = 0; i < targetLeaves_.size(); i++) {
				if (targetLeaves_[i] == mode_x) { return true; }
			}
			return false;
		}

		/// Check whether the "part"-th SPO acts on mode "mode_x"
		bool isActive(size_t part, size_t mode_x) const {
			return (targetLeaves_[part] == mode_x);
		}
	*/
/*	/// Check whether a Potential operator is included in the MLO. Important for CDVR.
	[[nodiscard]] bool hasV() const { return hasV_; }

	/// Get a reference to the PotentialOperator
	PotentialOperator& v() {
		return v_;
	}

	/// Get a reference to the PotentialOperator
	[[nodiscard]] const PotentialOperator& v() const {
		return v_;
	}

	/// Set the PotentialOperator
	void setV(const PotentialOperator& V);
*/
	/// Return vector of all active_ modes in this operator
//	[[nodiscard]] const vector<size_t>& targetLeaves() const { return targetLeaves_; }

	void print(ostream& os = cout) const {
		os << size() << " operators in MLO" << endl;
	}

	[[nodiscard]] const vector<size_t>& targetLeaves() const { return targetLeaves_; }

protected:
	void applyFactor(TensorTree<T>& Psi, size_t i) const;
	/// These are the SPOs
	vector<shared_ptr<LeafOperator<T>>> leafOperators_;
	/// These are the modes the SPOs act on
	vector<size_t> targetLeaves_;
	/// The potential operator
//	PotentialOperator v_;
	/// Is there a PotentialOperator?
	bool hasV_;
};

template<typename T>
using PO = ProductOperator<T>;

typedef PO<complex<double>> POcd;

typedef PO<double> POd;

