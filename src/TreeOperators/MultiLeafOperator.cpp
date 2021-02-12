#include "TreeOperators/MultiLeafOperator.h"


template <typename T>
MultiLeafOperator<T>::MultiLeafOperator()
	: hasV_(false) {
	targetLeaves_.clear();
}

template <typename T>
Tensor<T> MultiLeafOperator<T>::apply(Tensor<T> Phi,
	const Leaf& leaf) const {
	Tensor<T> hPhi(Phi.shape());
	size_t mode_x = leaf.Mode();
	const LeafInterface& grid = leaf.PrimitiveGrid();
	bool switchbool = true;

	// Applying the MLO uses switching of the result Tensor to increase performance.
	for (size_t l = 0; l < leafOperators_.size(); ++l) {
		if (mode_x != targetLeaves_[l]) { continue; }

		shared_ptr<LeafOperator<T>> spo = operator[](l);
		// apply it
		if (switchbool) {
			spo->apply(grid, hPhi, Phi);
		} else {
			spo->apply(grid, Phi, hPhi);
		}
		switchbool = !switchbool;
	}

	if (switchbool) {
		return Phi;
	} else {
		return hPhi;
	}
}

template <typename T>
Tensor<T> MultiLeafOperator<T>::apply(Tensor<T> Acoeffs,
	const vector<int>& list, const LeafInterface& grid) const {
	Tensor<T> hAcoeff(Acoeffs.shape());
	bool switchbool = true;
	// Applying the MLO uses switching of the result Tensor to increase performance.
	for (size_t l = 0; l < list.size(); l++) {
		// get the active_ part in the MLO
		int part = list[l];
		shared_ptr<LeafOperator<T>> spo = operator[](part);

		// apply it
		if (switchbool) {
			spo->apply(grid, hAcoeff, Acoeffs);
		} else {
			spo->apply(grid, Acoeffs, hAcoeff);
		}
		switchbool = !switchbool;
	}

	if (switchbool) {
		return Acoeffs;
	} else {
		return hAcoeff;
	}
}

template<typename T>
TensorTree<T> MultiLeafOperator<T>::apply(TensorTree<T> Psi,
	const Tree& tree) const {
	applyReference(Psi, tree);
	return Psi;
}

template<typename T>
void MultiLeafOperator<T>::applyReference(TensorTree<T>& Psi, const Tree& tree) const {
	for (const Node& node : tree) {
		if (node.isBottomlayer()) {
			const Leaf& phy = node.getLeaf();
			const LeafInterface& grid = phy.PrimitiveGrid();
			size_t coord = phy.Mode();

			// build list with active_ parts
			vector<int> activelayerparts;
			for (size_t l = 0; l < leafOperators_.size(); l++) {
				if (targetLeaves_[l] == coord) {
					activelayerparts.push_back(l);
				}
			}

			Tensor<T>& Acoeff = Psi[node];
			Acoeff = apply(Acoeff, activelayerparts, grid);
		}
	}
}


template <typename T>
void MultiLeafOperator<T>::setV(const PotentialOperator& V_) {
    v_ = V_;
    hasV_ = true;
}

template <typename T>
MultiLeafOperator<T> operator*(const MultiLeafOperator<T>& A,
	const MultiLeafOperator<T>& B) {
	MultiLeafOperator<T> M = B;

	for (size_t i = 0; i < A.size(); i++) {
		M.push_back(A[i], A.mode(i));
	}

	return M;
}


typedef complex<double> cd;
template class MultiLeafOperator<complex<double>>;
template MultiLeafOperator<cd> operator*(const MultiLeafOperator<cd>& A, const MultiLeafOperator<cd>& B);

template class MultiLeafOperator<double>;
