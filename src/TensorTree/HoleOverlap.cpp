//
// Created by Roman on 3/29/2019.
//

#include "HoleOverlap.h"

template<typename T>
HoleOverlap<T>::HoleOverlap(const TTBasis& basis){
	Initialize(basis);
}

template<typename T>
HoleOverlap<T>::HoleOverlap(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const DenseOverlap<T>& S,
		const TTBasis& basis){
	Initialize(basis);
	Calculate(Psi1, Psi2, S, basis);
}

template<typename T>
void HoleOverlap<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		size_t n = tdim.getntensor();
		attributes.emplace_back(Matrix<T>(n, n));
	}
}

template<typename T>
void HoleOverlap<T>::CalculateLayer(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const DenseOverlap<T>& S, const Node& node) {

	const Node& parent = node.Up();
	const Tensor<T>& Bra = Psi[parent];
	Tensor<T> Ket = Chi[parent];

	// Transform sublying basis
	auto child_idx = (size_t) node.ChildIdx();
	for (size_t k = 0; k < parent.nChildren(); ++k) {
		if (k != child_idx) {
			const Node& child = parent.Down(k);
			const FactorMatrix<T>& h = S[child];
			Ket = h * Ket;
		}
	}

	// Transform with upper hole matrix
	if (!parent.IsToplayer()) {
		const Matrix<T>& holemat = this->operator[](parent);
		Ket = multStateAB(holemat, Ket);
	}

	// Calculate hole-product and save to attributes
	this->operator[](node) = HoleProduct(Bra, Ket, child_idx);
}

template<typename T>
void HoleOverlap<T>::Calculate(const TensorTree<T>& Psi1,
		const TensorTree<T>& Psi2, const DenseOverlap<T>& S,
		const TTBasis& basis) {
	// Calculate the holematrices for two different wavefunctions.
	// @TODO: This routines offers many ways for improvement. Go ahead!

	// Top-down_ swipe. Note that the top-node is excluded from swipe!
	for (int n = (int) basis.nNodes() - 2; n >= 0; --n) {
		const Node& node = basis.GetNode(n);
		CalculateLayer(Psi1, Psi2, S, node);
	}
}

template <typename T>
void HoleOverlap<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		node.info();
		(*this)[node].print(os);
	}
}

template class HoleOverlap<complex<double>>;
template class HoleOverlap<double>;

