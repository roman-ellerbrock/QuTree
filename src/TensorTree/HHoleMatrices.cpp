#include "HHoleMatrices.h"

template<typename T>
HHoleMatrices<T>::HHoleMatrices(const TensorTree<T>& Psi,
	const HMatrices<T>& hmat, const MPO<T>& M, const TTBasis& basis)
	: HHoleMatrices(M, basis) {
	Calculate(Psi, Psi, hmat, basis);
}

template<typename T>
void HHoleMatrices<T>::Initialize(const TTBasis& basis) {
	attributes.clear();
	for (const Node*const node_ptr : Active()) {
		const Node& node = *node_ptr;
		size_t dim = node.TDim().getntensor();
		attributes.emplace_back(Matrix<T>(dim, dim));
	}
}

// Calculate Hole-Matrices
template<typename T>
void HHoleMatrices<T>::Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
	const HMatrices<T>& hmat, const TTBasis& basis) {

	// Swipe top-down_ but exclude topnode
	int sub_topnode = Active().size() - 2;
	for (int n = sub_topnode; n >= 0; --n) {
		const Node& node = Active().MCTDHNode(n);
		assert(Active(node));

		const Node& parent = node.Up();
		Tensor<T> hKet = hmat.ApplyHole(Ket[parent], node);
		if (!parent.IsToplayer()) {
			hKet = multStateAB(this->operator[](parent), hKet);
		}
		operator[](node) = HoleProduct(Bra[parent], hKet, node.ChildIdx());
	}
}

template<typename T>
Tensor<T> HHoleMatrices<T>::Apply(const Tensor<T>& Phi, const Node& node) const {
	if (node.IsToplayer()) {
		return Phi;
	} else {
		assert(Active(node));
		return multStateAB(this->operator[](node), Phi);
	}
}

template<typename T>
void HHoleMatrices<T>::print(TTBasis& basis, ostream& os) {
	for (const Node& node : basis) {
		if (!node.IsToplayer()) {
			node.info(os);
			this->operator[](node).print(os);
		}
	}
}

template class HHoleMatrices<complex<double>>;
template class HHoleMatrices<double>;

