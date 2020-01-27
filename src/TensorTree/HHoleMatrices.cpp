#include "HHoleMatrices.h"

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
		Tensorcd hKet = hmat.ApplyHole(Ket[parent], node);
		if (!parent.IsToplayer()) {
			hKet = multStateAB(operator[](parent), hKet);
		}
		operator[](node) = HoleProduct(Bra[parent], hKet, node.ChildIdx());
	}
}

template<typename T>
Tensorcd HHoleMatrices<T>::Apply(const Tensorcd& Phi, const Node& node) const {
	if (node.IsToplayer()) {
		return Phi;
	} else {
		assert(Active(node));
		return multStateAB(operator[](node), Phi);
	}
}

template<typename T>
void HHoleMatrices<T>::print(TTBasis& basis, ostream& os) {
	for (const Node& node : basis) {
		if (!node.IsToplayer()) {
			node.info(os);
			operator[](node).print(os);
		}
	}
}


