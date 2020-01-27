#include "HMatrices.h"


template<typename T>
void HMatrices<T>::Calculate(const TensorTree<T>& Bra, const TensorTree<T>& Ket,
	const MPO<T>& M, const TTBasis& basis) {
	for (size_t n = 0; n < Active().size(); ++n) {
		const Node& node = Active().MCTDHNode(n);
		CalculateLayer(Bra[node], Ket[node], M, node);
	}
}

template<typename T>
SPO<T> HMatrices<T>::CalculateUpper(const Tensor<T>& Bra, const Tensor<T>& Ket,
	const Node& node) {
	// @TODO: Optimize with switchbool trick
	// Swipe through children and apply active children's SPOs.
	Tensor<T> hKet(Ket);
	for (size_t l = 0; l < node.nChildren(); l++) {
		const Node& child = node.Down(l);
		if (!Active(child)) { continue; }
		const SPO<T>& h = operator[](child);
		hKet = h * hKet;
	}

	// calculate overlap
	size_t ChildIdx = node.ChildIdx();
	Matrixcd resultmat = Bra.DotProduct(hKet);
	return SPO<T>(resultmat, ChildIdx);
}

template<typename T>
SPO<T> HMatrices<T>::CalculateBottom(const Tensor<T>& Bra, const Tensor<T>& Ket,
	const MPO<T>& M, const Node& node,
	const Leaf& phys) {

	Tensor<T> MKet = M.ApplyBottomLayer(Ket, phys);
	int ChildIdx = node.ChildIdx();
	Matrixcd resultmat = Bra.DotProduct(MKet);
	return SPO<T>(resultmat, ChildIdx);
}

template<typename T>
void HMatrices<T>::CalculateLayer(const Tensor<T>& Bra, const Tensor<T>& Ket,
	const MPO<T>& M, const Node& node) {
	if (!Active(node)) { return; }

	if (node.IsBottomlayer()) {
		operator[](node) = CalculateBottom(Bra, Ket, M, node, node.PhysCoord());
	} else {
		operator[](node) = CalculateUpper(Bra, Ket, node);
	}
}

template<typename T>
Tensor<T> HMatrices<T>::Apply(const Tensor<T>& Phi, const MPO<T>& M,
	const Node& node) const {
	if (!Active(node)) { return Phi; }
	if (node.IsBottomlayer()) {
		const Leaf& phys = node.PhysCoord();
		return M.ApplyBottomLayer(Phi, phys);
	} else {
		return ApplyUpper(Phi, node);
	}
}

template<typename T>
Tensor<T> HMatrices<T>::ApplyUpper(Tensor<T> Phi, const Node& node) const {
	Tensor<T> hPhi(Phi.Dim());
	bool switchbool = true;
	for (size_t k = 0; k < node.nChildren(); ++k) {
		const Node& child = node.Down(k);
		if (!Active(child)) { continue; }
		const SPO<T>& h = operator[](child);
		if (switchbool) {
			multAB(hPhi, h, Phi, true);
		} else {
			multAB(Phi, h, hPhi, true);
		}
		switchbool = !switchbool;
	}
	if (switchbool) {
		return Phi;
	} else {
		return hPhi;
	}
}

template<typename T>
Tensor<T> HMatrices<T>::ApplyHole(Tensor<T> Phi, const Node& hole_node) const {
	// @TODO: Optimize with switchbool trick
	assert(!hole_node.IsToplayer());
	const Node& parent = hole_node.Up();
	size_t drop = hole_node.ChildIdx();

	for (size_t k = 0; k < parent.nChildren(); ++k) {
		const Node& child = parent.Down(k);
		if ((child.ChildIdx() == drop) || (!Active(child))) { continue; }
		const SPO<T>& h = operator[](child);
		const TensorDim& tdim = Phi.Dim();
		Phi = h * Phi;
	}
	return Phi;
}

template<typename T>
void HMatrices<T>::print(const TTBasis& basis, ostream& os) const {
	for (const Node& node : basis) {
		if (!Active(node)) { continue; }
		node.info(os);
	}
}

