#include "DenseOverlap.h"

template<typename T>
DenseOverlap<T>::DenseOverlap(const TTBasis& basis) {
	Initialize(basis);
}

template<typename T>
DenseOverlap<T>::DenseOverlap(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {
	Initialize(basis);
	Calculate(Psi, Chi, basis);
}

template<typename T>
void DenseOverlap<T>::Initialize(const TTBasis& basis) {
	// Clear the overlaps for reinitialization
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		const size_t dim = tdim.getntensor();
		const auto k = (size_t) node.ChildIdx();
		FactorMatrix<T> mat(dim, k);
		attributes.push_back(mat);
	}
}

template<typename T>
FactorMatrix<T> DenseOverlap<T>::Calculate(const TensorTree<T>& Psi,
	const TensorTree<T>& Chi, const TTBasis& basis) {

	for (const Node& node : basis) {
		CalculateLayer(Psi[node], Chi[node], node);
	}

	const Node& topnode = basis.TopNode();
	return this->operator[](topnode);
}

template<typename T>
void DenseOverlap<T>::CalculateLayer(const Tensor<T>& Phi,
	Tensor<T> AChi, const Node& node) {
	// Get references to the ACoefficients at each node
	if (!node.IsBottomlayer()) {
		for (int k = 0; k < node.nChildren(); k++) {
			// Get overlap-matrix from down_ under
			const Node& child = node.Down(k);
			const FactorMatrix<T>& spo = this->operator[](child);
			// Apply it to the right-hand side
			AChi = spo * AChi;
		}
	}
	int kchild = node.ChildIdx();
	Matrix<T> resultmat = Phi.DotProduct(AChi);
	this->operator[](node) = FactorMatrix<T>(resultmat, kchild);
}

template<typename T>
Tensor<T> DenseOverlap<T>::TransformTensor(const Tensor<T>& Phi,
	const Node& node) const {
	if (!node.IsBottomlayer()) {
		Tensor<T> hPhi(Phi);
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.Down(k);
			const FactorMatrix<T>& s = this->operator[](child);
			hPhi = s * hPhi;
		}
		return hPhi;
	} else {
		return Phi;
	}
}

template class DenseOverlap<complex<double>>;
template class DenseOverlap<double>;


