//
// Created by Roman Ellerbrock on 2/1/20.
//
#include "SpectralDecompositionTree.h"

template <typename T>
SpectralDecompositionTree<T>::SpectralDecompositionTree(const TTBasis& basis) {
	Initialize(basis);
}

template<typename T>
SpectralDecompositionTree<T>::SpectralDecompositionTree(const HoleMatrixTree<T>& H,
	const TTBasis& basis) {
	Initialize(basis);
	Calculate(H, basis);
}

template <typename T>
void SpectralDecompositionTree<T>::Initialize(const TTBasis & basis) {
	attributes.clear();
	for (const Node& node : basis) {
		const TensorDim& tdim = node.TDim();
		size_t dim = tdim.GetNumTensor();
		auto x = SpectralDecomposition<T>(Matrix<T>(dim, dim), Vectord(dim));
		attributes.emplace_back(x);
	}
}

template <typename T>
void SpectralDecompositionTree<T>::Calculate(const HoleMatrixTree<T> & H,
	const TTBasis & basis) {
	for (const Node& node : basis) {
		const Matrix<T>& mat = H[node];
		auto y = Diagonalize(mat);
		this->operator[](node) = y;
	}
}

template<typename T>
HoleMatrixTree<T> SpectralDecompositionTree<T>::Invert(const TTBasis& basis, double eps) {
	assert(attributes.size() == basis.nNodes());
	HoleMatrixTree<T> Inv_Hole(basis);
	for (const Node& node : basis) {
		const auto& x = this->operator[](node);
		const Vectord& vec = x.second;
//		auto inv_ev = Inverse(vec, eps);
		SpectralDecomposition<T> y(x.first, x.second);
		Inv_Hole[node] = BuildMatrix(y);
	}
	return Inv_Hole;
}

template <typename T>
void SpectralDecompositionTree<T>::print(const TTBasis& basis) const {
	for (const Node& node : basis) {
		node.info();
		const SpectralDecomposition<T>& x = this->operator[](node);
		x.second.print();
	}
}

template class SpectralDecompositionTree<complex<double>>;
template class SpectralDecompositionTree<double>;

