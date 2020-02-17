//
// Created by Roman Ellerbrock on 2/1/20.
//
#include "SpectralDecompositionTree.h"

template<typename T>
SpectralDecompositionTree<T>::SpectralDecompositionTree(const Tree& tree) {
	Initialize(tree);
}

template<typename T>
SpectralDecompositionTree<T>::SpectralDecompositionTree(const MatrixTree<T>& H,
	const Tree& tree) {
	Initialize(tree);
	Calculate(H, tree);
}

template<typename T>
void SpectralDecompositionTree<T>::Initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorDim& tdim = node.TDim();
		size_t dim = tdim.LastActive();
		auto x = SpectralDecomposition<T>(Matrix<T>(dim, dim), Vectord(dim));
		attributes_.emplace_back(x);
	}
}

template<typename T>
void SpectralDecompositionTree<T>::Calculate(const MatrixTree<T>& H,
	const Tree& tree) {
	for (const Node& node : tree) {
		const Matrix<T>& mat = H[node];
		auto y = Diagonalize(mat);
		this->operator[](node) = y;
	}
}

template<typename T>
MatrixTree<T> SpectralDecompositionTree<T>::Invert(const Tree& tree, double eps) {
	assert(attributes_.size() == tree.nNodes());
	MatrixTree<T> Inv_Hole(tree);
	for (const Node& node : tree) {
		const auto& x = this->operator[](node);
		Inv_Hole[node] = BuildInverse(x, eps);
	}
	return Inv_Hole;
}

template<typename T>
void SpectralDecompositionTree<T>::print(const Tree& tree) const {
	for (auto it = tree.rbegin(); it != tree.rend(); it++) {
		const Node& node = *it;
		if (!node.IsToplayer()) {
			node.info();
			const SpectralDecomposition<T>& x = this->operator[](node);
			x.second.print();
		}
	}
}

template
class SpectralDecompositionTree<complex<double>>;

template
class SpectralDecompositionTree<double>;
