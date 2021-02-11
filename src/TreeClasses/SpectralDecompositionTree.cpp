//
// Created by Roman Ellerbrock on 2/1/20.
//
#include "TreeClasses/SpectralDecompositionTree.h"

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
		const TensorShape& tdim = node.shape();
		size_t dim = tdim.lastDimension();
		auto x = SpectralDecomposition<T>(Matrix<T>(dim, dim), Vectord(dim));
		attributes_.emplace_back(x);
	}
}

template<typename T>
void SpectralDecompositionTree<T>::Calculate(const MatrixTree<T>& H,
	const Tree& tree) {
	for (const Node& node : tree) {
		const Matrix<T>& mat = H[node];
		auto y = diagonalize(mat);
		this->operator[](node) = y;
	}
}

template<typename T>
MatrixTree<T> SpectralDecompositionTree<T>::Invert(const Tree& tree, double eps) {
	assert(attributes_.size() == tree.nNodes());
	MatrixTree<T> Inv_Hole(tree);
	for (const Node& node : tree) {
		const auto& x = this->operator[](node);
		Inv_Hole[node] = toMatrix(inverse(x, eps));
	}
	return Inv_Hole;
}

template<typename T>
void SpectralDecompositionTree<T>::print(const Tree& tree) const {
	for (auto it = tree.rbegin(); it != tree.rend(); it++) {
		const Node& node = *it;
		if (!node.isToplayer()) {
			node.info();
			const SpectralDecomposition<T>& x = this->operator[](node);
			auto ev = x.second;
			double norm = 0.;
			for (size_t i = 0; i < ev.dim(); ++i) {
				norm += ev(i);
			}
			norm = max(1e-50, norm);
			ev /= norm;
			ev.print();
		}
	}
}

template<typename T>
MatrixTree<T> to_matrixtree(const SpectralDecompositionTree<T>& X, const Tree& tree) {
	MatrixTree<T> mattree(tree);
	for (const Node& node : tree) {
		mattree[node] = toMatrix(X[node]);
	}
	return mattree;
}

template<typename T>
void CanonicalTransformation(TensorTree<T>& Psi, const Tree& tree, bool orthogonal) {

	auto rho = TreeFunctions::Contraction(Psi, tree, orthogonal);
	SpectralDecompositionTree<T> spec(rho, tree);

	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			Matrix<T> p = spec[node].first;
			p = p.adjoint();
			Psi[node] = matrixTensor(p, Psi[node], node.nChildren());
//			Psi[node] = multATB(p, Psi[node], node.nChildren());
			const Node& parent = node.parent();
			Psi[parent] = tensorMatrix(Psi[parent], p.adjoint(), node.childIdx());
//			Psi[parent] = multATB(p, Psi[parent], node.childIdx());
		}
	}
}

template<typename T>
SpectralDecompositionTree<T> sqrt(SpectralDecompositionTree<T> X) {
	for (auto& x : X) {
		x = sqrt(x);
	}
	return X;
}

template<typename T>
MatrixTree<T> sqrt(MatrixTree<T> X, const Tree& tree) {
	SpectralDecompositionTree<T> x(X, tree);
	x = sqrt(x);
	return to_matrixtree(x, tree);
}

template<typename T>
SpectralDecompositionTree<T> inverse(SpectralDecompositionTree<T> X, double eps) {
	for (SpectralDecomposition<T>& x : X) {
		x = inverse(x, eps);
	}
	return X;
}

template<typename T>
MatrixTree<T> inverse(MatrixTree<T> X, const Tree& tree, double eps) {
	SpectralDecompositionTree<T> x(X, tree);
	x = inverse(x, eps);
	return to_matrixtree(x, tree);
}

template SpectralDecompositionTreecd sqrt(SpectralDecompositionTreecd X);
template SpectralDecompositionTreed sqrt(SpectralDecompositionTreed X);
template MatrixTreecd sqrt(MatrixTreecd X, const Tree& tree);
template MatrixTreed sqrt(MatrixTreed X, const Tree& tree);

template SpectralDecompositionTreecd inverse(SpectralDecompositionTreecd X, double eps);
template SpectralDecompositionTreed inverse(SpectralDecompositionTreed X, double eps);
template MatrixTreecd inverse(MatrixTreecd X, const Tree& tree, double eps);
template MatrixTreed inverse(MatrixTreed X, const Tree& tree, double eps);

template MatrixTreecd to_matrixtree(const SpectralDecompositionTreecd& X, const Tree& tree);
template MatrixTreed to_matrixtree(const SpectralDecompositionTreed& X, const Tree& tree);

template void CanonicalTransformation(TensorTreecd& Psi, const Tree& tree, bool orthogonal);
template void CanonicalTransformation(TensorTreed& Psi, const Tree& tree, bool orthogonal);

template
class SpectralDecompositionTree<complex<double>>;

template
class SpectralDecompositionTree<double>;

