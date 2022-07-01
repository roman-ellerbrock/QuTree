//
// Created by Roman Ellerbrock on 4/2/22.
//
#include "TensorNetwork/TensorTreeFactory.h"


template<typename T>
TensorTree<T> matrixTree(const Tree& tree, const SubTreeParameters& par) {
	SubTree subTree(tree, par);
	TensorTree<T> Psi(subTree);
	for (const Node* node : Psi.nodes_) {
		Psi[node] = Tensor<T>({1});
	}

	for (const Edge* edge : Psi.edges_) {
		Psi[edge] = Tensor<T>(edge->shape());
	}
	return Psi;
}

typedef complex<double> cd;
typedef double d;

template TensorTree<cd> matrixTree<cd>(const Tree& tree, const SubTreeParameters&);
template TensorTree<d> matrixTree<d>(const Tree& tree, const SubTreeParameters&);


template <typename T>
TensorTree<T> matrixTree(const Tree& tree, const ProductOperator<T>& P) {
	return matrixTree<T>(tree, P.targetLeaves());
}

template TensorTree<cd> matrixTree(const Tree& tree, const ProductOperator<cd>& P);
template TensorTree<d> matrixTree(const Tree& tree, const ProductOperator<d>& P);


template <typename T>
vector<TensorTree<T>> matrixTree(const Tree& tree, const SumOfProductsOperator<T>& H) {
	vector<TensorTree<T>> Svec;
	for (const ProductOperator<T>& P : H) {
		Svec.emplace_back(matrixTree<T>(tree, P));
	}
	return Svec;
}

template vector<TensorTree<cd>> matrixTree(const Tree& tree, const SumOfProductsOperator<cd>& H);
template vector<TensorTree<d>> matrixTree(const Tree& tree, const SumOfProductsOperator<d>& H);


template <typename T>
TensorTree<T> tensorTree(TensorTree<T> Psi, const TensorTree<T>& mat,
	function<Tensor<T>(const TensorShape&)> f) {

	for (const Node* node : mat.nodes_) {
		Psi[node] = f(node->shape_);
	}

	for (const Edge* edge : mat.edges_) {
		const Node& node = edge->from();
		Psi[edge] = Psi[node];
	}
	Psi.normalize();

	return Psi;
}

template TensorTree<cd> tensorTree(TensorTree<cd> Psi, const TensorTree<cd>& mat,
	function<Tensor<cd>(const TensorShape&)> f);
template TensorTree<d> tensorTree(TensorTree<d> Psi, const TensorTree<d>& mat,
	function<Tensor<d>(const TensorShape&)> f);
