//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"

template class TensorTree<double>;
template class TensorTree<complex<double>>;

template <typename T>
TensorTree<T> sizedTensorTree(const Tree& tree) {
	TensorTree<T> psi(tree);
	for (const Node& node : tree.nodes()) {
		psi[node] = Tensor<T>(node.shape_);
	}

	for (const Edge& edge: tree.upEdges()) {
		psi[edge] = Tensor<T>(edge.shape());
	}

	for (const Edge& edge: tree.downEdges()) {
		psi[edge] = Tensor<T>(edge.shape());
	}
	return psi;
}

template TensorTree<complex<double>> sizedTensorTree<complex<double>>(const Tree& tree);
template TensorTree<double> sizedTensorTree<double>(const Tree& tree);


template <typename T>
TensorTree<T> occupiedTensorTree(const Tree& tree) {
	TensorTree<T> psi(tree);
	for (const Node& node : tree.nodes()) {
		psi[node] = random<T>(node.shape_);
	}

	for (const Edge& edge: tree.upEdges()) {
		psi[edge] = random<T>(edge.shape());
	}

	for (const Edge& edge: tree.downEdges()) {
		psi[edge] = random<T>(edge.shape());
	}
	return psi;
}

template TensorTree<complex<double>> occupiedTensorTree<complex<double>>(const Tree& tree);
template TensorTree<double> occupiedTensorTree<double>(const Tree& tree);

