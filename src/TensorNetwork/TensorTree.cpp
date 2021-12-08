//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"

template<typename T>
TensorTree<T>::TensorTree(const Tree& tree, function<Tensor<T>(const TensorShape&)> gen)
	:TensorTree(tree)
{
	for (const Node& node : tree.nodes()) {
		(*this)[node] = gen(node.shape_);
	}

}

template<typename T>
void TensorTree<T>::normalize() {

	for (const Edge& edge: this->edges()) {
//		(*this)[edge] = gen(edge.from().shape_);
	}

}


template class TensorTree<double>;
template class TensorTree<complex<double>>;

template <typename T>
TensorTree<T> sizedTensorTree(const Tree& tree) {
	TensorTree<T> psi(tree);
	for (const Node& node : tree.nodes()) {
		psi[node] = Tensor<T>(node.shape_);
	}

	for (const Edge& edge: tree.edges()) {
		psi[edge] = Tensor<T>(edge.from().shape_);
	}

	return psi;
}



template TensorTree<complex<double>> sizedTensorTree<complex<double>>(const Tree& tree);
template TensorTree<double> sizedTensorTree<double>(const Tree& tree);


template <typename T>
TensorTree<T> randomTensorTree(const Tree& tree) {
	TensorTree<T> psi(tree);
	for (const Node& node : tree.nodes()) {
		psi[node] = random<T>(node.shape_);
	}

	for (const Edge& edge: tree.upEdges()) {
		psi[edge] = random<T>(edge.from().shape_);
	}

	return psi;
}

template TensorTree<complex<double>> randomTensorTree<complex<double>>(const Tree& tree);
template TensorTree<double> randomTensorTree<double>(const Tree& tree);

