//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"
#include "Tensor/TensorLapack.h"

template<typename T>
TensorTree<T>::TensorTree(const Tree& tree, function<Tensor<T>(const TensorShape&)> gen)
	:TensorTree(tree) {
	for (const Node& node : tree.nodes()) {
		(*this)[node] = gen(node.shape_);
	}
	normalize();
}

template<typename T>
void TensorTree<T>::normalize(double eps) {

	for (const Edge& edge: this->edges()) {
		const Tensor<T>& phi = (*this)[edge.from()];
		size_t k = edge.outIdx();
		(*this)[edge] = ::normalize(phi, k, eps);
	}
}

template<typename T>
void TensorTree<T>::print() const {
	cout << "Nodes:\n";
	for (const Node& node : Tree::nodes()) {
		node.info();
		(*this)[node].print();
	}
	cout << "up-Edges:\n";
	for (const Edge& edge : Tree::upEdges()) {
		(*this)[edge].print();
	}
	cout << "down-Edges:\n";
	for (const Edge& edge : Tree::downEdges()) {
		(*this)[edge].print();
	}
}

template
class TensorTree<double>;

template
class TensorTree<complex<double>>;

