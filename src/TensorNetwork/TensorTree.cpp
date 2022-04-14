//
// Created by Roman Ellerbrock on 12/3/21.
//

#include "TensorNetwork/TensorTree.h"
#include "Tensor/TensorLapack.h"

template<typename T>
TensorTree<T>::TensorTree(const Tree& tree,
	function<Tensor<T>(const TensorShape&)> gen)
	:SubTreeAttribute<Tensor<T>>(tree) {

	for (const Node* node : this->nodes_) {
		(*this)[*node] = gen(node->shape_);
	}

	normalize();
}

template<typename T>
void TensorTree<T>::normalize(double eps) {

	for (const Edge* edge: this->edges_) {
		const Tensor<T>& phi = (*this)[edge->from()];
		size_t k = edge->outIdx();
		(*this)[edge] = ::normalize(phi, k, eps);
	}
}

template<typename T>
void TensorTree<T>::print() const {
	cout << "Nodes:\n";
	for (const Node* node : SubTree::nodes_) {
		node->info();
		(*this)[node].print();
	}
	cout << "Edges:\n";
	for (const Edge* edge : SubTree::edges_) {
		edge->info();
		(*this)[edge].print();
	}
}

template
class TensorTree<double>;

template
class TensorTree<complex<double>>;

template<typename T>
TensorTree<T>& operator+=(TensorTree<T>& A, const TensorTree<T>& add) {
}

template<typename T>
TensorTree<T>& operator*=(TensorTree<T>& A, T factor);

