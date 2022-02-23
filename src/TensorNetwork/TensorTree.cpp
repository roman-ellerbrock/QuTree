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
TensorTree<T> matrixTree(const Tree& tree, const vector<size_t>& idx) {
	SubTree subTree(tree, idx);
	TensorTree<T> Psi(subTree);
	for (const Edge* edge : Psi.edges_) {
		Psi[edge] = Tensor<T>(edge->shape());
	}
	return Psi;
}

typedef complex<double> cd;
typedef double d;

template TensorTree<cd> matrixTree<cd>(const Tree& tree, const vector<size_t>&);
template TensorTree<d> matrixTree<d>(const Tree& tree, const vector<size_t>&);
