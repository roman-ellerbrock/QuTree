//
// Created by Roman Ellerbrock on 2/16/21.
//

#include "TreeClasses/SparseTensorTree.h"

template<typename T>
void SparseTensorTree<T>::initialize(const Tree& tree) {
	attributes_.clear();
	for (const Node *const node_ptr : sparseTree()) {
		attributes_.emplace_back(Tensor<T>(node_ptr->shape()));
	}
}

template class SparseTensorTree<complex<double>>;

template class SparseTensorTree<double>;

