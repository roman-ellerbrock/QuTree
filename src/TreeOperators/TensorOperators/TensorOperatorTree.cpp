//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTree.h"

TensorOperatorTree::TensorOperatorTree(const Tree &tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		if (node.isBottomlayer()) {
			if (node.shape().order() != 3) {
				cerr << "Error: Cannot initialize TensorOperatorTree from tree:\n";
				cerr << "Order at bottomlayer does not equal 3." << endl;
			}
		}
		attributes_.emplace_back(Tensorcd(node.shape()));
	}
}

void TensorOperatorTree::Occupy(const Tree& tree) {
	for (const Node& node : tree) {
		Tensorcd& Phi = operator[](node);
		Phi.Zero();
		for (size_t i = 0; i < Phi.shape().lastDimension(); ++i) {
			Phi(i, i) = 1.;
		}
	}
}

void TensorOperatorTree::print(const Tree& tree) const {
	for (const Node& node : tree) {
		node.info();
		attributes_[node.Address()].print();
	}
}

void TensorOperatorTree::setLeafOperator(const LeafMatrixcd& M,
	size_t operator_idx, const Node& node) {

	const TensorShape& shape = node.shape();
	const Matrixcd& m = M.matrix();
	assert(m.Dim1() * m.Dim2() == shape.lastBefore());
	assert(operator_idx < shape.lastDimension());

	Tensorcd& h = operator[](node);
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		h(i, operator_idx) = m[i];
	}
}



