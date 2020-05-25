//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTree.h"

TensorOperatorTree::TensorOperatorTree(const MLOcd& M,
	const Tree& tree) : TensorOperatorTree(tree) {
	Occupy(tree);
	vector<size_t> idxs(tree.nLeaves());

	for (size_t k = 0; k < M.size(); ++k) {
		size_t mode = M.Mode(k);
		const Leaf& leaf = tree.GetLeaf(mode);
		const auto& node = (const Node&) leaf.Up();

		const shared_ptr<LeafOperatorcd>& h = M[k];
		Matrixcd hmat = toMatrix(*h, leaf);
		hmat.print();
		setLeafOperator(hmat, idxs[mode], node);
		idxs[mode]++;
	}

}

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
		if (node.isBottomlayer()) {
			size_t dim = sqrt(node.shape().lastBefore());
			setLeafOperator(IdentityMatrixcd(dim), 0, node);
		}
	}
}

void TensorOperatorTree::print(const Tree& tree) const {
	for (const Node& node : tree) {
		node.info();
		attributes_[node.Address()].print();
	}
}

void TensorOperatorTree::setLeafOperator(const Matrixcd& m,
	size_t operator_idx, const Node& node) {

	const TensorShape& shape = node.shape();
	assert(m.Dim1() * m.Dim2() == shape.lastBefore());
	assert(operator_idx < shape.lastDimension());

	Tensorcd& h = operator[](node);
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		h(i, operator_idx) = m[i];
	}
}

