//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTree.h"

TensorOperatorTree::TensorOperatorTree(const MLOd& M,
	const Tree& tree)
	: TensorOperatorTree(tree) {
	occupy(tree);
	vector<size_t> idxs(tree.nLeaves(), 0);

	for (size_t k = 0; k < M.size(); ++k) {
		size_t mode = M.mode(k);
		const Leaf& leaf = tree.getLeaf(mode);
		const auto& node = (const Node&) leaf.parent();

		const shared_ptr<LeafOperator<double>>& h = M[k];
		Matrixd hmat = toMatrix(*h, leaf);

		setLeafOperator(hmat, idxs[mode], node);
		idxs[mode]++;
	}
}

TensorOperatorTree::TensorOperatorTree(const SOPd& S,
	const Tree& tree)
	: TensorOperatorTree(tree) {
	occupy(tree);
	assert(S.size() > 0);
}

TensorOperatorTree::TensorOperatorTree(const Tree& tree) {
	attributes_.clear();
	for (const Node& node : tree) {
		const TensorShape& shape = node.shape();
		if (node.isBottomlayer()) {
			vector<size_t> dim({shape.lastBefore(), shape.lastBefore(), shape.lastDimension()});
			TensorShape oshape(dim);
			attributes_.emplace_back(Tensord(oshape));
		} else {
			attributes_.emplace_back(Tensord(shape));
		}
	}
	occupy(tree);
}

void TensorOperatorTree::occupy(const Tree& tree) {
	for (const Node& node : tree) {
		Tensord& Phi = operator[](node);
		Phi.zero();
		for (size_t i = 0; i < Phi.shape().lastDimension(); ++i) {
			Phi(i, i) = 1.;
		}
		if (node.isBottomlayer()) {
			if (!(Phi.shape().order() == 3)) {
				cerr << "Wrong order in tensor operator tree.\n";
				exit(1);
			}
			size_t dim = sqrt(Phi.shape().lastBefore());
			setLeafOperator(identityMatrixd(dim), 0, node);
		}
	}
}

void TensorOperatorTree::print(const Tree& tree) const {
	for (const Node& node : tree) {
		node.info();
		attributes_[node.address()].print();
	}
}

void TensorOperatorTree::setLeafOperator(const Matrixd& m,
	size_t operator_idx, const Node& node) {

	const TensorShape& shape = node.shape();
	assert(m.dim1() * m.dim2() == shape.lastBefore());
	assert(operator_idx < shape.lastDimension());

	Tensord& h = operator[](node);
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		h(i, operator_idx) = m[i];
	}
}

