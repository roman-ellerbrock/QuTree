//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators//MatrixListTree.h"


void MatrixListTree::Initialize(const Tree& tree) {
	attributes_.resize(tree.nNodes());
}

void MatrixListTree::print(const Tree& tree, ostream& os) const {

	for (const Node& node : tree) {
		node.info();
		const MatrixList& hs = operator[](node);
		for (const Matrixcd& h : hs) {
			h.print();
		}
	}
}

