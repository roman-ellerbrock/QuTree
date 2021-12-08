//
// Created by Roman Ellerbrock on 12/3/21.
//
#include "TensorNetwork/TensorTree.h"
#include "UnitTest++/UnitTest++.h"
#include "Tree/TreeFactory.h"

SUITE(TensorTree) {

	class Trees {
	public:
		Trees() {
			tree_ = balancedTree(4, 3, 2);
		}

		Tree tree_;
	};

	TEST_FIXTURE(Trees, construct) {
		TensorTreecd Psi(tree_);
			CHECK_EQUAL(7, Psi.nodes_.size());
			CHECK_EQUAL(6, Psi.upEdges_.size());
			CHECK_EQUAL(6, Psi.downEdges_.size());
	}

	TEST_FIXTURE(Trees, sized) {
		TensorTreecd Psi(tree_, deltacd);
		for (const Node& node : tree_) {
				CHECK_EQUAL(node.shape_, Psi[node].shape_);
		}

		for (const Edge& edge : tree_.edges()) {
//				CHECK_EQUAL(edge.from().shape_, Psi[edge].shape_);
		}
	}

	TEST_FIXTURE(Trees, symmetricTTN) {

	}

}

