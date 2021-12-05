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
		TensorTreed Psi = sizedTensorTreed(tree_);
		TensorShape s0({2, 2, 1});
		TensorShape s1({2, 2, 2});
		TensorShape s2({3, 2});
		vector<TensorShape> shape({s0, s1, s2, s2, s1, s2, s2});
		size_t i = 0;
		for (const Node& node : tree_) {
				CHECK_EQUAL(shape[i++], Psi[node].shape_);
		}

		TensorShape eshape({2, 2});
		for (const Edge& edge : tree_.downEdges()) {
				CHECK_EQUAL(eshape, Psi[edge].shape_);
		}

		for (const Edge& edge : tree_.upEdges()) {
				CHECK_EQUAL(eshape, Psi[edge].shape_);
		}
	}

	TEST_FIXTURE(Trees, a) {

	}

}

