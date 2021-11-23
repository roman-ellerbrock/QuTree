//
// Created by Roman Ellerbrock on 11/19/21.
//
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tree/Leaf.h"
#include "Tree/Node.h"

SUITE (Nodes) {

	class lf {
	public:
		lf() {
			par_ = BasisParameters({1., 0., 0., 1.,
									10, 2, 0});
			leaf_.initialize(par_);
		}

		BasisParameters par_;
		Leaf leaf_;
	};

	// ==================================================
	// ==== Leaf ========================================
	// ==================================================

	TEST_FIXTURE (lf, Leaf) {
			CHECK_EQUAL(true, leaf_.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == leaf_.api_.basis()->par_);
	}

	TEST_FIXTURE(lf, leaf_copy) {
		Leaf l2 (leaf_);
			CHECK_EQUAL(true, l2.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == l2.api_.basis()->par_);
	}

	TEST_FIXTURE(lf, leaf_equal) {
		Leaf l2 = leaf_;
			CHECK_EQUAL(true, l2.parent_ == nullptr);
			CHECK_EQUAL(true, par_ == l2.api_.basis()->par_);
	}


	// ==================================================
	// ==== Node ========================================
	// ==================================================

	TEST(Node) {
		Node node;
	}


}
