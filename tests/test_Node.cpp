//
// Created by Roman Ellerbrock on 11/19/21.
//
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tree/Leaf.h"


SUITE(Nodes) {

	TEST(Leaf) {
		BasisParameters par({1., 0., 0., 1.,
							 10, 2, 0});
		Leaf leaf;
		leaf.initialize(par);
			CHECK_EQUAL(true, leaf.parent_ == nullptr);
			CHECK_EQUAL(0, leaf.api_.basis()->par_.mode_);
			CHECK_EQUAL(2, leaf.api_.basis()->par_.type_);
	}



/*	TEST(Node_construct_default) {
		Node node;
			CHECK_EQUAL(0, node.nChildren());
			CHECK_EQUAL(1, node.isRoot());
	}*/

}
