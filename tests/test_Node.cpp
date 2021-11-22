//
// Created by Roman Ellerbrock on 11/19/21.
//
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tree/Leaf.h"
#include "Tree/Node.h"


SUITE(Nodes) {

	TEST(Leaf) {
		Leaf leaf;
	}

/*	TEST(Node_construct_default) {
		Node node;
			CHECK_EQUAL(0, node.nChildren());
			CHECK_EQUAL(1, node.isRoot());
	}*/

}
