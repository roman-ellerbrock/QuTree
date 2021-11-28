//
// Created by Roman Ellerbrock on 11/19/21.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/LeafArray.h"
#include "Tree/Tree.h"

SUITE (Tree) {

	class root {
	public:
		root() {
			BasisParameters par = BasisParameters({1., 0., 0., 1.,
												   10, 2, 0});
			Node inter;
			Node bottom;
			Leaf leaf;
			leaf.initialize(par);
			bottom.push_back(leaf);
			inter.push_back(bottom);
			inter.push_back(bottom);

			root_.push_back(inter);
			root_.push_back(inter);

			Node *p = &root_;
			size_t addr = 0;
			while (p) {
				p->address_ = addr++;
				p = sweep(p);
			}
		}

		Node root_;
	};

	TEST_FIXTURE (root, init) {
		Tree tree;
		tree.setRoot(root_);
		Node *p = &root_;
		for (const Node& node : tree) {
			CHECK_EQUAL(true, *p == node);
			if (p) { p = sweep(p); }
		}
	}

	TEST_FIXTURE(root, ){

	}
}