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
			string file("1	-2"
						"	 	1	-2"
						"	 		1	-1"
						" 				1	0	0"
						"	 		1	-1"
						"	 			1	0	1"
						"	 	1	-2"
						"	 		1	-1"
						"	 			1	0	2"
						"	 		1	-1"
						"	 			1	0	3");
			stringstream is(file);
			root_ = readNode(is);

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

	class tree : root {
	public:
		tree() {
			tree_.setRoot(root_);
		}

		Tree tree_;
	};

	TEST_FIXTURE(tree, copy){
		Tree tree(tree_);
		CHECK_EQUAL(true, tree == tree_);
	}

	TEST_FIXTURE(tree, assign){
		Tree tree = tree_;
			CHECK_EQUAL(true, tree == tree_);
	}

	TEST_FIXTURE(tree, mconstruct){
		Tree tree2(tree_);
		Tree tree(move(tree2));
			CHECK_EQUAL(true, tree == tree_);
	}

	TEST_FIXTURE(tree, massign){
		Tree tree2 = tree_;
		Tree tree =move(tree2);
			CHECK_EQUAL(true, tree == tree_);
	}

	TEST_FIXTURE(tree, read) {
		/// prepare tree_
		for (size_t l = 0; l < tree_.nLeaves(); ++l) {
			Leaf& leaf = tree_.leafArray()[l];
			leaf.par().par0_ = 1.;
			leaf.par().par1_ = 0.;
			leaf.par().par2_ = 0.;
			leaf.par().par3_ = 1.;
		}

		string file("1	-2\n"
					"	 	1	-2\n"
					"	 		1	-1\n"
					" 				1	0	0\n"
					"	 		1	-1\n"
					"	 			1	0	1\n"
					"	 	1	-2\n"
					"	 		1	-1\n"
					"	 			1	0	2\n"
					"	 		1	-1\n"
					"	 			1	0	3\n"
					"1.	0.	0.	1.\n"
					"1.	0.	0.	1.\n"
					"1.	0.	0.	1.\n"
					"1.	0.	0.	1.\n");
		stringstream is(file);
		Tree tree(is);
			CHECK_EQUAL(tree_, tree);
	}
}