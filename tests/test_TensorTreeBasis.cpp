//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"

SUITE (TensorTreeBasis) {
	TEST (TensorTreeBasis_Generator) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorShape tdim_top({n_node, n_node, 1});
		TensorShape tdim_upper({n_node, n_node, n_node});
		TensorShape tdim_bottom({n_leaf, n_node});

		for (size_t n_modes = 2; n_modes < 18; ++n_modes) {
			Tree tree = TreeFactory::balancedTree(n_modes, n_leaf, n_node);

			for (const Node& node : tree) {
				const TensorShape& tdim = node.shape();
				if (node.isToplayer()) {
						CHECK_EQUAL(tdim_top, tdim);
				} else if (node.isBottomlayer()) {
						CHECK_EQUAL(tdim_bottom, tdim);
				} else {
						CHECK_EQUAL(tdim_upper, tdim);
				}
			}
		}
	}

	TEST (TensorTreeBasis_FileIO) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorShape tdim_top({n_node, n_node, 1});
		TensorShape tdim_upper({n_node, n_node, n_node});
		TensorShape tdim_bottom({n_leaf, n_node});
		size_t n_modes = 13;

		Tree tree = TreeFactory::balancedTree(n_modes, n_leaf, n_node);
		{
			ofstream os("TTBasis.IO.tmp.dat");
			tree.write(os);
			os.close();
		}
		Tree tree2("TTBasis.IO.tmp.dat");
			CHECK_EQUAL(tree2.nNodes(), tree.nNodes());
	}

	TEST (TensorTreeBasis_Reindexing) {
		size_t n_modes = 9;
		Tree tree = TreeFactory::balancedTree(n_modes, 2, 4);

		map<size_t, size_t> Map;
		for (size_t k = 0; k < n_modes; ++k) {
			Map[k] = n_modes - 1 - k;
		}
		tree.reindexLeafModes(Map);

		size_t k = 0;
		for (const Node& node : tree) {
			if (node.isBottomlayer()) {
				const Leaf& leaf = node.getLeaf();
					CHECK_EQUAL(k++, leaf.mode());
			}
		}
	}

	TEST (TensorTreeBasis_Train) {
		size_t nLeaves = 12;
		auto tree = TreeFactory::unbalancedTree(nLeaves, 4, 2, 6);
			CHECK_EQUAL(2 * nLeaves - 1, tree.nNodes());
	}

	TEST (TensorTreeBasis_Copy) {
		/// Construct a tree and check that it works
		Tree tree = TreeFactory::balancedTree(12, 4, 3);
			CHECK_EQUAL(true, tree.isWorking());

		{
			/// Copy-constructor test
			Tree tree_copy_construct(tree);
				CHECK_EQUAL(true, tree_copy_construct.isWorking());

			/// Move constructor
			Tree tree_move_construct(move(tree_copy_construct));
				CHECK_EQUAL(true, tree_move_construct.isWorking());
		}

		{
			/// Copy-asignment test
			Tree tree_copy_asign = tree;
				CHECK_EQUAL(true, tree_copy_asign.isWorking());

			/// Move asignment operator
			Tree tree_move_asign = move(tree_copy_asign);
				CHECK_EQUAL(true, tree_move_asign.isWorking());
		}
	}
}