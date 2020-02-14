//
// Created by Roman Ellerbrock on 2/3/20.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeHandling/Tree.h"
#include "TreeHandling/TTBasisFactory.h"

SUITE (TensorTreeBasis) {
	TEST (TensorTreeBasis_Generator) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);

		for (size_t n_modes = 2; n_modes < 18; ++n_modes) {
			Tree tree(n_modes, n_leaf, n_node);

			for (const Node& node : tree) {
				const TensorDim& tdim = node.TDim();
				if (node.IsToplayer()) {
						CHECK_EQUAL(tdim_top, tdim);
				} else if (node.IsBottomlayer()) {
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
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);
		size_t n_modes = 13;

		Tree tree(n_modes, n_leaf, n_node);
		{
			ofstream os("TTBasis.IO.tmp.dat");
			tree.Write(os);
			os.close();
		}
		Tree tree2("TTBasis.IO.tmp.dat");
			CHECK_EQUAL(tree2.nNodes(), tree.nNodes());
	}

	TEST (TensorTreeBasis_Reindexing) {
		size_t n_modes = 9;
		Tree tree(n_modes, 2, 4);

		map<size_t, size_t> Map;
		for (size_t k = 0; k < n_modes; ++k) {
			Map[k] = n_modes - 1 - k;
		}
		tree.ReindexLeafModes(Map);

		size_t k = 0;
		for (const Node& node : tree) {
			if (node.IsBottomlayer()) {
				const Leaf& leaf = node.PhysCoord();
					CHECK_EQUAL(k++, leaf.Mode());
			}
		}
	}

	TEST (TensorTreeBasis_Train) {
		size_t nLeaves = 12;
		auto tree = TTBasisFactory::TensorTrain(nLeaves, 4, 2, 6);
			CHECK_EQUAL(2 * nLeaves - 1, tree.nNodes());
	}

	TEST (TensorTreeBasis_Copy) {
		/// Construct a tree and check that it works
		Tree tree(12, 4, 3);
			CHECK_EQUAL(true, tree.IsWorking());

		{
			/// Copy-constructor test
			Tree tree_copy_construct(tree);
				CHECK_EQUAL(true, tree_copy_construct.IsWorking());

			/// Move constructor
			Tree tree_move_construct(move(tree_copy_construct));
				CHECK_EQUAL(true, tree_move_construct.IsWorking());
		}

		{
			/// Copy-asignment test
			Tree tree_copy_asign = tree;
				CHECK_EQUAL(true, tree_copy_asign.IsWorking());

			/// Move asignment operator
			Tree tree_move_asign = move(tree_copy_asign);
				CHECK_EQUAL(true, tree_move_asign.IsWorking());
		}
	}
}