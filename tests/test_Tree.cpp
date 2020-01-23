//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TensorTreeBasis.h"

SUITE (Tree) {

	TEST(TensorTreeBasisGenerator) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		size_t n_modes = 7;
		TensorDim tdim_top({n_node, n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);

		TensorTreeBasis basis(n_modes, n_leaf, n_node);

		for (const Node& node : basis) {
			const TensorDim& tdim = node.TDim();
			if (node.IsToplayer()) {
				CHECK_EQUAL(tdim, tdim_top);
			} else if (node.IsBottomlayer()) {
				CHECK_EQUAL(tdim, tdim_bottom);
			} else {
				CHECK_EQUAL(tdim, tdim_upper);
			}
		}
	}

}

