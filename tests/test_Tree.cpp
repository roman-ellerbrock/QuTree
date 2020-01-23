//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TensorTreeBasis.h"

SUITE (Tree) {

	TEST(TensorTreeBasisGenerator) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);

		for (size_t n_modes = 2; n_modes < 18; ++n_modes) {
			TensorTreeBasis basis(n_modes, n_leaf, n_node);

			for (const Node& node : basis) {
				const TensorDim& tdim = node.TDim();
				if (node.IsToplayer()) {
						CHECK_EQUAL(tdim_top, tdim);
				} else if (node.IsBottomlayer()) {
						CHECK_EQUAL(tdim_bottom, tdim);
				} else {
						CHECK_EQUAL(tdim_upper, tdim);
				}
			}

//				CHECK_EQUAL(13, basis.nNodes());
		}
	}

	TEST(TensorTreeBasisFileIO) {
		size_t n_leaf = 4;
		size_t n_node = 2;
		TensorDim tdim_top({n_node, n_node}, 1);
		TensorDim tdim_upper({n_node, n_node}, n_node);
		TensorDim tdim_bottom({n_leaf}, n_node);
		size_t n_modes = 13;

		TensorTreeBasis basis(n_modes, n_leaf, n_node);
		{
			ofstream os("TTBasis.IO.tmp.dat");
			basis.Write(os);
			os.close();
		}

		TTBasis basis2("TTBasis.IO.tmp.dat");
		CHECK_EQUAL(basis2.nNodes(), basis.nNodes());

	}

}

