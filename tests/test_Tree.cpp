//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "Tree.h"

SUITE (Tree) {
/*	TEST(TreeNodeTests) {
		TreeNode leaf(2);
		TreeNode bottom(2);
		bottom.push_back(leaf);

		TreeNode upper = Double(bottom);
		for (size_t l = 0; l < 3; ++l) {
			upper = Double(upper);
		}
		upper.MakeRoot();
		upper.GenInput(cout);

		// write out stuff.
		auto next = upper.NextNode();
		size_t i = 0;
		while (next != &upper) {
			print(next->GetPath());
			next = upper.NextNode();
			i++;
			if (i == 100) {
				getchar();
			}
		}
		print(upper.GetPath());
	}
 */

	TEST (TreeClass) {
		TensorDimTree tree(3, 2, 4, 2);
		for (const TreeNode node : tree) {
			node.TDim().print();
		}

		for (size_t mode = 0; mode < tree.nLeaves(); ++mode) {
			const auto& node = tree.Leaf(mode);
			cout << "mode: " << node.content_.f << endl;
		}
	}
}

