//
// Created by Roman Ellerbrock on 11/26/20.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/SparseTree.h"
#include "TreeShape/TreeFactory.h"

/**
 * 	for
 *	-> Sparse tree (tail, noinvert)
 * 	-> Sparse tree (notail, noinvert)
 * 	-> Sparse tree (tail, invert)
 * 	-> Sparse tree (notail, invert)
 * 	do
 * 		test: sparse-iterator
 * 		test: dense-iterator + active-test
 *	enddo
 **/


SUITE (SparseTree) {

	void checkNumberActiveNodes(const SparseTree& stree, const Tree& tree, size_t n_expect) {
		/// Count number of nodes in sparse tree. Check whether number of
		/// nodes marked as active matches number of nodes in the
		size_t num1 = 0;
		for (const Node* node_ptr : stree) {
			num1++;
		}
			CHECK_EQUAL(n_expect, num1);
		size_t num2 = 0;
		for (const Node& node : tree) {
			if (stree.Active(node)) {
				num2++;
			}
		}
			CHECK_EQUAL(n_expect, num2);
	}

	TEST (TestTail) {
		auto tree = TreeFactory::BalancedTree(8, 2, 2);
		{
			vector<size_t> leaves({0, 2});
			SparseTree stree1(leaves, tree, true, false);
			checkNumberActiveNodes(stree1, tree, 6);
			SparseTree stree2(leaves, tree, false, false);
			checkNumberActiveNodes(stree2, tree, 5);
		}
		{
			vector<size_t> leaf({3});
			SparseTree stree1(leaf, tree, true, false);
			checkNumberActiveNodes(stree1, tree, 4);
			SparseTree stree2(leaf, tree, false, false);
			checkNumberActiveNodes(stree2, tree, 1);
		}
	}
}

