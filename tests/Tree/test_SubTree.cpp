//
// Created by Roman Ellerbrock on 1/21/22.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/LeafArray.h"
#include "Tree/Tree.h"
#include "Tree/SubTree.h"

SUITE (SubTree) {

	class SomeTree {
	public:
		SomeTree() {
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
			stringstream ss(file);
			tree_ = Tree(ss);
		}

		~SomeTree() = default;

		Tree tree_;
	};

	TEST_FIXTURE (SomeTree, init) {
		SubTree subtree(tree_);

		size_t i = 0;
		for (const Node* node : subtree.nodes_) {
			const Node& other = tree_.nodes()[i++];
			CHECK_EQUAL(&other, node);
		}
	}

	TEST_FIXTURE (SomeTree, initPartial) {
		vector<size_t> idx = {0, 1};
		SubTree stree(tree_, idx);

		const Node& node0 = tree_.nodeArray()[0];
		const Node& node1 = tree_.nodeArray()[1];
		const Node& node2 = tree_.nodeArray()[2];
		const Node& node3 = tree_.nodeArray()[3];
			CHECK_EQUAL(node0.address_, stree.nodes_[0]->address_);
			CHECK_EQUAL(node1.address_, stree.nodes_[1]->address_);
			CHECK_EQUAL(node2.address_, stree.nodes_[2]->address_);
			CHECK_EQUAL(node3.address_, stree.nodes_[3]->address_);
	}

	TEST_FIXTURE (SomeTree, initPartial2) {
		vector<size_t> idx = {1, 3};
		SubTree stree(tree_, idx);

		const Node& node0 = tree_.nodeArray()[0];
		const Node& node1 = tree_.nodeArray()[1];
		const Node& node3 = tree_.nodeArray()[3];
		const Node& node4 = tree_.nodeArray()[4];
		const Node& node6 = tree_.nodeArray()[6];
			CHECK_EQUAL(true, stree.nodes_.contains(node0.address_));
			CHECK_EQUAL(true, stree.nodes_.contains(node1.address_));
			CHECK_EQUAL(true, stree.nodes_.contains(node3.address_));
			CHECK_EQUAL(true, stree.nodes_.contains(node4.address_));
			CHECK_EQUAL(true, stree.nodes_.contains(node6.address_));

		const Node& node2 = tree_.nodeArray()[2];
		const Node& node5 = tree_.nodeArray()[5];
			CHECK_EQUAL(false, stree.nodes_.contains(node2.address_));
			CHECK_EQUAL(false, stree.nodes_.contains(node5.address_));
	}

	TEST_FIXTURE (SomeTree, testEdges) {
		vector<size_t> idx = {1, 3};
		SubTree stree(tree_, idx);

		vector<size_t> edgeidx = {0, 2, 3, 5, 6, 8, 9, 11};
		for (auto i : edgeidx) {
			const Edge& edge = tree_.edges()[i];
			CHECK_EQUAL(true, stree.edges_.contains(edge.address()));
		}
	}

	TEST_FIXTURE(SomeTree, testPreEdges) {
		vector<size_t> idx = {1, 3};
		SubTree stree(tree_, idx);

		for (const Edge* edge : stree.edges_) {
			CHECK_EQUAL(true, stree.edges_.contains(edge->address()));
			const vector<Edge> preEdges = stree.preEdges(edge);
			for (const Edge& e : preEdges) {
				CHECK_EQUAL(true, stree.edges_.contains(e.address()));
			}
		}

	}

}

