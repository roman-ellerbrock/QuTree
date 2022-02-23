//
// Created by Roman Ellerbrock on 11/30/21.
//

#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tree/EdgeArray.h"

SUITE (Edges) {


	class tree4 {
	public:
		tree4() {
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
			nodes_ = NodeArray(root_);

			child_ = root_.child_[0];
			up_ = Edge(child_, root_);
			down_ = Edge(root_, child_);
		}

		Node root_;
		Node child_;
		NodeArray nodes_;
		Edge down_;
		Edge up_;
	};

	TEST_FIXTURE (tree4, from) {
			CHECK_EQUAL(child_, up_.from());
			CHECK_EQUAL(root_, down_.from());
	}

	TEST_FIXTURE (tree4, to) {
			CHECK_EQUAL(root_, up_.to());
			CHECK_EQUAL(child_, down_.to());
	}

	TEST_FIXTURE (tree4, up) {
			CHECK_EQUAL(root_, up_.up());
			CHECK_EQUAL(root_, down_.up());
	}

	TEST_FIXTURE (tree4, down) {
			CHECK_EQUAL(child_, up_.down());
			CHECK_EQUAL(child_, down_.down());
	}

	TEST_FIXTURE (tree4, inIdx) {
			CHECK_EQUAL(0, up_.inIdx());
			CHECK_EQUAL(2, down_.inIdx());
	}

	TEST_FIXTURE (tree4, outIdx) {
			CHECK_EQUAL(2, up_.outIdx());
			CHECK_EQUAL(0, down_.outIdx());
	}

	TEST_FIXTURE (tree4, address) {
		EdgeArray edges(nodes_);
		int i = 0;
		vector<int> add(edges.size(), -1);
		for (const Edge& edge : edges.down_) {
				CHECK_EQUAL(i, edge.address());
				CHECK_EQUAL(-1, add[edge.address()]);
				i += 2;
			add[edge.address()] = 1;
		}
		i--;
		for (const Edge& edge : edges.up_) {
				CHECK_EQUAL(i, edge.address());
				CHECK_EQUAL(-1, add[edge.address()]);
				i -= 2;
			add[edge.address()] = 2;
		}
	}

	TEST_FIXTURE (tree4, incomingEdge) {
		Edge expected(root_, root_.child_[1]);
		Edge edge = outgoingEdge(root_, 1);
			CHECK_EQUAL(expected, edge);
	}

	TEST_FIXTURE (tree4, outgoingEdge) {
		Edge expected(root_.child_[1], root_);
		Edge edge = incomingEdge(root_, 1);
			CHECK_EQUAL(expected, edge);
	}

	TEST_FIXTURE (tree4, incomingEdges) {
		vector<Edge> expected = {
			Edge(child_.child_[0], child_),
			Edge(child_.child_[1], child_),
			Edge(root_, child_),
		};
			CHECK_EQUAL(expected, incomingEdges(child_));
	}

	TEST_FIXTURE (tree4, outgoingEdges) {
		vector<Edge> expected = {
			Edge(child_, child_.child_[0]),
			Edge(child_, child_.child_[1]),
			Edge(child_, root_),
		};
			CHECK_EQUAL(expected, outgoingEdges(child_));
	}

	TEST_FIXTURE (tree4, preUpEdge) {
		vector<Edge> expected = {
			Edge(child_.child_[0], child_),
			Edge(child_.child_[1], child_)
		};
		Edge edge(child_, root_);
			CHECK_EQUAL(expected, preUpEdges(edge));
	}

	TEST_FIXTURE (tree4, preDownEdge) {
		vector<Edge> expected = {
			Edge(child_.child_[0], child_),
			Edge(root_, child_)
		};
		Edge edge(child_, child_.child_[1]);
			CHECK_EQUAL(expected, preDownEdges(edge));
	}

	TEST_FIXTURE (tree4, preEdges) {
		vector<Edge> expected = {
			Edge(child_.child_[0], child_),
			Edge(child_.child_[1], child_)
		};
		Edge edge(child_, root_);
			CHECK_EQUAL(expected, preEdges(edge));
	}

	TEST_FIXTURE (tree4, preEdges2) {
		vector<Edge> expected = {
			Edge(child_.child_[0], child_),
			Edge(root_, child_)
		};
		Edge edge(child_, child_.child_[1]);
			CHECK_EQUAL(expected, preEdges(edge));
	}
}

