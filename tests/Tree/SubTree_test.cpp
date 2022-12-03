//
// Created by Roman Ellerbrock on 1/21/22.
//

#include <gtest/gtest.h>
#include "Tree/LeafArray.h"
#include "Tree/Tree.h"
#include "Tree/SubTree.h"

class SomeTree : public ::testing::Test
{
public:
	SomeTree()
	{
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

TEST_F(SomeTree, init)
{
	SubTree subtree(tree_);

	size_t i = 0;
	for (const Node *node : subtree.nodes_)
	{
		const Node &other = tree_.nodes()[i++];
		EXPECT_EQ(&other, node);
	}
}

TEST_F(SomeTree, initPartial)
{
	vector<size_t> idx = {0, 1};
	SubTree stree(tree_, idx);

	const Node &node0 = tree_.nodeArray()[0];
	const Node &node1 = tree_.nodeArray()[1];
	const Node &node2 = tree_.nodeArray()[2];
	const Node &node3 = tree_.nodeArray()[3];
	EXPECT_EQ(node0.address_, stree.nodes_[0]->address_);
	EXPECT_EQ(node1.address_, stree.nodes_[1]->address_);
	EXPECT_EQ(node2.address_, stree.nodes_[2]->address_);
	EXPECT_EQ(node3.address_, stree.nodes_[3]->address_);
}

TEST_F(SomeTree, initPartial2)
{
	vector<size_t> idx = {1, 3};
	SubTree stree(tree_, idx);

	const Node &node0 = tree_.nodeArray()[0];
	const Node &node1 = tree_.nodeArray()[1];
	const Node &node3 = tree_.nodeArray()[3];
	const Node &node4 = tree_.nodeArray()[4];
	const Node &node6 = tree_.nodeArray()[6];
	EXPECT_EQ(true, stree.nodes_.contains(node0.address_));
	EXPECT_EQ(true, stree.nodes_.contains(node1.address_));
	EXPECT_EQ(true, stree.nodes_.contains(node3.address_));
	EXPECT_EQ(true, stree.nodes_.contains(node4.address_));
	EXPECT_EQ(true, stree.nodes_.contains(node6.address_));

	const Node &node2 = tree_.nodeArray()[2];
	const Node &node5 = tree_.nodeArray()[5];
	EXPECT_EQ(false, stree.nodes_.contains(node2.address_));
	EXPECT_EQ(false, stree.nodes_.contains(node5.address_));
}

TEST_F(SomeTree, testEdges)
{
	vector<size_t> idx = {1, 3};
	SubTree stree(tree_, idx);

	vector<size_t> edgeidx = {0, 2, 3, 5, 6, 8, 9, 11};
	EXPECT_EQ(edgeidx.size(), stree.edges_.size());
	for (auto i : edgeidx)
	{
		const Edge &edge = tree_.edges()[i];
		EXPECT_EQ(true, stree.edges_.contains(edge.address()));
	}
}

TEST_F(SomeTree, testPreEdges)
{
	vector<size_t> idx = {1, 3};
	SubTree stree(tree_, idx);

	for (const Edge *edge : stree.edges_)
	{
		EXPECT_EQ(true, stree.edges_.contains(edge->address()));
		const vector<Edge> preEdges = stree.preEdges(edge);
		for (const Edge &e : preEdges) {
			EXPECT_EQ(true, stree.edges_.contains(e.address()));
		}
	}
}

TEST_F(SomeTree, testInEdges)
{
	vector<size_t> idx = {1, 3};
	SubTree stree(tree_, idx);
	for (const Node *node : stree.nodes_)
	{
		const vector<Edge> inEdges = stree.incomingEdges(node);
		for (const Edge &e : inEdges)
		{
			EXPECT_EQ(true, stree.edges_.contains(e.address()));
		}
	}
}
