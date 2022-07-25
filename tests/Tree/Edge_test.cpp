//
// Created by Roman Ellerbrock on 11/30/21.
//

#include <gtest/gtest.h>
#include <iostream>
#include "Tree/EdgeArray.h"

class tree4 : public ::testing::Test
{
public:
	tree4()
	{
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

TEST_F(tree4, from)
{
	EXPECT_EQ(child_, up_.from());
	EXPECT_EQ(root_, down_.from());
}

TEST_F(tree4, to)
{
	EXPECT_EQ(root_, up_.to());
	EXPECT_EQ(child_, down_.to());
}

TEST_F(tree4, up)
{
	EXPECT_EQ(root_, up_.up());
	EXPECT_EQ(root_, down_.up());
}

TEST_F(tree4, down)
{
	EXPECT_EQ(child_, up_.down());
	EXPECT_EQ(child_, down_.down());
}

TEST_F(tree4, inIdx)
{
	EXPECT_EQ(0, up_.inIdx());
	EXPECT_EQ(2, down_.inIdx());
}

TEST_F(tree4, outIdx)
{
	EXPECT_EQ(2, up_.outIdx());
	EXPECT_EQ(0, down_.outIdx());
}

TEST_F(tree4, address)
{
	EdgeArray edges(nodes_);
	int i = 0;
	vector<int> add(edges.size(), -1);
	for (const Edge &edge : edges.down_)
	{
		EXPECT_EQ(i, edge.address());
		EXPECT_EQ(-1, add[edge.address()]);
		i += 2;
		add[edge.address()] = 1;
	}
	i--;
	for (const Edge &edge : edges.up_)
	{
		EXPECT_EQ(i, edge.address());
		EXPECT_EQ(-1, add[edge.address()]);
		i -= 2;
		add[edge.address()] = 2;
	}
}

TEST_F(tree4, incomingEdge)
{
	Edge expected(root_, root_.child_[1]);
	Edge edge = outgoingEdge(root_, 1);
	EXPECT_EQ(expected, edge);
}

TEST_F(tree4, outgoingEdge)
{
	Edge expected(root_.child_[1], root_);
	Edge edge = incomingEdge(root_, 1);
	EXPECT_EQ(expected, edge);
}

TEST_F(tree4, incomingEdges)
{
	vector<Edge> expected = {
		Edge(child_.child_[0], child_),
		Edge(child_.child_[1], child_),
		Edge(root_, child_),
	};
	EXPECT_EQ(expected, incomingEdges(child_));
}

TEST_F(tree4, outgoingEdges)
{
	vector<Edge> expected = {
		Edge(child_, child_.child_[0]),
		Edge(child_, child_.child_[1]),
		Edge(child_, root_),
	};
	EXPECT_EQ(expected, outgoingEdges(child_));
}

TEST_F(tree4, preUpEdge)
{
	vector<Edge> expected = {
		Edge(child_.child_[0], child_),
		Edge(child_.child_[1], child_)};
	Edge edge(child_, root_);
	EXPECT_EQ(expected, preUpEdges(edge));
}

TEST_F(tree4, preDownEdge)
{
	vector<Edge> expected = {
		Edge(child_.child_[0], child_),
		Edge(root_, child_)};
	Edge edge(child_, child_.child_[1]);
	EXPECT_EQ(expected, preDownEdges(edge));
}

TEST_F(tree4, preEdges)
{
	vector<Edge> expected = {
		Edge(child_.child_[0], child_),
		Edge(child_.child_[1], child_)};
	Edge edge(child_, root_);
	EXPECT_EQ(expected, preEdges(edge));
}

TEST_F(tree4, preEdges2)
{
	vector<Edge> expected = {
		Edge(child_.child_[0], child_),
		Edge(root_, child_)};
	Edge edge(child_, child_.child_[1]);
	EXPECT_EQ(expected, preEdges(edge));
}

TEST_F(tree4, postEdges)
{
	vector<Edge> expected = {
		Edge(child_, child_.child_[1]),
		Edge(child_, root_)};
	Edge edge(child_.child_[0], child_);
	EXPECT_EQ(expected, postEdges(edge));
}

TEST_F(tree4, postEdges2)
{
	vector<Edge> expected = {
		Edge(root_.child_[1], root_.child_[1].child_[0]),
		Edge(root_.child_[1], root_.child_[1].child_[1]),
	};
	Edge edge(root_, root_.child_[1]);
	EXPECT_EQ(expected, postEdges(edge));
}
