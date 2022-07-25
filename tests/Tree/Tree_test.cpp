//
// Created by Roman Ellerbrock on 11/19/21.
//

#include <gtest/gtest.h>
#include "Tree/LeafArray.h"
#include "Tree/Tree.h"

class root : public ::testing::Test
{
public:
	root()
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

		Node *p = &root_;
		size_t addr = 0;
		while (p)
		{
			p->address_ = addr++;
			p = sweep(p);
		}
	}

	Node root_;
};

TEST_F(root, init)
{
	Tree tree;
	tree.setRoot(root_);
	Node *p = &root_;
	for (const Node &node : tree)
	{
		EXPECT_EQ(true, *p == node);
		if (p)
		{
			p = sweep(p);
		}
	}
}

class tree : public root
{
public:
	tree()
	{
		tree_.setRoot(root_);
	}

	Tree tree_;
};

TEST_F(tree, copy)
{
	Tree tree(tree_);
	EXPECT_EQ(true, tree == tree_);
}

TEST_F(tree, assign)
{
	Tree tree = tree_;
	EXPECT_EQ(true, tree == tree_);
}

TEST_F(tree, mconstruct)
{
	Tree tree2(tree_);
	Tree tree(move(tree2));
	EXPECT_EQ(true, tree == tree_);
}

TEST_F(tree, massign)
{
	Tree tree2 = tree_;
	Tree tree = move(tree2);
	EXPECT_EQ(true, tree == tree_);
}

TEST_F(tree, read)
{
	/// prepare tree_
	for (size_t l = 0; l < tree_.nLeaves(); ++l)
	{
		Leaf &leaf = tree_.leafArray()[l];
		leaf.par().par0_ = 1.;
		leaf.par().par1_ = 0.;
		leaf.par().par2_ = 0.;
		leaf.par().par3_ = 1.;
	}

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
	stringstream is(file);
	Tree tree(is);
	EXPECT_EQ(tree_, tree);
}

/**
	Tree tree
	for (const Node& node : tree.nodes()) {}
	for (const Node& node : tree.nodes.reverse()) {}
	for (const Edge& edge : tree.edges()) {}
	for (const Edge& edge : tree.edges().reverse()) {}

	Node node;
	for (const Edge& edge : tree.incoming(node)) {}
	for (const Edge& edge : tree.outgoing(node)) {}

	Edge e;
	for (const Edge& edge : tree.prior(e)) {}


	/// Usage examples:
	Psi[node] = Tensorcd();
	S[edge] = Tensorcd();
**/
