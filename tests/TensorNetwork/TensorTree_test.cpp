//
// Created by Roman Ellerbrock on 12/3/21.
//
#include <gtest/gtest.h>
#include "TensorNetwork/contractions.h"
#include "TensorNetwork/TensorTree.h"
#include "Tree/TreeFactory.h"

class Trees : public ::testing::Test
{
public:
	Trees()
	{
		tree_ = balancedTree(4, 2, 2);
	}

	Tree tree_;
	double eps_{1e-10};
};


TEST_F(Trees, constructorSize)
{
	TensorTreecd Psi(tree_, randomcd);
	EXPECT_EQ(7, Psi.nodes_.size());
	EXPECT_EQ(12, Psi.edges_.size());
}

TEST_F(Trees, TensorShapes)
{
	TensorTreecd Psi(tree_, deltacd);
	for (const Node &node : tree_)
	{
		EXPECT_EQ(node.shape_, Psi[node].shape_);
	}

	for (const Edge &edge : tree_.edges())
	{
		EXPECT_EQ(edge.from().shape_, Psi[edge].shape_);
	}
}

TEST_F(Trees, constructorTensor)
{
	TensorTreecd Psi(tree_, deltacd);
	for (const Node &node : tree_)
	{
		EXPECT_NEAR(0., residual(deltacd(node.shape_), Psi[node]), eps_);
	}

	for (const Edge &edge : tree_.edges())
	{
		const Tensorcd &phi = Psi[edge];
		auto delta = contraction(phi, phi, edge);
		EXPECT_NEAR(0., residual(delta, identitycd(delta.shape_)), eps_);
	}
}

TEST_F(Trees, consistency)
{
	TensorTreecd Psi(tree_, randomcd);
	for (const Node *node : Psi.nodes_)
	{
		Psi[node] = randomcd(node->shape_);
	}
	for (const Edge *edge : Psi.edges_)
	{
		Psi[edge] = Psi[edge->from()];
	}

	TensorTreecd Chi(tree_, randomcd);
	for (const Node *node : Psi.nodes_)
	{
		Chi[node] = Psi[node];
	}
	Chi.normalize();
}

TEST_F(Trees, matrixTree)
{
	auto mat = matrixTreecd(tree_);
	for (const Node &node : tree_)
	{
		EXPECT_EQ(true, mat.nodes_.contains(node.address_));
		/// @TODO: define how node tensor should look like. Until now there it isn't defined
	}

	for (const Edge &edge : tree_.edges())
	{
		EXPECT_EQ(true, mat.edges_.contains(edge.address()));
		EXPECT_EQ(edge.shape(), mat[edge].shape_);
	}
}

TEST_F(Trees, matrixTreeEdges)
{
	auto mat = matrixTreecd(tree_, SubTreeParameters({}, false));
	for (const Node &node : tree_)
	{
		EXPECT_EQ(false, mat.nodes_.contains(node.address_));
	}

	for (const Edge &edge : tree_.edges())
	{
		EXPECT_EQ(true, mat.edges_.contains(edge.address()));
	}
}

TEST_F(Trees, matrixTreeNodes)
{
	auto mat = matrixTreecd(tree_, SubTreeParameters({}, true, false));
	for (const Node &node : tree_)
	{
		EXPECT_EQ(true, mat.nodes_.contains(node.address_));
	}

	for (const Edge &edge : tree_.edges())
	{
		EXPECT_EQ(false, mat.edges_.contains(edge.address()));
	}
}

TEST_F(Trees, matrixTreePreEdges)
{
	auto mat = matrixTreecd(tree_);
	for (const Edge *edge : mat.edges_)
	{
		EXPECT_EQ(edge->shape(), mat[edge].shape_);
	}

	for (const Edge *e : mat.edges_)
	{
		vector<Edge> pre = mat.preEdges(e);
		for (const Edge pe : pre)
		{
			EXPECT_EQ(true, mat.edges_.contains(pe.address()));
			EXPECT_EQ(e->shape(), mat[e].shape_);
		}
	}
}

TEST_F(Trees, matrixTreePreEdgesSparse)
{
	vector<size_t> idx = {1, 3};
	auto mat = matrixTreecd(tree_, idx);
	for (const Edge *e : mat.edges_)
	{
		auto pre = mat.preEdges(e);
		for (const Edge &edge : pre)
		{
			EXPECT_EQ(true, mat.edges_.contains(edge.address()));
			EXPECT_EQ(edge.shape(), mat[edge].shape_);
		}
	}
}
