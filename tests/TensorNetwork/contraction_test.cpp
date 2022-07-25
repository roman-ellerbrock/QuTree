//
// Created by Roman Ellerbrock on 12/14/21.
//

#include <gtest/gtest.h>
#include "TensorNetwork/contractions.h"
#include "Tree/TreeFactory.h"
#include "Util/GateOperators.h"

class Trees : public ::testing::Test
{
public:
	Trees()
	{
		tree_ = balancedTree(4, 2, 2);
		psi_ = TensorTreecd(tree_, identitycd);

		chi_ = TensorTreecd(tree_, deltacd);
		/// Manually set basis
		for (const Edge *edge : chi_.edges_)
		{
			auto shape = edge->from().shape_;
			size_t idx = edge->outIdx();
			shape = transpose(shape, idx);
			auto I = identitycd(shape);
			I = transpose(I, idx, true);
			chi_[edge] = I;
		}
	}

	Tree tree_;
	TensorTreecd psi_;
	TensorTreecd chi_;
	double eps_{1e-10};
};

TEST_F(Trees, dotproduct)
{
	TensorTreecd S = matrixTreecd(tree_);
	contraction(S, psi_, psi_);
	for (const Edge *edge : S.edges_)
	{
		EXPECT_NEAR(0., residual(S[edge], identitycd({2, 2})), eps_);
	}
}

TEST_F(Trees, residual)
{
	EXPECT_NEAR(0., residual(psi_, psi_, tree_), eps_);
}

TEST_F(Trees, consistency)
{
	Tree tree_ = balancedTree(2, 2, 2);
	TensorTreecd Psi(tree_, randomcd);
	for (const Node *node : Psi.nodes_)
	{
		Psi[node] = randomcd(node->shape_);
	}
	for (const Edge *edge : Psi.edges_)
	{
		Psi[edge] = Psi[edge->from()];
	}

	TensorTreecd Chi(Psi);
	Chi.normalize();

	TensorTreecd S = matrixTreecd(tree_);
	contraction(S, Psi, Chi);
	TensorTreecd psipsi = matrixTreecd(tree_);
	contraction(psipsi, Psi, Psi);

	ProductOperatorcd P;
	for (const Node *node : Chi.nodes_)
	{
		auto sChi = Chi[node];
		apply(sChi, S, P, node);
		auto sPsi = Psi[node];
		apply(sPsi, psipsi, P, node);
		EXPECT_NEAR(0., residual(sChi, sPsi), eps_);
	}
}

TEST_F(Trees, consistency2)
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

	TensorTreecd Chi(Psi);
	Chi.normalize();

	TensorTreecd S = matrixTreecd(tree_);
	contraction(S, Psi, Chi);
	TensorTreecd psipsi = matrixTreecd(tree_);
	contraction(psipsi, Psi, Psi);

	ProductOperatorcd P;
	for (const Node *node : Chi.nodes_)
	{
		auto sChi = Chi[node];
		apply(sChi, S, P, node);
		auto sPsi = Psi[node];
		apply(sPsi, psipsi, P, node);
		EXPECT_NEAR(0., residual(sChi, sPsi), eps_);
	}
}

TEST_F(Trees, ProductOperator)
{
	ProductOperatorcd P(sigma_x(), 0);
	TensorTreecd S = matrixTreecd(tree_, P);
	contraction(S, psi_, psi_, P);
	for (const Edge *edge : S.edges_)
	{
		if (edge->isUpEdge())
		{
			EXPECT_NEAR(0., residual(S[edge], sigma_x()), eps_);
		}
		else
		{
			EXPECT_NEAR(0., residual(S[edge], identitycd({2, 2})), eps_);
		}
	}
}

TEST_F(Trees, ProductOperator2)
{
	ProductOperatorcd P(sigma_x(), 0);
	P.push_back(identitycd({2, 2}), 2);
	TensorTreecd S = matrixTreecd(tree_, P);
	contraction(S, psi_, psi_, P);
	auto PPsi = P.apply(psi_);
	TensorTreecd S2 = matrixTreecd(tree_, P);
	contraction(S2, psi_, PPsi);
	/// Check consistency between two methods
	for (const Edge *edge : S.edges_)
	{
		EXPECT_NEAR(0., residual(S[edge], S2[edge]), eps_);
	}
}

TEST_F(Trees, SumOfProductOperators)
{
	/// Check whether matrixTree(tree_, SOP) provides correct matrix-Trees.
	SOPcd cnot = CNot(1, 3);
	auto Svec = matrixTreecd(tree_, cnot);

	auto S0 = matrixTreecd(tree_, cnot[0]);
	auto S1 = matrixTreecd(tree_, cnot[1]);
	vector<TensorTreecd> SvecRef = {S0, S1};

	EXPECT_EQ(2, Svec.size());
	for (size_t l = 0; l < Svec.size(); ++l)
	{
		const auto &s = Svec[l];
		const auto &sRef = SvecRef[l];
		EXPECT_EQ(sRef.nodes_.size(), s.nodes_.size());
		EXPECT_EQ(sRef.edges_.size(), s.edges_.size());
	}
}

TEST_F(Trees, SOPcontraction)
{
	/// Calculate matrixTree for SOP and compare to manually calculated version
	SOPcd cnot = CNot(1, 3);
	auto Svec = matrixTreecd(tree_, cnot);
	contraction(Svec, psi_, psi_, cnot);

	auto S0 = matrixTreecd(tree_, cnot[0]);
	contraction(S0, psi_, psi_, cnot[0]);
	for (const Node *node : S0.nodes_)
	{
		EXPECT_NEAR(0., residual(S0[node], Svec[0][node]), eps_);
	}

	auto S1 = matrixTreecd(tree_, cnot[1]);
	contraction(S1, psi_, psi_, cnot[1]);
	for (const Node *node : S1.nodes_)
	{
		EXPECT_NEAR(0., residual(S1[node], Svec[1][node]), eps_);
	}
}

TEST_F(Trees, applyPO)
{
	ProductOperatorcd P(sigma_x(), 0);
	P.push_back(sigma_x(), 2);

	auto hmat = matrixTreecd(tree_, P);
	contraction(hmat, psi_, psi_, P);
	TensorTreecd hpsi(psi_);
	apply(hpsi, hmat, P);
}

TEST_F(Trees, applySOPedge)
{
	SOPcd cnot = CNot(1, 3);
	auto Svec = matrixTreecd(tree_, cnot);
	contraction(Svec, chi_, chi_, cnot);

	TensorTreecd chi0(chi_);
	for (const Edge *edge : Svec[0].edges_)
	{
		apply(chi_[edge], Svec, cnot, edge);
	}

	EXPECT_NEAR(0., residual(chi0, chi_, tree_), eps_);

	/// @TODO: flip control bit or apply Hadamard and see what happens
	/// Apply sigma_x to 1st node
	/*		const Leaf& leaf = tree_.leafArray()[1];
			const Node& node = leaf.parent();
			auto edges = outgoingEdges(node);
			const Edge& edge = edges.front();
			chi0[edge] = matrixTensor(sigma_x(), chi0[edge], node.parentIdx());
			chi_ = chi0;
			contraction(Svec, chi_, chi_, cnot);
			spsi = chi_;
			for (const Edge *edge : Svec[0].edges_) {
				apply(chi_, spsi, Svec, cnot, edge);
			}
			chi_.print();*/
}
