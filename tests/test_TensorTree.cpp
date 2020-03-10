//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/TensorTree_Implementation.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"

SUITE (TensorTree) {

	TEST (TensorTree_FILE_IO) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 4);
		TensorTreecd T(tree);
		T.Write("TT.tmp.dat");
		TensorTreecd Q(tree);
		ifstream is("TT.tmp.dat");
		is >> Q;
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : tree) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_RandomGenerate) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 2);
		TensorTreecd T(tree);
		mt19937 gen(2468);
		T.FillRandom(gen, tree, false);
		string filename("TT.RNG.dat");
		T.Write(filename);
		TensorTreecd Q(filename);
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : tree) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_Train) {
		size_t nLeaves = 12;
		auto tree = TreeFactory::UnbalancedTree(nLeaves, 4, 2, 6);
			CHECK_EQUAL(2 * nLeaves - 1, tree.nNodes());
		mt19937 gen(2468);
		TensorTreecd Psi(gen, tree);
			CHECK_EQUAL(tree.nNodes(), Psi.size());
	}



	TEST (TreeTest) {
		Tree tree = TreeFactory::BalancedTree(12, 5, 2);
		TensorTreecd Psi(tree);
		for (const Tensorcd& A : Psi) {
			TensorShape shape = A.shape();
		}
		for (const Node& node : tree) {
			Tensorcd& A = Psi[node];
		}

	}

}

