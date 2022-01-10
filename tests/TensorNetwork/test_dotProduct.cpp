//
// Created by Roman Ellerbrock on 12/14/21.
//

#include "UnitTest++/UnitTest++.h"
#include "TensorNetwork/DotProduct.h"
#include "Tree/TreeFactory.h"


SUITE (dotProduct) {
	class Trees {
	public:
		Trees() {
			tree_ = balancedTree(4, 3, 2);
			psi_ = TensorTreecd(tree_, deltacd);
		}

		Tree tree_;
		TensorTreecd psi_;
	};

	double eps = 1e-10;

	TEST_FIXTURE (Trees, dotproduct) {
		TensorTreecd S(tree_);
		dotProduct(S, psi_, psi_);
		for (const Node& node : tree_.nodes()) {
				CHECK_CLOSE(0., residual(S[node], Tensorcd()), eps);
		}
		for (const Edge& edge : tree_.edges()) {
				CHECK_CLOSE(0., residual(S[edge], identitycd({2, 2})), eps);
		}
	}

	TEST_FIXTURE(Trees, fullContraction) {
		TensorTreecd S(tree_);
		dotProduct(S, psi_, psi_);
		TensorTreecd tr = fullContraction(psi_, S, psi_);
		Tensorcd h({1});
		h(0) = 1;
		for (const Node& node : psi_) {
			CHECK_CLOSE(0., residual(h, tr[node]), eps);
		}
	}

}