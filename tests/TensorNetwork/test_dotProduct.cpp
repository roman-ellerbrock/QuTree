//
// Created by Roman Ellerbrock on 12/14/21.
//

#include "UnitTest++/UnitTest++.h"
#include "TensorNetwork/contractions.h"
#include "Tree/TreeFactory.h"
#include "Util/GateOperators.h"


SUITE (contraction) {
	class Trees {
	public:
		Trees() {
			tree_ = balancedTree(4, 2, 2);
			psi_ = TensorTreecd(tree_, deltacd);
		}

		Tree tree_;
		TensorTreecd psi_;
	};

	double eps = 1e-10;

	TEST_FIXTURE (Trees, dotproduct) {
		TensorTreecd S = matrixTreecd(tree_, {});
		contraction(S, psi_, psi_);
		for (const Edge* edge : S.edges_) {
				CHECK_CLOSE(0., residual(S[edge], identitycd({2, 2})), eps);
		}
	}

	TEST_FIXTURE(Trees, CNotMatrices) {
		ProductOperatorcd P(sigma_x(), 0);
		TensorTreecd S = matrixTreecd(tree_, P.targetLeaves());
		contraction(S, psi_, psi_, P);
	}

}
