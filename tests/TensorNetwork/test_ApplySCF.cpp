//
// Created by Roman Ellerbrock on 4/14/22.
//

#include "UnitTest++/UnitTest++.h"
#include "TensorNetwork/ApplySCF.h"
#include "TensorNetwork/contractions.h"
#include "Tree/TreeFactory.h"
#include "Util/GateOperators.h"


SUITE (ApplySCF) {
	class Trees {
	public:
		Trees() {
			tree_ = balancedTree(4, 2, 2);
			psi_ = TensorTreecd(tree_, identitycd);

			chi_ = TensorTreecd(tree_, deltacd);
			/// Manually set basis
			for (const Edge *edge : chi_.edges_) {
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
	};

	double eps = 1e-10;

	TEST_FIXTURE (Trees, applyCNot) {
		SOPcd cnot = CNot(1, 3);
		auto Svec = matrixTreecd(tree_, cnot);
		contraction(Svec, chi_, chi_, cnot);


	}
}
