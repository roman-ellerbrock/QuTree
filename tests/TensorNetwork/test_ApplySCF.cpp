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
			psi_ = TensorTreecd(tree_, deltacd);

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
		SOPcd cnot = CNot(0, 2);
		auto mat = matrixTreecd(tree_, cnot);

		ProductOperatorcd P(hadamard(), 1);
		chi_ = P.apply(chi_);

		TensorTreecd hchi(tree_, randomcd);

		contraction(mat, hchi, chi_, cnot);

		apply(hchi, mat, chi_, cnot, 10);

		hchi.print();

	}
}
