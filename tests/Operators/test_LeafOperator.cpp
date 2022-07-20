//
// Created by Roman Ellerbrock on 12/21/21.
//

#include "UnitTest++/UnitTest++.h"
#include "Operator/LeafFunction.h"
#include "Operator/LeafMatrix.h"
#include "Tensor/Tensor"
#include "Tree/TreeFactory.h"

SUITE(LeafOperators) {
	class Operators {
	public:
		Operators()
			: tree_(balancedTree(3, 3, 2)),
		A_(identitycd({2, 2})),
		hA_({2, 2}),
		basis_(tree_.leafArray()[0].basis_) {

		}
		~Operators() = default;

		Tree tree_;
		BasisAPI& basis_;
		Tensorcd A_;
		Tensorcd hA_;
	};

	double eps = 1e-10;

	void fun1(const BasisAPI& basis, Tensorcd& hA, const Tensorcd& A) {
		Tensorcd h{2, 2};
		h(0, 1) = 1;
		h(1, 0) = 1;

		matrixTensor(hA, h, A, 0);
	}

	TEST_FIXTURE(Operators, LeafFunction) {
		function<void(const BasisAPI& basis, Tensorcd&, const Tensorcd& A)> fo(fun1);
		LeafFunctioncd G(fun1);
		G.apply(basis_, hA_, A_);
			CHECK_CLOSE(0., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(1., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(1., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

	TEST_FIXTURE(Operators, LeafMatrix) {
		Matrixcd h({2, 2});
		h(1, 0) = 1.;
		LeafMatrixcd G(h);
		G.apply(basis_, hA_, A_);
			CHECK_CLOSE(0., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(1., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

}