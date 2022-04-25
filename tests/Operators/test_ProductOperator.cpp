//
// Created by Roman Ellerbrock on 12/23/21.
//
#include "UnitTest++/UnitTest++.h"
#include "Tensor/Tensor"
#include "Tree/TreeFactory.h"
#include "Operator/ProductOperator.h"

SUITE (ProductOperator) {

	class Operators {
	public:
		Operators()
			: tree_(balancedTree(3, 2, 2)),
			  A_(identitycd({2, 2})),
			  hA_({2, 2}),
			  basis_(tree_.leafArray()[0].basis_),
			  x_({2, 2}) {
			x_(1, 0) = 1.;
			x_(0, 1) = 1.;
			Psi_ = TensorTreecd(tree_, deltacd);
		}

		~Operators() = default;

		Tree tree_;
		BasisAPI& basis_;
		Tensorcd A_;
		Tensorcd hA_;
		Matrixcd x_;
		TensorTreecd Psi_;
	};

	double eps = 1e-10;

	TEST_FIXTURE (Operators, ProductOperator) {
		POcd P;
		P.push_back(x_, 1);
		P.applyReference(Psi_);
		const Leaf* leaf = Psi_.leaves_[1];
		const Tensorcd& hA_ = Psi_[leaf->parent()];
			CHECK_CLOSE(0., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(1., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

	TEST_FIXTURE (Operators, doubleOperator) {
		POcd P;
		P.push_back(x_, 1);
		P.push_back(x_, 1);
		P.applyReference(Psi_);
		const Leaf* leaf = Psi_.leaves_[1];
		const Tensorcd& hA_ = Psi_[leaf->parent()];
			CHECK_CLOSE(1., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

	TEST_FIXTURE (Operators, reindexed) {
		for (size_t l = 0; l < tree_.leafArray().size(); ++l){
			Leaf& leaf = tree_.leafArray()[l];
			leaf.par().mode_ = 2 - l;
		}
		tree_.update();

		POcd P;
		P.push_back(x_, 0);
		P.applyReference(Psi_);
		const Leaf* leaf = Psi_.leaves_[0];
		const Tensorcd& hA_ = Psi_[leaf->parent()];
			CHECK_CLOSE(0., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(1., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

	TEST_FIXTURE (Operators, operatorTimes) {
		POcd P;
		P.push_back(x_, 1);
		POcd Q = P * P;

		Q.applyReference(Psi_);
		const Leaf* leaf = Psi_.leaves_[0];
		const Tensorcd& hA_ = Psi_[leaf->parent()];
			CHECK_CLOSE(1., abs(hA_(0, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 0)), eps);
			CHECK_CLOSE(0., abs(hA_(0, 1)), eps);
			CHECK_CLOSE(0., abs(hA_(1, 1)), eps);
	}

}

