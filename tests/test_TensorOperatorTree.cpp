//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeOperators/TensorOperators/MatrixListTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/TensorTreeFunctions.h"
#include "TreeOperators/TensorOperators/TensorOperatorTreeFunctions.h"
#include <UnitTest++/UnitTest++.h>

SUITE (TensorOperatorTree) {

	double eps = 1e-7;

	class HelperFactory {
	public:
		HelperFactory() {
			Initialize();
		}

		~HelperFactory() = default;
		mt19937 rng_;
		Tree tree_;
		Tree optree_;
		TensorTreecd Psi_;
		TensorOperatorTree H_;
		LeafMatrixcd leafI_;
		LeafMatrixcd leafX_;

		void Initialize() {
			rng_ = mt19937(1990);
			tree_ = TreeFactory::BalancedTree(12, 2, 2);
			optree_ = TreeFactory::OperatorTree(tree_);
			Psi_ = TensorTreecd(rng_, tree_);

			Matrixcd I = identityMatrixcd(2);
			Matrixcd X(2, 2);
			X(0, 1) = 1.;
			X(1, 0) = 1.;
			leafI_ = LeafMatrixcd(I);
			leafX_ = LeafMatrixcd(X);

			H_ = TensorOperatorTree(optree_);
			H_.occupy(optree_);
			for (const Node& node : optree_) {
				if (node.isBottomlayer()) {
					H_.setLeafOperator(leafI_, 0, node);
					H_.setLeafOperator(leafX_, 1, node);
				}
			}
			H_.print(tree_);
		}
	};

	TEST (OperatorTree) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 2);
		Tree optree = TreeFactory::OperatorTree(tree);
		for (const Node& node : optree) {
			if (!node.isToplayer()) {
					CHECK_EQUAL(3, node.shape().order());
					CHECK_EQUAL(2, node.shape().lastDimension());
					CHECK_EQUAL(4, node.shape().lastBefore());
			}
		}
	}

/*	TEST (TensorOperatorTree) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 2);
		Tree optree = TreeFactory::OperatorTree(tree);
		TensorOperatorTree H(tree);
		H.Occupy(tree);
		for (const Node node : tree) {
			if (node.isBottomlayer()) {
				Tensorcd& h = H[node];
				CHECK_CLOSE(1., abs(h(0, 0)), 1e-7);
				CHECK_CLOSE(0., abs(h(1, 0)), 1e-7);
				CHECK_CLOSE(0., abs(h(2, 0)), 1e-7);
				CHECK_CLOSE(1., abs(h(3, 0)), 1e-7);
			}
		}
		H.print(tree);
	}*/

/*
	TEST (SetLeafOperator) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 2);
		Tree optree = TreeFactory::OperatorTree(tree);

		Matrixcd I = IdentityMatrixcd(2);
		Matrixcd X(2, 2);
		LeafMatrixcd leafI(I);
		LeafMatrixcd leafX(X);

		TensorOperatorTree H(optree);
		for (const Node& node : optree) {
			if (node.isBottomlayer()) {
				H.setLeafOperator(leafI, 0, node);
				H.setLeafOperator(leafX, 1, node);
			}
		}
	}

	TEST (MatrixListTree) {
		Tree tree = TreeFactory::BalancedTree(12, 2, 2);
		mt19937 gen(1990);
		TensorTreecd Psi(gen, tree);

		Tree optree = TreeFactory::OperatorTree(tree);

		Matrixcd I = IdentityMatrixcd(2);
		Matrixcd X(2, 2);
		X(0, 1) = 1.;
		X(1, 0) = 1.;
		LeafMatrixcd leafI(I);
		LeafMatrixcd leafX(X);

		TensorOperatorTree H(optree);
		for (const Node& node : optree) {
			if (node.isBottomlayer()) {
				H.setLeafOperator(leafI, 0, node);
				H.setLeafOperator(leafX, 1, node);
			}
		}

		MatrixListTree Hrep = TreeFunctions::Represent(Psi, H, tree);
	}

	TEST_FIXTURE (HelperFactory, CNot) {
		auto Hrep = TreeFunctions::Represent(Psi_, H_, tree_);
		for (const Node& node : tree_) {
			Matrixcd I = IdentityMatrixcd(node.shape().lastDimension());
			const MatrixList& hs = Hrep[node];
				CHECK_CLOSE(0., Residual(I, hs[0]), 1e-6);
		}

		MatrixListTree Hmean = TreeFunctions::Contraction(Psi_, H_, Hrep, tree_);

		for (const Node& node : tree_) {
			if (!node.isToplayer()) {
				const MatrixList& hmeans = Hmean[node];
				size_t dim = hmeans[0].Dim1();
					CHECK_CLOSE(1., abs(hmeans[0](0, 0)), eps);
					CHECK_CLOSE(0., abs(hmeans[0](1, 0)), eps);
					CHECK_CLOSE(0., abs(hmeans[0](0, 1)), eps);
					CHECK_CLOSE(0., abs(hmeans[0](1, 1)), eps);
					CHECK_CLOSE(0., Residual(hmeans[1], Matrixcd(dim, dim)), eps);
			}
		}
	}

	TEST_FIXTURE (HelperFactory, Compress) {
		MatrixTreecd rho = TreeFunctions::Contraction(H_, optree_, true);
		SpectralDecompositionTreecd X(rho, optree_);
		CanonicalTransformation(H_, optree_, true);
		TreeFunctions::Adjust(H_, optree_, X, 1e-7);
		for (const Node& node : optree_) {
			const TensorShape& shape = node.shape();
				CHECK_EQUAL(3, shape.order());
			if (node.isBottomlayer()) {
					CHECK_EQUAL(4, shape.totalDimension());
				Tensorcd Phi(shape);
				Phi(0) = 1.;
				Phi(3) = 1.;
					CHECK_CLOSE(0., Residual(H_[node], Phi), eps);
			} else {
					CHECK_EQUAL(1, shape.totalDimension());
				Tensorcd Phi(shape);
				Phi(0) = 1.;
					CHECK_CLOSE(0., Residual(H_[node], Phi), eps);
			}
		}
	}

	TEST_FIXTURE (HelperFactory, MLOInit) {
		MLOcd M;
		M.push_back(leafX_, 0);
		M.push_back(leafX_, 1);
		M.push_back(leafX_, tree_.nLeaves() - 1);
		M.push_back(leafI_, tree_.nLeaves() - 1);

		TensorOperatorTree H(M, optree_);

		const Leaf& leaf = tree_.GetLeaf(tree_.nLeaves() - 1);
		const auto& node = (const Node&) leaf.Up();
		const Tensorcd& Phi = H[node];
			CHECK_CLOSE(0., abs(Phi(0)), eps);
			CHECK_CLOSE(1., abs(Phi(1)), eps);
			CHECK_CLOSE(1., abs(Phi(2)), eps);
			CHECK_CLOSE(0., abs(Phi(3)), eps);

			CHECK_CLOSE(1., abs(Phi(0, 1)), eps);
			CHECK_CLOSE(0., abs(Phi(1, 1)), eps);
			CHECK_CLOSE(0., abs(Phi(2, 1)), eps);
			CHECK_CLOSE(1., abs(Phi(3, 1)), eps);
	}
*/
}
