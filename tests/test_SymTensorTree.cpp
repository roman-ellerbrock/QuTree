//
// Created by Roman Ellerbrock on 1/31/21.
//

#include <gtest/gtest.h>
#include "TreeClasses/SymTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Extension.h"
#include "TreeClasses/SpectralDecompositionTree.h"

/**
 * Rationale: Test Suite for Symmetric Tensor Trees
 * We are checking whether
 * - all contractions of weighted tensors are diagonal
 * - upwards normalized tensors are orthonormal
 * - downwards normalized tensors are orthonormal
 * - sqrt(rho) * up/down-normalized tensors equal weighted tensors
 * - regularization: setting unoccupied orbitals without changing the TT
 */
double feps = 1e-7;

class TTFactory : public ::testing::Test {
protected:
    TTFactory() {
        tree_ = TreeFactory::balancedTree(10, 4, 3);
        mt19937 gen(34676949);
        psi_ = SymTensorTree(gen, tree_, false);
        chi_ = SymTensorTree(gen, tree_, true);

        /// Operator initialization
        auto I = &LeafInterface::identity;
        for (size_t l = 0; l < tree_.nLeaves(); ++l) { I_.push_back(I, l); }
        stree_ = make_shared<SparseTree>(SparseTree(I_, tree_, false));
        SparseMatrixTreecd x1(stree_, tree_);
        SparseMatrixTreecd x2(stree_, tree_);
    }

    ~TTFactory() = default;

    Tree tree_;
    SymTensorTree psi_;
    SymTensorTree chi_;

    MLOcd I_;
    shared_ptr<SparseTree> stree_;
};

TEST (SymTensorTree, Generate) {
    mt19937 gen(34676949);
    Tree tree = TreeFactory::balancedTree(10, 2, 3);
    SymTensorTree psi(gen, tree, true);
}

TEST (SymTensorTree, RegularizeTensor) {
    double delta = 1E-6;
    mt19937 gen(1283);
    TensorShape shape({2, 3, 4});
    Tensorcd A(shape);
    Tensor_Extension::generate(A, gen);
    for (size_t k = 0; k < shape.lastBefore(); ++k) {
        A(k) = 0.;
    }
    auto s = contraction(A, A, shape.lastIdx());
    ASSERT_NEAR(0, abs(s(0, 0)), delta);

    A = Tensor_Extension::regularize(A, shape.lastIdx(), delta);
    s = contraction(A, A, A.shape().lastIdx());
    ASSERT_EQ(true, (sqrt(abs(s(0, 0))) > delta / 2.));
}

TEST (SymTensorTree, DISABLED_NormalizeTensor) {
    mt19937 gen(1283);
    TensorShape shape({2, 3, 4});
    Tensorcd A(shape);
    Tensor_Extension::generate(A, gen);

    for (size_t k = 0; k < shape.order(); ++k) {
        auto Anorm = Tensor_Extension::normalize(A, k, feps);
        auto s = contraction(Anorm, Anorm, k);
        ASSERT_NEAR(0., residual(s, identityMatrixcd(s.dim1())), feps);
    }
}

/*	TEST_F (TTFactory, Normalization) {
		MatrixTreecd Sup(tree_); // wazzzz suuuuuuuup???!!!

		TreeFunctions::contractionUp(Sup, psi_, psi_, tree_);
		for (const Node& node : tree_) {
			if (!node.isToplayer()) {
				ASSERT_NEAR(0., residual(Sup[node], identityMatrixcd(Sup[node].dim1())), feps);
			}
		}

		MatrixTreecd Sdown(tree_); // wazzzz suuuuuuuup???!!!
		TreeFunctions::contractionDown(Sdown, psi_, psi_, Sup, tree_);
		for (const Node& node : tree_) {
			if (!node.isToplayer() && !node.isBottomlayer()) {
#					ASSERT_NEAR(0., residual(Sdown[node], identityMatrixcd(Sdown[node].dim1())), feps);
			}
		}
	}*/

TEST_F (TTFactory, DISABLED_Orthogonal) {
    mt19937 gen(1239);
    SymTensorTree psi(gen, tree_, false);
    SymTensorTree opsi = psi;
    opsi.orthogonal(tree_);
    auto s = TreeFunctions::dotProduct(psi, psi, tree_);
    for (auto x : s) {
        ASSERT_NEAR(0., abs(x - 1.), feps);
    }
}

TEST_F (TTFactory, Convert_checkUp) {
    mt19937 gen(1239);
    TensorTreecd Psi(gen, tree_, false);

    SymTensorTree spsi(Psi, tree_);
    TensorTreecd Psi2 = spsi.up_;
    const Node& top = tree_.topNode();
    Psi2[top] = spsi.weighted_[top];
    ASSERT_NEAR(0., TreeFunctions::residual(Psi, Psi2, tree_), feps);
}

TEST_F (TTFactory, Convert_checkDown) {
    mt19937 gen(1239);
    TensorTreecd Psi(gen, tree_, false);

    SymTensorTree spsi(Psi, tree_);
}

TEST_F (TTFactory, Apply) {
}

