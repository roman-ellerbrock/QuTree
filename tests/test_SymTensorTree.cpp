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

class TTFactory: public ::testing::Test {
protected:
	TTFactory() {
		tree_ = TreeFactory::balancedTree(10, 4, 3);
		gen_ = mt19937(34676949);
		psi_ = SymTensorTree(gen_, tree_, false);
		chi_ = SymTensorTree(gen_, tree_, true);

		/// Operator initialization
		auto I = &LeafInterface::identity;
		for (size_t l = 0; l < tree_.nLeaves(); ++l) { I_.push_back(I, l); }
		stree_ = make_shared<SparseTree>(SparseTree(I_, tree_, false));
		SparseMatrixTreecd x1(stree_, tree_);
		SparseMatrixTreecd x2(stree_, tree_);

		/// operator for testing matrix representations
		Matrixcd X(2, 2);
		X(0, 0) = 0.5;
		X(1, 0) = 0.5;
		X(0, 1) = -0.5;
		X(1, 1) = 0.5;
		for (size_t i = 0; i < tree_.nLeaves(); ++i) {
			LeafMatrixcd x(X);
			M_.push_back(x, i);
		}
	}

	~TTFactory() = default;

	Tree tree_;
	SymTensorTree psi_;
	SymTensorTree chi_;
	mt19937 gen_;

	MLOcd I_, M_;
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

TEST_F (TTFactory, Orthogonal) {
	mt19937 gen(1239);
	SymTensorTree psi(gen, tree_, false);
	auto s = TreeFunctions::dotProduct(psi, psi, tree_);
	for (auto x: s) {
		ASSERT_NEAR(0., abs(x - 1.), feps);
	}
}

TEST_F (TTFactory, mixedDotProductLocal) {
	TensorTreecd Psi(gen_, tree_);
	SymTensorTree psi(Psi, tree_);
	psi.up_ = Psi;
	psi.weighted_[tree_.topNode()] = Psi[tree_.topNode()];

	auto ssym = TreeFunctions::mixedDotProduct(Psi, psi, tree_);
	auto s = TreeFunctions::dotProduct(Psi, Psi, tree_);
	for (const Node& node : tree_) {
		ASSERT_NEAR(0., residual(ssym[node], s[node]), feps);
	}
}

TEST_F (TTFactory, mixedDotProduct) {
	TensorTreecd Psi(gen_, tree_);
	SymTensorTree psi(Psi, tree_);
	auto ssym = TreeFunctions::mixedDotProduct(Psi, psi, tree_);
		ASSERT_NEAR(0., residual(ssym[tree_.topNode()], identityMatrixcd(1)), feps);
}

TEST_F (TTFactory, mixedRho) {
	TensorTreecd Psi(gen_, tree_);
	SymTensorTree spsi(Psi, tree_);

	auto r = TreeFunctions::mixedRho(Psi, spsi, tree_);
	auto rho = TreeFunctions::contraction(Psi, tree_, true);

	for (const Node& node: tree_) {
		if (node.isToplayer()) { continue; }
		auto r2 = r[node] * r[node].adjoint();
		ASSERT_NEAR(0., residual(rho[node], r2), feps);
	}
}

TEST_F (TTFactory, symMatrices) {
	TensorTreecd Psi(gen_, tree_);
	SymTensorTree spsi(Psi, tree_);
	auto s = TreeFunctions::mixedDotProduct(Psi, spsi, tree_);

	SparseMatrixTreecd xmat(M_, tree_);
	SparseMatrixTreecd xhole(M_, tree_);
	SymMatrixTree sMatPair(xmat, xhole);


	TreeFunctions::represent(xmat, M_, Psi, Psi, tree_);

	TreeFunctions::symRepresent(sMatPair, spsi, spsi, M_, tree_);
	SparseMatrixTreecd smat = sMatPair.first;

	for (const Node& node : tree_) {
		if (node.isToplayer()) { continue; }
		smat[node] = s[node] * smat[node] * s[node].adjoint();
		ASSERT_NEAR(0., residual(xmat[node], smat[node]), feps);
	}

}

TEST_F (TTFactory, symHoles) {
	TensorTreecd Psi(gen_, tree_);
	SymTensorTree spsi(Psi, tree_);
	auto r = TreeFunctions::mixedRho(Psi, spsi, tree_);

	SparseMatrixTreecd xmat(M_, tree_);
	SparseMatrixTreecd xhole(M_, tree_);
	SymMatrixTree smatpair(xmat, xhole);

	TreeFunctions::represent(xmat, M_, Psi, Psi, tree_);
	TreeFunctions::contraction(xhole, Psi, Psi, xmat, tree_);

	TreeFunctions::symRepresent(smatpair, spsi, spsi, M_, tree_);
	SparseMatrixTreecd shole = smatpair.second;

	for (const Node& node : tree_) {
		if (node.isToplayer()) { continue; }
		shole[node] = r[node] * shole[node] * r[node].adjoint();
		ASSERT_NEAR(0., residual(shole[node], xhole[node]), feps);
	}

}
