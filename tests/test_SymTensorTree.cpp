//
// Created by Roman Ellerbrock on 1/31/21.
//

#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/SymTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Functions.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Util/GateOperators.h"
#include "TreeClasses/TreeIO.h"

SUITE (SymTensorTree) {
	/**
	 * Rationale: Test Suite for Symmetric Tensor Trees
	 * We are checking whether
	 * - all contractions of weighted tensors are diagonal
	 * - upwards normalized tensors are orthonormal
	 * - downwards normalized tensors are orthonormal
	 * - sqrt(rho) * up/down-normalized tensors equal weighted tensors
	 * - regularization: setting unoccupied orbitals without changing the TT
	 */
	double eps = 1e-7;

	class TTFactory {
	public:
		TTFactory() {
			tree_ = TreeFactory::balancedTree(4, 2, 3);
			mt19937 gen(34676949);
			psi_ = SymTensorTree(gen, tree_, true);
			chi_ = SymTensorTree(gen, tree_, false);

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

	TEST (Generate) {
		mt19937 gen(34676949);
		Tree tree = TreeFactory::balancedTree(10, 2, 3);
		SymTensorTree psi(gen, tree, true);
	}

	TEST (NormalizeTensor) {
		mt19937 gen(1283);
		TensorShape shape({2, 3, 4});
		Tensorcd A(shape);
		Tensor_Extension::generate(A, gen);
		for (size_t i = 0; i < A.shape().lastBefore(); ++i) {
			A(i) = 0;
		}

		for (size_t k = 0; k < shape.order(); ++k) {
			auto Anorm = Tensor_Extension::normalize(A, k, 1e-7);
			auto s = contraction(Anorm, Anorm, k);
				CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim1())), eps);
		}
	}

	TEST_FIXTURE (TTFactory, Normalization) {
		MatrixTreecd Sup(tree_); // wazzzz suuuuuuuup???!!!

		TreeFunctions::contractionUp(Sup, psi_, psi_, tree_);
		for (const Node& node : tree_) {
			if (!node.isToplayer()) {
					CHECK_CLOSE(0., residual(Sup[node], identityMatrixcd(Sup[node].dim1())), eps);
			}
		}

		for (const Node& node : tree_) {
			if (!node.isToplayer()) {
				Matrixcd s = contraction(psi_.down_[node], psi_.down_[node], node.childIdx());
					CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim1())), eps);
			}
		}

		MatrixTreecd Sdown(tree_);
		TreeFunctions::contractionDown(Sdown, psi_, psi_, Sup, tree_);
		for (const Node& node : tree_) {
			if (!node.isToplayer()) {
					CHECK_CLOSE(0., residual(Sdown[node], identityMatrixcd(Sdown[node].dim1())), eps);
			}
		}
	}

	TEST_FIXTURE (TTFactory, Represent) {
		SOPcd cnot = CNot(0, tree_.nLeaves() - 1);
		SymMatrixTrees mats(cnot, tree_);
		TreeFunctions::symRepresent(mats, psi_, psi_, cnot, tree_);
//		mats.print(tree_);
	}

	TEST_FIXTURE(TTFactory, Apply) {
		SOPcd cnot = CNot(0, tree_.nLeaves() - 1);
		SymMatrixTrees mats(cnot, tree_);
		SymTensorTree HPsi = psi_;
		double delta = 1e-7;
//		TreeIO::output(psi_.toConventional(tree_), tree_);

		TreeFunctions::applySCF(HPsi, mats, psi_, cnot, tree_, delta, 1, &cout);
//		TreeIO::output(HPsi.toConventional(tree_), tree_);
//		HPsi.print(tree_);

	}
}
