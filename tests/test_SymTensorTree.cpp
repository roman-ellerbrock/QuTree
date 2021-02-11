//
// Created by Roman Ellerbrock on 1/31/21.
//

#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/SymTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Extension.h"
#include "TreeClasses/SpectralDecompositionTree.h"

SUITE (SymTensorTree) {
	double eps = 1e-7;

	class TTFactory {
		public:
		TTFactory() {
			tree_ = TreeFactory::BalancedTree(10, 2, 3);
			mt19937 gen(34676949);
			psi_ = TensorTreecd(gen, tree_, true);
			chi_ = TensorTreecd(gen, tree_, false);
			spsi_ = SymTensorTree(psi_, tree_);
			schi_ = SymTensorTree(chi_, tree_);

			/// Operator initialization
			auto I = &LeafInterface::Identity;
			for (size_t l = 0; l < tree_.nLeaves(); ++l) { I_.push_back(I, l); }
			stree_ = make_shared<SparseTree>(SparseTree(I_, tree_, false));
			SparseMatrixTreecd x1(stree_, tree_);
			SparseMatrixTreecd x2(stree_, tree_);
		}

		~TTFactory() = default;

		Tree tree_;
		TensorTreecd psi_;
		TensorTreecd chi_;
		SymTensorTree spsi_;
		SymTensorTree schi_;

		MLOcd I_;
		shared_ptr<SparseTree> stree_;
	};

	TEST (TensorFlatten) {
		TensorShape shape({2, 3, 4});
		mt19937 gen(19239123);
		Tensorcd A(shape);
		uniform_real_distribution<double> dist(-1., 1.);
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			double x = dist(gen);
			double y = dist(gen);
			A(I) = complex<double>(x, y);
		}

		size_t k = 1;
		auto Aflat = toMatrix(A, k);
		auto A2 = toTensor(Aflat, A.shape(), k);
		CHECK_CLOSE(0., residual(A, A2), eps);
	}

	TEST_FIXTURE(TTFactory, canonicalTensors) {
		auto psis = {psi_, chi_};
		/// Check whether all density matrices are diagonal
		for (const auto& Psi : psis) {
			SymTensorTree sPsi(Psi, tree_);
			/// Check down
			for (const Node& node : tree_) {
				const Tensorcd& w = sPsi.weighted_[node];
				if (node.isBottomlayer()) { continue; }
				for (size_t k = 0; k < node.nChildren(); ++k) {
					auto rho = contraction(w, w, k);
					for (size_t l = 0; l < rho.dim1(); ++l) {
						rho(l, l) = 0.;
					}
					CHECK_CLOSE(0., rho.frobeniusNorm(), eps);
				}
			}
		}
	}

	TEST_FIXTURE (TTFactory, normalized) {
		SymTensorTree sPsi(psi_, tree_);

		/// Check normalization
		for (const Node& node : tree_) {
			if (node.isToplayer()) { continue; }
			const Tensorcd& phi = sPsi.up_[node];
			auto s = phi.dotProduct(phi);
			CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim2())), eps);
		}

		for (const Node& node : tree_) {
			if (node.isToplayer()) { continue; }
			const Tensorcd& phi = sPsi.down_[node];
			auto s = contraction(phi, phi, node.childIdx());
			CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim2())), eps);
		}
	}

	TEST_FIXTURE (TTFactory, canonicalRepresentation) {
		/// Check Bottom-up overlap
		auto psis = {psi_, chi_};
		for (const auto& psi : psis) {
			SymTensorTree sPsi(psi, tree_);
			/// Check upwards coherence: sqrt(rho)*up_ = weighted
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Tensorcd& w = sPsi.weighted_[node];
				auto rho = contraction(w, w, node.parentIdx());
				auto B = calculateB(sPsi.weighted_[node], node.parentIdx());
//				auto w2 = MBatrixTensor(B, sPsi.up_[node], node.parentIdx());
				auto w2 = tensorMatrix(sPsi.up_[node], B, node.parentIdx());
					CHECK_CLOSE(0., residual(sPsi.weighted_[node], w2), eps);
			}

			/// Check downwards coherence: sqrt(rho)*down_ = weighted
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Node& parent = node.parent();
				auto B = calculateB(sPsi.weighted_[parent], node.childIdx());
//				auto w2 = MatrixTensor(B, sPsi.down_[node], node.childIdx());
				auto w2 = tensorMatrix(sPsi.down_[node], B, node.childIdx());
				CHECK_CLOSE(0., residual(sPsi.weighted_[parent], w2), eps);
			}
		}
	}

	TEST_FIXTURE(TTFactory, Overlap) {
		auto psis = {psi_, chi_};

		for (auto psi : psis) {
			CanonicalTransformation(psi, tree_, true);
			auto rho = TreeFunctions::Contraction(psi, tree_, true);
			SymTensorTree sPsi(psi, tree_);
			auto psiup = sPsi.bottomUpNormalized(tree_);
			auto S = TreeFunctions::DotProduct(psi, psiup, tree_);
			auto s = S[tree_.TopNode()];
				CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim1())), eps);
		}
	}

	TEST_FIXTURE(TTFactory, OverlapNonCanonical) {
		auto psis = {psi_, chi_};

		for (const auto& psi : psis) {
			auto rho = TreeFunctions::Contraction(psi, tree_, true);
			SymTensorTree sPsi(psi, tree_);
			auto psiup = sPsi.bottomUpNormalized(tree_);
			auto S = TreeFunctions::DotProduct(psi, psiup, tree_);
			auto s = S[tree_.TopNode()];
				CHECK_CLOSE(0., residual(s, identityMatrixcd(s.dim1())), eps);
		}
	}

	TEST_FIXTURE(TTFactory, symOverlap) {
		SymTensorTree spsi(psi_, tree_);
		SymTensorTree schi(chi_, tree_);
		auto S = TreeFunctions::symDotProduct(spsi, schi, tree_);
		auto s_top = S[tree_.TopNode()].trace();
			CHECK_CLOSE(-0.00557989, real(s_top), eps);
			CHECK_CLOSE(0., imag(s_top), eps);
		for (const Node& node : tree_) {
			CHECK_CLOSE(0., abs(s_top - S[node].trace()), eps);
		}
	}

	TEST_FIXTURE(TTFactory, symRepresent) {
		/**
		 * Rationale:
		 * - Represent operator using sym-routine and compare
		 *   to standard represent.
		 */
		SparseMatrixTreecd x1(stree_, tree_);
		SparseMatrixTreecd x2(stree_, tree_);
		SymMatrixTree mat({x1, x2});
//		TreeFunctions::symRepresent(mat, spsi_, schi_, I_, tree_);

//		auto hmat = TreeFunctions::Represent(I_, psi_, chi_, tree_);
		SparseMatrixTreecd hmat(stree_, tree_);
		auto S = TreeFunctions::DotProduct(psi_, chi_, tree_);
//		cout << "s:\n";
//		S.print(tree_);
//		getchar();
		TreeFunctions::Represent(hmat, I_, psi_, chi_, tree_);
//		SparseMatrixTreecd hhole(stree_, tree_);
//		TreeFunctions::Contraction(hhole, psi_, chi_, hmat, tree_);

		for (const Node* node_ptr : *stree_) {
			const Node& node = *node_ptr;
//			node.info();
//			mat.first[node].print();
//			hmat[node].print();
//			CHECK_CLOSE(0., residual(mat.first[node], hmat[node]), eps);
		}
	}

	TEST_FIXTURE(TTFactory, symRepresentOverlap) {
		/**
		 * Rationale:
		 * - Calculate overlap via representing identity-operators and
		 *   Applying them to wavefunction, then calculating overlap
		 */
		auto S = TreeFunctions::symDotProduct(spsi_, schi_, tree_);
		SparseMatrixTreecd x1(stree_, tree_);
		SparseMatrixTreecd x2(stree_, tree_);
		SymMatrixTree mat({x1, x2});
		TreeFunctions::symRepresent(mat, spsi_, schi_, I_, tree_);

		auto Hschi_ = schi_;
		for (const Node& node : tree_) {
			Hschi_.weighted_[node] = TreeFunctions::symApply(
				schi_.weighted_[node], mat, I_, node);
		}

		for (const Node& node : tree_) {
			auto s = spsi_.weighted_[node].dotProduct(Hschi_.weighted_[node]);
			CHECK_CLOSE(0., residual(s, S[node]), eps);
		}
	}

}
