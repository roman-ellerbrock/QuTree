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
		TTFactory()
		{
			tree_ = TreeFactory::BalancedTree(10, 2, 3);
			mt19937 gen(34676949);
			psi_ = TensorTreecd(gen, tree_, true);
			chi_ = TensorTreecd(gen, tree_, false);
		}

		~TTFactory() =
		default;

		Tree tree_;
		TensorTreecd psi_;
		TensorTreecd chi_;
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
		CHECK_CLOSE(0., Residual(A, A2), eps);
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
					auto rho = Contraction(w, w, k);
					for (size_t l = 0; l < rho.Dim1(); ++l) {
						rho(l, l) = 0.;
					}
					CHECK_CLOSE(0., rho.FrobeniusNorm(), eps);
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
			auto s = phi.DotProduct(phi);
			CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim2())), eps);
		}

		for (const Node& node : tree_) {
			if (node.isToplayer()) { continue; }
			const Tensorcd& phi = sPsi.down_[node];
			auto s = Contraction(phi, phi, node.childIdx());
			CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim2())), eps);
		}
	}

/*	TEST_FIXTURE (TTFactory, Weighting) {
		Aint gonna work: check eigenvalues instead?
		auto psis = {psi_, chi_};
		for (const auto& Psi : psis) {
			SymTensorTree sPsi(Psi, tree_);
			auto rho = TreeFunctions::Contraction(Psi, tree_, true);

			/// Check Down (rho-comparison)
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Node& parent = node.parent();
				const Tensorcd& w = sPsi.weighted_[parent];
				auto rho2 = Contraction(w, w, node.childIdx());
				CHECK_CLOSE(0., Residual(rho[node], rho2), eps);
			}

			/// Check up
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Tensorcd& w = sPsi.weighted_[node];
				auto rho2 = Contraction(w, w, node.parentIdx());
				CHECK_CLOSE(0., Residual(rho[node], rho2), eps);
			}
		}
	}
*/
	TEST_FIXTURE (TTFactory, canonicalRepresentation) {
		/// Check Bottom-up overlap
		auto psis = {psi_, chi_};
//		auto psis = {psi_};
		for (const auto& psi : psis) {
			SymTensorTree sPsi(psi, tree_);
			/// Check upwards coherence: sqrt(rho)*up_ = weighted
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Tensorcd& w = sPsi.weighted_[node];
				auto rho = Contraction(w, w, node.parentIdx());
				auto B = calculateB(sPsi.weighted_[node], node.parentIdx());
//				auto w2 = MBatrixTensor(B, sPsi.up_[node], node.parentIdx());
				auto w2 = TensorMatrix(sPsi.up_[node], B, node.parentIdx());
					CHECK_CLOSE(0., Residual(sPsi.weighted_[node], w2), eps);
			}

			/// Check downwards coherence: sqrt(rho)*down_ = weighted
			for (const Node& node : tree_) {
				if (node.isToplayer()) { continue; }
				const Node& parent = node.parent();
				auto B = calculateB(sPsi.weighted_[parent], node.childIdx());
//				auto w2 = MatrixTensor(B, sPsi.down_[node], node.childIdx());
				auto w2 = TensorMatrix(sPsi.down_[node], B, node.childIdx());
				CHECK_CLOSE(0., Residual(sPsi.weighted_[parent], w2), eps);
			}
		}
	}

	TEST_FIXTURE(TTFactory, Overlap) {
		auto psis = {chi_};
//		auto psis = {psi_, chi_};

		for (auto psi : psis) {
			CanonicalTransformation(psi, tree_, true);
			SymTensorTree sPsi(psi, tree_);
			auto psiup = sPsi.bottomUpNormalized(tree_);
			cout << "psi:\n";
			psi.print(tree_);
			cout << "sPsi.bottomUp:\n";
			//psiup.print(tree_);
			for (const Node& node : tree_) {
				node.info();
				sPsi.weighted_[node].print();
			}
			auto S = TreeFunctions::DotProduct(psi, psiup, tree_);
			getchar();
		}
	}

/*	TEST (normalizeTensor) {
		TensorShape shape({2, 3, 4});
		mt19937 gen(19239123);
		/// Generate Random complex tensor
		Tensorcd A(shape);
		uniform_real_distribution<double> dist(-1., 1.);
		for (size_t real = 0; real < 2; ++real) {
			cout << "REAL == " << real << endl;
			for (size_t I = 0; I < shape.totalDimension(); ++I) {
				double x = dist(gen);
				double y = dist(gen);
				if (real == 0) {
					A(I) = x;
				} else {
					A(I) = complex<double>(x, y);
				}
			}
			/// Test for one contraction
			size_t k = 1;
			Matrixcd B = calculateB(A, k);
			{
				/// 1.) Test SLE
				/// At the moment this solves Anorm * B = A for Anorm.
				/// Should I solve: B * Anorm = A?
				auto Anorm = solveSLE(B, A, k);
				auto s = Contraction(Anorm, Anorm, k);
					CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim1())), eps);
				auto A2 = TensorMatrix(Anorm, B, k); /// Do I want this ot the version above?
				auto r = Residual(A, A2);
					CHECK_CLOSE(0., r, eps);
			}
			{
				/// 2.) Test Unitary SLE
				/// At the moment this solves Anorm * B = A for Anorm.
				/// Should I solve: B * Anorm = A?
				auto Anorm = normalizedTensorSVD(B, A, k);
				auto s = Contraction(Anorm, Anorm, k);
					CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim1())), eps);
//				auto A2 = TensorMatrix(Anorm, B, k); /// Do I want this ot the version above?
				auto A2 = MatrixTensor(B, Anorm, k); /// Do I want this ot the version above?
				auto r = Residual(A, A2);
					CHECK_CLOSE(0., r, eps);
			}
			canonicalTensor(A);
		}
	}

	TEST (CreateSym) {
		Tree tree = TreeFactory::BalancedTree(10, 2, 3);
		mt19937 gen(34676949);
		for (size_t l = 0; l < 2; ++l) {
			bool delta = true;
			if (l) delta = false;
			TensorTreecd Psi(gen, tree, delta);

			SymTensorTree sPsi(Psi, tree);
			/// Check top-down against QR scheme
			for (const Node& node : tree) {
				/// Top-down normalized
				if (!node.isToplayer()) {
					const Node& parent = node.parent();
					const auto& w = sPsi.weighted_[parent];
					auto Anorm = normalizedTensor(w, node.childIdx());
					auto u = Contraction(Anorm, sPsi.down_[node], node.childIdx());
					auto s = u * u.Adjoint();
						CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim1())), eps);
				}
			}
			for (const Node& node : tree) {
				const auto& w = sPsi.weighted_[node];
				auto Anorm = normalizedTensor(w, node.parentIdx());
				auto u = Anorm.DotProduct(sPsi.up_[node]);
				auto s = u * u.Adjoint();
					CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim1())), eps);
			}

			/// Check that correct wavefunction comes back out again
			TensorTreecd Chi = sPsi.bottomUpNormalized(tree);
			auto S = TreeFunctions::DotProduct(Chi, Psi, tree);
				CHECK_CLOSE(0., Residual(S[tree.TopNode()], IdentityMatrixcd(1)), eps);

			/// Check Dot
			for (const Node& node : tree) {
				const auto& Phi = sPsi.weighted_[node];
				Matrixcd s = Phi.DotProduct(Phi);
				double abss = abs(s.Trace());
					CHECK_CLOSE(1., abss, eps);
			}
		}
	}
 */
}
