//
// Created by Roman Ellerbrock on 2/12/20.
//

#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/MatrixTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Util/RandomMatrices.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/TreeTransformation.h"
#include "TreeClasses/TensorTreeFunctions.h"

SUITE (MatrixTree) {
	double eps = 1e-8;

	TEST (Constructor) {
		Tree tree = TreeFactory::balancedTree(7, 5, 4);
		MatrixTreecd M(tree);
			CHECK_EQUAL(tree.nNodes(), M.size());
	}

	TEST (IO) {
		Tree tree = TreeFactory::balancedTree(7, 5, 4);
		MatrixTreecd M(tree);
		mt19937 gen(1988);
		for (const Node& node : tree) {
			Matrixcd& m = M[node];
			m = RandomMatrices::gue(m.dim1(), gen);
		}

		string filename("MatrixTree.dat");
		M.write(filename);
		MatrixTreecd Mread(filename);
		for (const Node& node : tree) {
			double delta = residual(M[node], Mread[node]);
				CHECK_CLOSE(0., delta, eps);
		}
	}
}

SUITE (MatrixTreeFunctions) {
	using namespace TreeFunctions;

	double eps = 1e-8;

	TEST (dotProduct) {
		mt19937 gen(1923);
		Tree tree = TreeFactory::balancedTree(7, 5, 4);
		TensorTreecd Psi(gen, tree);
		MatrixTreecd S = dotProduct(Psi, Psi, tree);
		for (const Node& node : tree) {
			const Matrixcd& s = S[node];
			Matrixcd Identity = identityMatrix<complex<double>>(s.dim1());
			double r = residual(Identity, s);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (contraction) {
		mt19937 gen(1923);
		Tree tree = TreeFactory::balancedTree(7, 5, 4);
		TensorTreecd Psi(gen, tree, false);
		TensorTreecd Chi(gen, tree, false);
		MatrixTreecd S = dotProduct(Psi, Chi, tree);
		MatrixTreecd Rho = contraction(Psi, Chi, S, tree);
			CHECK_EQUAL(tree.nNodes(), Rho.size());
			CHECK_EQUAL(tree.nNodes(), S.size());
	}

	TEST (Density) {
		mt19937 gen(1923);
		Tree tree = TreeFactory::balancedTree(7, 5, 4);
		TensorTreecd Psi(gen, tree, true);
		MatrixTreecd Rho = contraction(Psi, tree, true);
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
				Matrixcd& rho = Rho[node];
					CHECK_EQUAL(rho.dim2(), rho.dim1());
				for (size_t j = 0; j < rho.dim2(); ++j) {
					for (size_t i = 0; i < rho.dim1(); ++i) {
						if (j == 0 && i == 0) {
								CHECK_CLOSE(1., abs(rho(j, i)), eps);
						} else {
								CHECK_CLOSE(0., abs(rho(j, i)), eps);
						}
					}
				}
			}
		}
	}

	TEST (SpectralDecompositionTree_Calc) {
		mt19937 gen(1993);
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		TensorTreecd Psi(gen, tree);
		MatrixTreecd Rho = TreeFunctions::contraction(Psi, tree, true);
		SpectralDecompositionTreecd X(Rho, tree);
			CHECK_EQUAL(Rho.size(), X.size());
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
					CHECK_CLOSE(1., X[node].second(1), eps);
			}
		}
	}

	TEST (SpectralDecompositionTree_Inverse) {
		Tree tree = TreeFactory::balancedTree(12, 4, 2);
		mt19937 gen(1993);
		MatrixTreecd H(tree);
		for (const Node& node : tree) {
			const TensorShape& dim = node.shape();
			auto mat = RandomMatrices::gue(dim.lastDimension(), gen);
			auto mat_dagger = mat.adjoint();
			H[node] = mat * mat_dagger;
		}

		SpectralDecompositionTreecd X(H, tree);
		auto H_inv = X.invert(tree, 1e-10);

		MatrixTreecd Identity(tree);
		for (const Node& node : tree) {
			Identity[node] = H_inv[node] * H[node];
		}
		for (const Node& node : tree) {
			const Matrixcd& I_test = Identity[node];
			auto r = residual(I_test, identityMatrix<complex<double>>(I_test.dim1()));
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (canonicalTransformation) {

		Tree tree = TreeFactory::balancedTree(12, 4, 2);
		mt19937 gen(1993);
		TensorTreecd Psi(gen, tree);

		canonicalTransformation(Psi, tree, true);
		auto rho = TreeFunctions::contraction(Psi, tree, true);
		double off = 0.;
		for (const auto& mat : rho) {
			for (size_t j = 0; j < mat.dim2(); ++j) {
				for (size_t i = 0; i < mat.dim1(); ++i) {
					if (i != j) { off += abs(mat(i, j)); }
				}
			}
		}
			CHECK_CLOSE(0., off, eps);
	}

	TEST (CanonicalTransformation2) {

		Tree tree = TreeFactory::balancedTree(12, 4, 2);
		mt19937 gen(1993);
		TensorTreecd Psi(gen, tree, false);
		auto rho = TreeFunctions::contraction(Psi, tree, true);
		canonicalTransformation(Psi, tree, true);
		rho = TreeFunctions::contraction(Psi, tree, true);

		auto S = TreeFunctions::dotProduct(Psi, Psi, tree);

		rho = TreeFunctions::contraction(Psi, tree, true);
		double off = 0.;
		for (const auto& mat : rho) {
			for (size_t j = 0; j < mat.dim2(); j++) {
				for (size_t i = 0; i < mat.dim1(); ++i) {
					if (i != j) { off += abs(mat(i, j)); }
				}
			}
		}
			CHECK_CLOSE(0., off, eps);
	}

}

SUITE(TreeTransformations) {
	double eps = 1e-8;

	TEST(ContractionNormalized) {
		Tree tree = TreeFactory::balancedTree(12, 4, 2);
		// Increase number of states
		Node& top = tree.topNode();
		TensorShape& shape = top.shape();
		shape[shape.lastIdx()] += 1;
		tree.update();
		mt19937 gen(1993);
		TensorTreecd Psi(gen, tree);

		auto edgePsi = TreeFunctions::contractionNormalization(Psi, tree, true);

		for (const Edge& e : tree.edges()) {
			auto phi = edgePsi[e];
			Matrixcd deltaij = contraction(phi, phi, e.upIdx());
			Matrixcd I = identityMatrix<complex<double>>(deltaij.dim1());
			auto r = residual(deltaij, I);
			CHECK_CLOSE(0., r, 1e-7);
		}
	}

	TEST(CompressTensorTree) {
		Tree tree = TreeFactory::balancedTree(12, 4, 2);
		mt19937 gen(1993);
		TensorTreecd Psi(gen, tree);
		canonicalTransformation(Psi, tree, true);
		MatrixTreecd Rho = TreeFunctions::contraction(Psi, tree, true);
		SpectralDecompositionTreecd X(Rho, tree);

		TreeFunctions::adjust(Psi, tree, X, 1e-7);
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			/// Check Node TensorShape
			if (!node.isBottomlayer()) {
				CHECK_EQUAL(1, shape.totalDimension());
			} else {
				CHECK_EQUAL(4, shape.totalDimension());
			}
			/// Check Tensor
			Tensorcd Phiacc(shape);
			Phiacc(0) = 1.;
				CHECK_CLOSE(0., residual(Psi[node], Phiacc), eps);
		}
	}
}