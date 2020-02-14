//
// Created by Roman Ellerbrock on 2/12/20.
//

#include "UnitTest++/UnitTest++.h"
#include "Tree/MatrixTree.h"
#include "Tree/SpectralDecompositionTree.h"
#include "Tree/MatrixTreeFunctions.h"
#include "Util/RandomMatrices.h"

SUITE (MatrixTree) {
	double eps = 1e-8;

	TEST (Constructor) {
		Tree tree(7, 5, 4);
		MatrixTreecd M(tree);
			CHECK_EQUAL(tree.nNodes(), M.size());
	}

	TEST (IO) {
		Tree tree(7, 5, 4);
		MatrixTreecd M(tree);
		mt19937 gen(1988);
		for (const Node& node : tree) {
			Matrixcd& m = M[node];
			m = RandomMatrices::GUE(m.Dim1(), gen);
		}

		string filename("MatrixTree.dat");
		M.Write(filename);
		MatrixTreecd Mread(filename);
		for (const Node& node : tree) {
			double delta = Residual(M[node], Mread[node]);
				CHECK_CLOSE(0., delta, eps);
		}
	}
}

SUITE (MatrixTreeFunctions) {
	using namespace MatrixTreeFunctions;

	double eps = 1e-8;

	TEST (DotProduct) {
		mt19937 gen(1923);
		Tree tree(7, 5, 4);
		TensorTreecd Psi(tree, gen);
		MatrixTreecd S = DotProduct(Psi, Psi, tree);
		for (const Node& node : tree) {
			const Matrixcd& s = S[node];
			Matrixcd Identity = IdentityMatrix<complex<double>>(s.Dim1());
			double r = Residual(Identity, s);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (Contraction) {
		mt19937 gen(1923);
		Tree tree(7, 5, 4);
		TensorTreecd Psi(tree, gen, false);
		TensorTreecd Chi(tree, gen, false);
		MatrixTreecd S = DotProduct(Psi, Chi, tree);
		MatrixTreecd Rho = Contraction(Psi, Chi, S, tree);
			CHECK_EQUAL(tree.nNodes(), Rho.size());
			CHECK_EQUAL(tree.nNodes(), S.size());
	}

	TEST (Density) {
		mt19937 gen(1923);
		Tree tree(7, 5, 4);
		TensorTreecd Psi(tree, gen, true);
		MatrixTreecd Rho = Contraction(Psi, tree, true);
		for (const Node& node : tree) {
			if (!node.IsToplayer()) {
				Matrixcd& rho = Rho[node];
					CHECK_EQUAL(rho.Dim2(), rho.Dim1());
				for (size_t j = 0; j < rho.Dim2(); ++j) {
					for (size_t i = 0; i < rho.Dim1(); ++i) {
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
		Tree tree(12, 2, 2);
		TensorTreecd Psi(tree, gen);
		MatrixTreecd Rho = MatrixTreeFunctions::Contraction(Psi, tree, true);
		SpectralDecompositionTreecd X(Rho, tree);
			CHECK_EQUAL(Rho.size(), X.size());
		for (const Node& node : tree) {
			if (!node.IsToplayer()) {
					CHECK_CLOSE(1., X[node].second(1), eps);
			}
		}
	}

	TEST (SpectralDecompositionTree_Inverse) {
		Tree tree(12, 4, 2);
		mt19937 gen(1993);
		MatrixTreecd H(tree);
		for (const Node& node : tree) {
			const TensorDim& dim = node.TDim();
			Matrixcd mat = RandomMatrices::GUE(dim.GetNumTensor(), gen);
			H[node] = FactorMatrixcd(mat, node.ChildIdx());
		}

		SpectralDecompositionTreecd X(H, tree);
		auto H_inv = X.Invert(tree, 1e-10);

		MatrixTreecd Identity(tree);
		for (const Node& node : tree) {
			Identity[node] = H_inv[node] * H[node];
		}
		for (const Node& node : tree) {
			const Matrixcd& I_test = Identity[node];
			auto r = Residual(I_test, IdentityMatrix<complex<double>>(I_test.Dim1()));
				CHECK_CLOSE(0., r, eps);
		}
	}

	// @TODO: Move Unit tests from FactorMatrixTree to here.
}
