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
		TTBasis basis(7, 5, 4);
		MatrixTreecd M(basis);
			CHECK_EQUAL(basis.nNodes(), M.size());
	}

	TEST (IO) {
		TTBasis basis(7, 5, 4);
		MatrixTreecd M(basis);
		mt19937 gen(1988);
		for (const Node& node : basis) {
			Matrixcd& m = M[node];
			m = RandomMatrices::GUE(m.Dim1(), gen);
		}

		string filename("MatrixTree.dat");
		M.Write(filename);
		MatrixTreecd Mread(filename);
		for (const Node& node : basis) {
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
		TTBasis basis(7, 5, 4);
		TensorTreecd Psi(basis, gen);
		MatrixTreecd S = DotProduct(Psi, Psi, basis);
		for (const Node& node : basis) {
			const Matrixcd& s = S[node];
			Matrixcd Identity = IdentityMatrix<complex<double>>(s.Dim1());
			double r = Residual(Identity, s);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (Contraction) {
		mt19937 gen(1923);
		TTBasis basis(7, 5, 4);
		TensorTreecd Psi(basis, gen, false);
		TensorTreecd Chi(basis, gen, false);
		MatrixTreecd S = DotProduct(Psi, Chi, basis);
		MatrixTreecd Rho = Contraction(Psi, Chi, S, basis);
			CHECK_EQUAL(basis.nNodes(), Rho.size());
			CHECK_EQUAL(basis.nNodes(), S.size());
	}

	TEST (Density) {
		mt19937 gen(1923);
		TTBasis basis(7, 5, 4);
		TensorTreecd Psi(basis, gen, true);
		MatrixTreecd Rho = Contraction(Psi, basis, true);
		Rho.print(basis);
		for (const Node& node : basis) {
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
		TensorTreeBasis basis(12, 2, 2);
		TensorTreecd Psi(basis, gen);
		MatrixTreecd Rho = MatrixTreeFunctions::Contraction(Psi, basis, true);
		SpectralDecompositionTreecd X(Rho, basis);
			CHECK_EQUAL(Rho.size(), X.size());
		for (const Node& node : basis) {
			if (!node.IsToplayer()) {
					CHECK_CLOSE(1., X[node].second(1), eps);
			}
		}
	}

	TEST (SpectralDecompositionTree_Inverse) {
		TensorTreeBasis basis(12, 4, 2);
		mt19937 gen(1993);
		MatrixTreecd H(basis);
		for (const Node& node : basis) {
			const TensorDim& dim = node.TDim();
			Matrixcd mat = RandomMatrices::GUE(dim.GetNumTensor(), gen);
			H[node] = FactorMatrixcd(mat, node.ChildIdx());
		}
		SpectralDecompositionTreecd X(H, basis);
		auto H_inv = X.Invert(basis, 1e-10);

		MatrixTreecd Identity(basis);
		for (const Node& node : basis) {
			Identity[node] = H_inv[node] * H[node];
		}
		for (const Node& node : basis) {
			const Matrixcd& I_test = Identity[node];
			auto r = Residual(I_test, IdentityMatrix<complex<double>>(I_test.Dim1()));
				CHECK_CLOSE(0., r, eps);
		}
	}
}
