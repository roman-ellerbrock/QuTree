//
// Created by Roman Ellerbrock on 2020-01-27.
//
#include "Core/RandomMatrices.h"
#include "Tree/FactorMatrixTree.h"
#include "Tree/HoleMatrixTree.h"
#include "Tree/SpectralDecompositionTree.h"
#include "UnitTest++/UnitTest++.h"

double eps = 1e-7;

SUITE (TensorTreeOverlaps) {

	TEST (FactorMatrixTree_IO) {
		TensorTreeBasis basis(12, 2, 2);
		FactorMatrixTreecd S(basis);
		string file2("DO.tmp.dat");
		{
			ofstream os(file2);
			S.Write(os);
			os.close();
		}
		FactorMatrixTreecd Q(file2);
			CHECK_EQUAL(S.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(S[node], Q[node]);
		}
	}

	TEST (FactorMatrixTree_Calculate) {
		/// Check that Overlap is working by calculating norm of wavefunction
		TensorTreeBasis basis(12, 2, 2);
		TensorTreecd T(string("TT.RNG.dat"));
		FactorMatrixTreecd S(T, T, basis);
		const FactorMatrixcd& s = S.Get();
			CHECK_CLOSE(1., abs(s[0]), eps);
	}

	TEST (SpectralDecompositionTree_Calc) {
		TensorTreeBasis basis(12, 2, 2);
		mt19937 gen(1993);
		TensorTreecd Psi(basis, gen);
		HoleMatrixTreecd Rho(Psi, basis);
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
		HoleMatrixTreecd H(basis);
	 	for (const Node& node : basis) {
	 		const TensorDim& dim = node.TDim();
			Matrixcd mat = RandomMatrices::GUE(dim.GetNumTensor(), gen);
			H[node] = FactorMatrixcd(mat, node.ChildIdx());
	 	}
		SpectralDecompositionTreecd X(H, basis);
		auto H_inv = X.Invert(basis, 1e-10);

		HoleMatrixTreecd Identity(basis);
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
