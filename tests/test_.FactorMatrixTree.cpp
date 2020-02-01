//
// Created by Roman Ellerbrock on 2020-01-27.
//
#include "Tree/FactorMatrixTree.h"
#include "Tree/HoleMatrixTree.h"
#include "UnitTest++/UnitTest++.h"

double eps = 1e-12;

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

	TEST (HoleMatrixTree_IO) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.dat");
		TensorTreecd T(filename);
		FactorMatrixTreecd S(T, T, basis);
		HoleMatrixTreecd Rho(T, T, S, basis);
		string file2("HO.tmp.dat");
		Rho.Write(file2);
		HoleMatrixTreecd Qh(file2);
			CHECK_EQUAL(Rho.size(), Qh.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(Rho[node], Qh[node]);
		}
	}

	TEST (HoleMatrixTree_Calc_persistance) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.dat");
		TensorTreecd T(filename);
		FactorMatrixTreecd S(T, T, basis);
		HoleMatrixTreecd Rho(T, T, S, basis);
		string correct_filename("TT.Holeoverlap.dat");
		Rho.Write(correct_filename);
		HoleMatrixTreecd Sigma(correct_filename);
			CHECK_EQUAL(Rho.size(), Sigma.size());
		for (const Node& node : basis) {
			auto r = Residual(Rho[node], Sigma[node]);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (HoleMatrixTree_Calc_Crosscheck) {
		TensorTreeBasis basis(12, 2, 2);
		mt19937 gen(1993);
		TensorTreecd Psi(basis, gen);
		HoleMatrixTreecd Rho(Psi, basis);
		FactorMatrixTreecd S(Psi, Psi, basis);
		HoleMatrixTreecd Rho2(Psi, Psi, S, basis);
			CHECK_EQUAL(Rho.size(), Rho2.size());
		for (const Node& node : basis) {
			auto d = abs(Rho[node](0, 0));
			if (!node.IsToplayer()) {
					CHECK_CLOSE(1., d, eps);
				auto r = Residual(Rho[node], Rho2[node]);
					CHECK_CLOSE(0., r, eps);
			}
		}
	}
}
