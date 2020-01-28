//
// Created by Roman Ellerbrock on 2020-01-27.
//
#include "TensorTree/DenseOverlap.h"
#include "TensorTree/HoleOverlap.h"
#include "UnitTest++/UnitTest++.h"

double eps = 1e-12;

SUITE (TensorTreeOverlaps) {

	TEST (DenseOverlap_IO) {
		TensorTreeBasis basis(12, 2, 2);
		DenseOverlapcd S(basis);
		string file2("DO.tmp.dat");
		{
			ofstream os(file2);
			S.Write(os);
			os.close();
		}
		DenseOverlapcd Q(file2);
			CHECK_EQUAL(S.size(), Q.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(S[node], Q[node]);
		}
	}

	TEST (TensorTree_Overlap) {
		/// Check that Overlap is working by calculating norm of wavefunction
		TensorTreeBasis basis(12, 2, 2);
		TensorTreecd T(string("TT.RNG.dat"));
		DenseOverlapcd S(T, T, basis);
		const FactorMatrixcd& s = S.Get();
			CHECK_CLOSE(1., abs(s[0]), eps);
	}

	TEST (TensorTree_HoleOverlap) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		HoleOverlapcd Rho(T, T, S, basis);
		string correct_filename("TT.Holeoverlap.dat");
		Rho.Write(correct_filename);
		HoleOverlapcd Sigma(correct_filename);
			CHECK_EQUAL(Rho.size(), Sigma.size());
		for (const Node& node : basis) {
			auto r = Residual(Rho[node], Sigma[node]);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST (HoleOverlap_IO) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.dat");
		TensorTreecd T(filename);
		DenseOverlapcd S(T, T, basis);
		HoleOverlapcd Rho(T, T, S, basis);
		string file2("HO.tmp.dat");
		{
			ofstream os(file2);
			Rho.Write(os);
			os.close();
		}
		HoleOverlapcd Qh(file2);
			CHECK_EQUAL(Rho.size(), Qh.size());
		for (const Node& node : basis) {
				CHECK_EQUAL(Rho[node], Qh[node]);
		}
	}
}
