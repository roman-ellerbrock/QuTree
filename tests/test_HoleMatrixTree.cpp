//
// Created by Roman Ellerbrock on 2/10/20.
//
#include "UnitTest++/UnitTest++.h"
#include "IOTree.h"

SUITE(HoleMatrixTree) {
	double eps = 1e-7;

	TEST (HoleMatrixTree_IO) {
		TensorTreeBasis basis(12, 2, 2);
		string filename("TT.RNG.dat");
		TensorTreecd T(filename);
		FactorMatrixTreecd S(T, T, basis);
		cout << "CalculateRho" << endl;
		HoleMatrixTreecd Rho(T, T, S, basis);
		cout << "CalculateedRho" << endl;
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

/*	TEST(IO) {
		mt19937 gen(1993);
		TTBasis basis(14, 6, 4);
		TensorTreecd Psi(basis, gen, false);
		IOTree::Occupancy(Psi, basis);
		HoleMatrixTreecd Rho(Psi, basis);
		IOTree::Leafs(Psi, Rho, basis);
	}*/


}


