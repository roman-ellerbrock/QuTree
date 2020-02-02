//
// Created by Roman Ellerbrock on 2020-01-27.
//
#include "SparseFactorMatrixTree.h"
#include "Tree/SparseHoleMatrixTree.h"
#include "UnitTest++/UnitTest++.h"


SUITE (HMatrices) {
	class HelperFactory {
	public:
		HelperFactory() = default;
		~HelperFactory() = default;
		FactorMatrixcd X;
		SPOMcd x;
		mt19937 rng;

		void Initialize() {
			// Generate an bit-flip operator and Fmatrix
			X = FactorMatrixcd(2, 1);
			X(0, 0) = 0.5;
			X(1, 1) = 0.5;
			x = SPOMcd(X);

			rng = mt19937(1993);
		}
	};

	TEST (TreeMarker) {
		TTBasis basis(7, 4, 2);
		vector<size_t> modes({3, 4});
		TreeMarker active(modes, basis);
			CHECK_EQUAL(7, active.size());
	}

	TEST_FIXTURE (HelperFactory, HMatrices_Calc) {
		Initialize();
		TTBasis basis(8, 2, 2);
		TensorTreecd Psi(basis, rng);
		MPOcd M(x, 0);
		M.push_back(x, 3);
		SparseFactorMatrixTreecd hmat(Psi, M, basis);
			CHECK_CLOSE(0.25, real(hmat[basis.TopNode()](0, 0)), 0E-12);
	}

	TEST_FIXTURE (HelperFactory, HMatrices_IO) {
		Initialize();

		// Create Basis
		TTBasis basis(7, 2, 2);
		// Create and occupy Tree
		TensorTreecd Psi(basis, rng, false);
		// Create a multiparticleoperator
		MPOcd M(x, 0);
		M.push_back(x, 3);
		// Build  matrix representation for operator
		SparseFactorMatrixTreecd hmat(Psi, M, basis);

		string filename("Hmat.dat");
		hmat.Write(filename);
		SparseFactorMatrixTreecd gmat(M, basis, filename);
			CHECK_EQUAL(hmat.Size(), gmat.Size());
		const auto& active = hmat.Active();
		for (const Node *node_ptr : active) {
			const Node& node = *node_ptr;
				CHECK_EQUAL(hmat[node], gmat[node]);
		}
	}

	TEST_FIXTURE (HelperFactory, HHoleMatrices) {
		Initialize();

		// Create Basis
		TTBasis basis(7, 2, 2);
		// Create and occupy Tree
		TensorTreecd Psi(basis, rng, false);
		// Create a multiparticleoperator acting on modes 0 and 3
		MPOcd M(x, 0);
		M.push_back(x, 3);
		// Build  matrix representation for operator
		SparseFactorMatrixTreecd hmat(Psi, M, basis);
		HoleMatrixTreecd hhole(Psi, hmat, M, basis);
		Psi.Write("Psi.dat");
		TensorTreecd Chi("Psi.dat");

		string filename("HHole.dat");
		hhole.Write(filename);
		HoleMatrixTreecd ghole(M, basis, filename);
			CHECK_EQUAL(hhole.Size(), ghole.Size());
		const auto& active = hmat.Active();
		for (const Node *node_ptr : active) {
			const Node& node = *node_ptr;
				CHECK_EQUAL(hhole[node], ghole[node]);
		}
	}
}


