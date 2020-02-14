//
// Created by Roman Ellerbrock on 2/13/20.
//
#include "UnitTest++/UnitTest++.h"
#include "SparseMatrixTreeFunctions.h"
#include "Util/RandomMatrices.h"

SUITE (SparseMatrixTree) {

	double eps = 1e-7;

	class HelperFactory {
	public:
		HelperFactory() {
			Initialize();
		}
		~HelperFactory() = default;
		mt19937 rng_;
		TTBasis basis_;
		TensorTreecd Psi_;
		MLOcd M_;

		void Initialize() {

			rng_ = mt19937(1993);
			basis_ = TTBasis(8, 2, 2);

			Psi_ = TensorTreecd(basis_, rng_);

			// Generate an bit-flip operator and Fmatrix
			FactorMatrixcd X(2, 1);
			X(0, 0) = 0.5;
			X(1, 1) = 0.5;
			LeafMatrixcd x(X);
			M_ = MLOcd(x, 0);
			M_.push_back(x, 3);
		}
	};

	TEST (TreeMarker) {
		TTBasis basis(7, 4, 2);
		vector<size_t> modes({3, 4});
		SubTree active(modes, basis);
			CHECK_EQUAL(7, active.size());
	}

	TEST_FIXTURE (HelperFactory, TreeMarker_NoTail) {
		/// Create TreeMarker omitting higher nodes in the tree after last branch
		SubTree active(M_.Modes(), basis_, false);
			CHECK_EQUAL(5, active.size());
	}

	TEST_FIXTURE (HelperFactory, Represent) {
		SparseMatrixTreecd hmat = SparseMatrixTreeFunctions::Represent(M_, Psi_, basis_);
			CHECK_CLOSE(0.25, real(hmat[basis_.TopNode()](0, 0)), 0E-12);
	}

	TEST_FIXTURE (HelperFactory, Contraction) {
		SparseMatrixTreecd mats = SparseMatrixTreeFunctions::Represent(M_, Psi_, basis_);
		SparseMatrixTreecd holes(M_, basis_);
		SparseMatrixTreeFunctions::Contraction(holes, Psi_, Psi_, mats, basis_);
		for (const Node& node : basis_) {

		}
	}

	TEST_FIXTURE (HelperFactory, Constructor) {
		SparseMatrixTreecd hmat(M_, basis_);
			CHECK_EQUAL(6, hmat.Size());
	}

	TEST_FIXTURE (HelperFactory, IO) {
		SparseMatrixTreecd hmat(M_, basis_);
		hmat.Write("SparseMatrixTree.dat");
		SparseMatrixTreecd gmat(M_, basis_);
		gmat.Read("SparseMatrixTree.dat");
		const SubTree& marker = gmat.Active();
			CHECK_EQUAL(hmat.Size(), gmat.Size());
		for (size_t i = 0; i < marker.size(); ++i) {
			const Node& node = marker.MCTDHNode(i);
			double r = Residual(hmat[node], gmat[node]);
				CHECK_CLOSE(0., r, eps);
		}
	}
}

