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
		mt19937 rng_;
		TTBasis basis_;
		TensorTreecd Psi_;
		MLOcd M_;

		void Initialize() {

			rng_ = mt19937(1993);
			basis_ = TTBasis(8, 2, 2);

			Psi_ = TensorTreecd (basis_, rng_);

			// Generate an bit-flip operator and Fmatrix
			FactorMatrixcd X(2, 1);
			X(0, 0) = 0.5;
			X(1, 1) = 0.5;
			LeafMatrixcd x(X);
			M_ = MLOcd (x, 0);
			M_.push_back(x, 3);
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
		SparseFactorMatrixTreecd hmat(Psi_, M_, basis_);
			CHECK_CLOSE(0.25, real(hmat[basis_.TopNode()](0, 0)), 0E-12);
	}

	TEST_FIXTURE (HelperFactory, HMatrices_IO) {
		Initialize();

		// Create a multiparticleoperator
		// Build  matrix representation for operator
		SparseFactorMatrixTreecd hmat(Psi_, M_, basis_);

		string filename("Hmat.dat");
		hmat.Write(filename);
		SparseFactorMatrixTreecd gmat(M_, basis_, filename);
			CHECK_EQUAL(hmat.Size(), gmat.Size());
		const auto& active = hmat.Active();
		for (const Node *node_ptr : active) {
			const Node& node = *node_ptr;
				CHECK_EQUAL(hmat[node], gmat[node]);
		}
	}

	TEST_FIXTURE (HelperFactory, HHoleMatrices) {
		Initialize();

		// Build  matrix representation for operator
		SparseFactorMatrixTreecd hmat(Psi_, M_, basis_);
		SparseHoleMatrixTreecd hhole(Psi_, hmat, M_, basis_);
		Psi_.Write("Psi.dat");

		TensorTreecd Chi("Psi.dat");
		string filename("HHole.dat");
		hhole.Write(filename);
		SparseHoleMatrixTreecd ghole(M_, basis_, filename);
			CHECK_EQUAL(hhole.Size(), ghole.Size());
		const auto& active = hmat.Active();
		for (const Node *node_ptr : active) {
			const Node& node = *node_ptr;
				CHECK_EQUAL(hhole[node], ghole[node]);
		}
	}

	TEST_FIXTURE (HelperFactory, TreeMarker_NoTail) {
		/// Create TreeMarker omitting higher nodes in the tree after last branch
		Initialize();

		SparseFactorMatrixTreecd hmat (Psi_, M_, basis_);
		SparseHoleMatrixTreecd hhole(Psi_, hmat, M_, basis_);
		TreeMarker active(cast_to_vector_size_t(M_.Modes()), basis_, false);
		CHECK_EQUAL(5, active.size());

		/// Calculate sparse
//		hhole.Calculate(Psi_, hmat, active);
	}
}


