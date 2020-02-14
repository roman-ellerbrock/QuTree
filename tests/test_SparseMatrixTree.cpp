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
		Tree tree_;
		TensorTreecd Psi_;
		MLOcd M_;

		void Initialize() {

			rng_ = mt19937(1993);
			tree_ = Tree(8, 2, 2);

			Psi_ = TensorTreecd(tree_, rng_);

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
		Tree tree(7, 4, 2);
		vector<size_t> modes({3, 4});
		SubTree active(modes, tree);
			CHECK_EQUAL(7, active.size());
	}

	TEST_FIXTURE (HelperFactory, TreeMarker_NoTail) {
		/// Create TreeMarker omitting higher nodes in the tree after last branch
		SubTree active(M_.Modes(), tree_, false);
			CHECK_EQUAL(5, active.size());
	}

	TEST_FIXTURE (HelperFactory, Represent) {
		SparseMatrixTreecd hmat = SparseMatrixTreeFunctions::Represent(M_, Psi_, tree_);
			CHECK_CLOSE(0.25, real(hmat[tree_.TopNode()](0, 0)), 0E-12);
	}

	TEST_FIXTURE (HelperFactory, Contraction) {
		SparseMatrixTreecd mats = SparseMatrixTreeFunctions::Represent(M_, Psi_, tree_);
		SparseMatrixTreecd holes(M_, tree_);
		SparseMatrixTreeFunctions::Contraction(holes, Psi_, Psi_, mats, tree_);
		for (const Node& node : tree_) {

		}
	}

	TEST_FIXTURE (HelperFactory, Constructor) {
		SparseMatrixTreecd hmat(M_, tree_);
			CHECK_EQUAL(6, hmat.Size());
	}

	TEST_FIXTURE (HelperFactory, IO) {
		SparseMatrixTreecd hmat(M_, tree_);
		hmat.Write("SparseMatrixTree.dat");
		SparseMatrixTreecd gmat(M_, tree_);
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

