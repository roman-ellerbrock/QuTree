//
// Created by Roman Ellerbrock on 2/13/20.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "Util/RandomMatrices.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/SparseTensorTree.h"

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
			tree_ = TreeFactory::balancedTree(8, 2, 2);

			Psi_ = TensorTreecd(rng_, tree_);

			// Generate an bit-flip operator and Fmatrix
			Matrixcd X(2, 2);
			X(0, 0) = 0.5;
			X(1, 1) = 0.5;
			LeafMatrixcd x(X);
			M_ = MLOcd(x, 0);
			M_.push_back(x, 3);
		}
	};

	TEST (TreeMarker) {
		Tree tree = TreeFactory::balancedTree(7, 4, 2);
		vector<size_t> modes({3, 4});
		SparseTree active(modes, tree);
			CHECK_EQUAL(7, active.size());
	}

	TEST_FIXTURE (HelperFactory, TreeMarker_NoTail) {
		/// Create TreeMarker omitting higher nodes in the tree after last branch
		SparseTree active(M_.targetLeaves(), tree_, false);
			CHECK_EQUAL(5, active.size());
	}

	TEST_FIXTURE (HelperFactory, Represent) {
		SparseMatrixTreecd mat = TreeFunctions::represent(M_, Psi_, tree_);
		const SparseTree& stree = mat.sparseTree();
		Matrixcd x(2, 2);
		x(0, 0) = x(1, 1) = 0.5;
			CHECK_CLOSE(0., residual(x, mat[stree.node(0)]), eps);
			CHECK_CLOSE(0., residual(x, mat[stree.node(1)]), eps);
			CHECK_CLOSE(0., residual(x, mat[stree.node(2)]), eps);
			CHECK_CLOSE(0., residual(x, mat[stree.node(3)]), eps);
		x(0, 0) = x(1, 1) = 0.25;
			CHECK_CLOSE(0., residual(x, mat[stree.node(4)]), eps);
		Matrixcd x1(1, 1);
		x1(0,0) = 0.;
			CHECK_CLOSE(0.0, residual(x1, mat[tree_.topNode()]), eps);
	}

	TEST_FIXTURE (HelperFactory, contraction) {
		SparseMatrixTreecd mats = TreeFunctions::represent(M_, Psi_, tree_);
		SparseMatrixTreecd holes(M_, tree_);
		TreeFunctions::contraction(holes, Psi_, Psi_, mats, tree_);
		const SparseTree& stree = holes.sparseTree();

		Matrixcd x(2, 2);
		x(0, 0) = 0.5;
			CHECK_CLOSE(0., residual(x, holes[stree.node(0)]), eps);
			CHECK_CLOSE(0., residual(x, holes[stree.node(1)]), eps);
			CHECK_CLOSE(0., residual(x, holes[stree.node(2)]), eps);
			CHECK_CLOSE(0., residual(x, holes[stree.node(3)]), eps);
		x(0, 0) = 1.;
			CHECK_CLOSE(0., residual(x, holes[stree.node(4)]), eps);
		Matrixcd x1(1, 1);
		x1(0, 0) = 1.;
			CHECK_CLOSE(0., residual(x1, holes[stree.node(5)]), eps);
	}

	TEST_FIXTURE (HelperFactory, Constructor) {
		SparseMatrixTreecd hmat(M_, tree_);
			CHECK_EQUAL(6, hmat.size());
	}

	TEST_FIXTURE (HelperFactory, IO) {
		SparseMatrixTreecd hmat(M_, tree_);
		hmat.write("SparseMatrixTree.dat");
		SparseMatrixTreecd gmat(M_, tree_);
		gmat.read("SparseMatrixTree.dat");
		const SparseTree& marker = gmat.sparseTree();
			CHECK_EQUAL(hmat.size(), gmat.size());
		for (size_t i = 0; i < marker.size(); ++i) {
			const Node& nodep = marker.node(i);
			double r = residual(hmat[nodep], gmat[nodep]);
				CHECK_CLOSE(0., r, eps);
		}
	}

	TEST_FIXTURE (HelperFactory, InverseTree) {
		SparseTree stree(M_, tree_);
		SparseTree itree(M_, tree_, true, true);
			CHECK_EQUAL(tree_.nNodes(), itree.size() + stree.size());
		for (const Node& node : tree_) {
				CHECK_EQUAL(true, stree.isActive(node) != itree.isActive(node));
				CHECK_EQUAL(true, stree.isActive(node) || itree.isActive(node));
		}
	}

	TEST_FIXTURE (HelperFactory, SparseTensorTree) {
		auto stree = make_shared<SparseTree>(M_.targetLeaves(), tree_, false);
		SparseTensorTreecd work(stree, tree_);
			CHECK_EQUAL(5, work.size());
		for (const Node *node_ptr : *stree) {
				CHECK_EQUAL(true, (node_ptr->shape() == work[*node_ptr].shape()));
		}
	}

	TEST_FIXTURE (HelperFactory, SparseTensorTree2) {
		SparseTree stree(M_.targetLeaves(), tree_, false);
		SparseTensorTreecd work(stree, tree_);
			CHECK_EQUAL(5, work.size());
		for (const Node *node_ptr : stree) {
				CHECK_EQUAL(true, (node_ptr->shape() == work[*node_ptr].shape()));
		}
	}
}

