//
// Created by Roman Ellerbrock on 5/23/20.
//

#include <UnitTest++/UnitTest++.h>
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeShape/LeafTypes/SpinGroup.h"
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "TreeOperators/TensorOperators/TTOrepresentation.h"
#include "TreeOperators/TensorOperators/TTOcontraction.h"

SUITE (TensorTreeOperator) {

	double eps = 1e-7;

	class HelperFactory {
	public:
		HelperFactory() {
			Initialize();
		}

		~HelperFactory() = default;
		mt19937 rng_;
		Tree tree_;
		Tree optree_;
		TensorTreed Psi_;
		TensorTreeOperator H_;
		LeafMatrixd leafI_;
		LeafMatrixd leafX_;

		void Initialize() {
			rng_ = mt19937(1990);
			tree_ = TreeFactory::balancedTree(12, 2, 2);
			optree_ = TreeFactory::operatorTree(tree_);
			Psi_ = TensorTreed(rng_, tree_);

			Matrixd I = identityMatrixd(2);
			Matrixd X(2, 2);
			X(0, 1) = 1.;
			X(1, 0) = 1.;
			leafI_ = LeafMatrixd(I);
			leafX_ = LeafMatrixd(X);

			H_ = TensorTreeOperator(optree_);
			H_.occupy(optree_);
			for (const Node& node : optree_) {
				if (node.isBottomlayer()) {
					H_.setLeafOperator(leafI_, 0, node);
					H_.setLeafOperator(leafX_, 1, node);
				}
			}
		}
	};

/*	TEST (OperatorTree) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		Tree optree = TreeFactory::operatorTree(tree);
		for (const Node& node : optree) {
			if (node.isBottomlayer()) {
					CHECK_EQUAL(2, node.shape().order());
					CHECK_EQUAL(4, node.shape().lastDimension());
					CHECK_EQUAL(4, node.shape().lastBefore());
			}
		}
	}*/

	TEST (TTNO) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		Tree optree = TreeFactory::operatorTree(tree);

		mt19937 gen(2348);
		TensorTreeOperator A(optree);
		A.occupy(optree, gen);
		for (const Node& node : optree) {
			const Tensord& B = A[node];
			for (size_t i = 0; i < node.shape().totalDimension(); ++i) {
				double r = abs(B[i]);
					CHECK_EQUAL((r > 1e-15), true);
			}
		}
	}

	TEST (TTNOrep) {
		SOPd S;
		Tree tree = TreeFactory::balancedTree(32, 2, 3);
		Tree optree = TreeFactory::operatorTree(tree);

		for (size_t l = 0; l < tree.nLeaves(); l++) {
			Matrixd sigma = JordanWigner::sigmaX();
			MLOd M(sigma, l);
			S.push_back(M, 1.);
		}

		mt19937 gen(time(NULL));
		TensorTreeOperator A(optree, gen);

		orthogonal(A, optree);
		orthonormal(A, optree);

		double err0 = error(A, S, optree);
			CHECK_EQUAL(1, (err0 > 1e-1));
		TensorTreeOperator B = contractSOP(A, S, 2, optree, nullptr);
		double err = error(B, S, optree);
			CHECK_EQUAL(1, (err < 1e-12));
	}

	TEST (TTNOrep_nonhermitian) {
		SOPd S;
		Tree tree = TreeFactory::balancedTree(6, 2, 3);
		Tree optree = TreeFactory::operatorTree(tree, 3);

		for (size_t l = 0; l < tree.nLeaves(); l++) {
			Matrixd sigma = JordanWigner::sigmaPlus();
			if (l % 2) { sigma = JordanWigner::sigmaMinus(); }
			MLOd M(sigma, l);
			for (size_t i = 0; i < l; ++i) {
				if (l != i) { M.push_back(JordanWigner::sigmaZ(), i); }
			}
			S.push_back(M, 1./(double)(l+1.));
		}

		mt19937 gen(time(NULL));
		TensorTreeOperator A(optree, gen);

		orthogonal(A, optree);
		orthonormal(A, optree);

		double err0 = error(A, S, optree);
			CHECK_EQUAL(1, (err0 > 1e-1));
		TensorTreeOperator B = contractSOP(A, S, 30, optree, nullptr);
		double err = error(B, S, optree);
			CHECK_EQUAL(1, (err < 1e-12));

		TTOrepresentation rep(tree, optree);
		TensorTreed Psi(gen, tree);
		rep.calculate(Psi, B, Psi, optree);
			CHECK_CLOSE(0., rep[optree.topNode()][0][0], eps);

		TTOcontraction con(tree, optree);
		con.calculate(Psi, B, rep, Psi, optree);
	}
}
