//
// Created by Roman Ellerbrock on 5/23/20.
//

#include "TreeOperators/TensorOperators/TensorOperatorTree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeOperators/TensorOperators/MatrixListTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/TensorTreeFunctions.h"
#include "TreeOperators/TensorOperators/TensorOperatorTreeFunctions.h"
#include <UnitTest++/UnitTest++.h>
#include "TreeOperators/TensorOperators/TTNOMatrixTree.h"
#include "TreeOperators/TensorOperators/TTNOHoleTree.h"
#include "TreeShape/LeafTypes/SpinGroup.h"
#include "TreeOperators/TensorOperators/contractSOP.h"

SUITE (TensorOperatorTree) {

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
		TensorOperatorTree H_;
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

			H_ = TensorOperatorTree(optree_);
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

	TEST(TTNO) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		Tree optree = TreeFactory::operatorTree(tree);

		mt19937 gen(2348);
		TensorOperatorTree A(optree);
		A.occupy(optree, gen);
		for (const Node& node : optree) {
			const Tensord& B = A[node];
			for (size_t i = 0; i < node.shape().totalDimension(); ++i) {
				double r = abs(B[i]);
				CHECK_EQUAL((r > 1e-15), true);
			}
		}

	}

	TEST(TTNOrep) {
		SOPd S;
		Tree tree = TreeFactory::balancedTree(16, 2, 3);
		Tree optree = TreeFactory::operatorTree(tree);

		optree.print();

		for (size_t l = 0; l < tree.nLeaves(); l++) {
			Matrixd sigma = JordanWigner::sigmaX();
			MLOd M(sigma, l);
			S.push_back(M, 1.);
		}
		cout << "S.size(): " << S.size() << endl;

		mt19937 gen(time(NULL));
		TensorOperatorTree A(optree, gen);

		orthogonal(A, optree);
		orthonormal(A, optree);

		cout << "rep: " << endl;
		TTNOMatrixTree rep(S, optree);
		rep.represent(A, S, optree);
		rep.print(optree);

		cout << "holes: " << endl;
		TTNOHoleTree hole(S, optree);
		hole.represent(A, rep, optree);
		hole.print(optree);

		TensorOperatorTree B = contractSOP(A, S, optree);
		cout << "Final result:\n";
		B.print(optree);
	}

}
