//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/TensorTree.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/TensorTreeFunctions.h"

SUITE (TensorTree) {

	TEST (TensorTree_FILE_IO) {
		Tree tree = TreeFactory::balancedTree(12, 2, 4);
		TensorTreecd T(tree);
		T.write("TT.tmp.dat", true);
		TensorTreecd Q(tree);
		ifstream is("TT.tmp.dat");
		is >> Q;
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : tree) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_RandomGenerate) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		TensorTreecd T(tree);
		mt19937 gen(2468);
		T.fillRandom(gen, tree, false);
		string filename("TT.RNG.dat");
		T.write(filename, true);
		TensorTreecd Q(filename);
			CHECK_EQUAL(T.size(), Q.size());
		for (const Node& node : tree) {
				CHECK_EQUAL(T[node], Q[node]);
		}
	}

	TEST (TensorTree_Train) {
		size_t nLeaves = 12;
		auto tree = TreeFactory::unbalancedTree(nLeaves, 4, 2, 6);
			CHECK_EQUAL(2 * nLeaves - 1, tree.nNodes());
		mt19937 gen(2468);
		TensorTreecd Psi(gen, tree);
			CHECK_EQUAL(tree.nNodes(), Psi.size());
	}

	TEST (TreeTest) {
		Tree tree = TreeFactory::balancedTree(12, 5, 2);
		TensorTreecd Psi(tree);
		for (const Tensorcd& A : Psi) {
			TensorShape shape = A.shape();
		}
		for (const Node& node : tree) {
			Tensorcd& A = Psi[node];
		}
	}

	TEST (DirectSum) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		mt19937 gen(234);
		TensorTreecd Psi(gen, tree);
		TensorTreecd Chi(gen, tree);
		TreeFunctions::sum(Psi, tree, Chi, true, true);

		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			if (node.isToplayer()) {
					CHECK_EQUAL(3, shape.order());
					CHECK_EQUAL(1, shape.lastDimension());
					CHECK_EQUAL(4 * 4, shape.lastBefore());
			} else if (node.isBottomlayer()) {
					CHECK_EQUAL(2, shape.order());
					CHECK_EQUAL(4, shape.lastDimension());
					CHECK_EQUAL(2, shape.lastBefore());
			} else {
					CHECK_EQUAL(3, shape.order());
					CHECK_EQUAL(4, shape.lastDimension());
					CHECK_EQUAL(4 * 4, shape.lastBefore());
			}
		}
	}

	/// Calculate grid size in CDVR and sCDVR
/*	TEST (GridSize) {
		for (size_t nleaves = 2; nleaves <= pow(2, 11); nleaves*=2){
			Tree tree = TreeFactory::balancedTree(nleaves, 10, 10);

			/// sCDVR
			size_t sgrid = 0;
			for (const Node& node: tree) {
				const TensorShape& shape = node.shape();
				size_t loc = shape.totalDimension();
				if (!node.isToplayer()) {
					loc += shape.lastDimension()*shape.lastDimension();
				}
				sgrid += loc;
			}
			/// CDVR
			size_t grid = 0;
			for (Node& node: tree) {
				const TensorShape& shape = node.shape();
				size_t loc = shape.lastBefore();
				if (node.isToplayer()) {
					grid += loc;
					break;
				}
				Node* p = &(node.parent());
				while(true) {
					const TensorShape& pshape = p->shape();
					loc *= pshape.lastBefore() / pshape[p->childIdx()];
					if (p->isToplayer()) { break; }
					p = &(p->parent());
				}
				grid += loc;
			}
			cout << nleaves << " " << sgrid << " " << grid << endl;
		}
		exit(0);
	}
 */
}

