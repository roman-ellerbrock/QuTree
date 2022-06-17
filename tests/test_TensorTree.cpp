//
// Created by Roman Ellerbrock on 2020-01-19.
//
#include <gtest/gtest.h>
#include "TreeClasses/TensorTree.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/TensorTreeFunctions.h"


TEST (TensorTree, TensorTree_FILE_IO) {
    Tree tree = TreeFactory::balancedTree(12, 2, 4);
    TensorTreecd T(tree);
    T.write("TT.tmp.dat", true);
    TensorTreecd Q(tree);
    ifstream is("TT.tmp.dat");
    is >> Q;
    ASSERT_EQ(T.size(), Q.size());
    for (const Node& node : tree) {
        ASSERT_EQ(T[node], Q[node]);
    }
}

TEST (TensorTree, TensorTree_RandomGenerate) {
    Tree tree = TreeFactory::balancedTree(12, 2, 2);
    TensorTreecd T(tree);
    mt19937 gen(2468);
    T.fillRandom(gen, tree, false);
    string filename("TT.RNG.dat");
    T.write(filename, true);
    TensorTreecd Q(filename);
    ASSERT_EQ(T.size(), Q.size());
    for (const Node& node : tree) {
        ASSERT_EQ(T[node], Q[node]);
    }
}

TEST (TensorTree, TensorTree_Train) {
    size_t nLeaves = 12;
    auto tree = TreeFactory::unbalancedTree(nLeaves, 4, 2, 6);
    ASSERT_EQ(2 * nLeaves - 1, tree.nNodes());
    mt19937 gen(2468);
    TensorTreecd Psi(gen, tree);
    ASSERT_EQ(tree.nNodes(), Psi.size());
}

TEST (TensorTree, TreeTest) {
    Tree tree = TreeFactory::balancedTree(12, 5, 2);
    TensorTreecd Psi(tree);
    for (const Tensorcd& A : Psi) {
        TensorShape shape = A.shape();
    }
    for (const Node& node : tree) {
        Tensorcd& A = Psi[node];
    }
}

TEST (TensorTree, DirectSum) {
    Tree tree = TreeFactory::balancedTree(12, 2, 2);
    mt19937 gen(234);
    TensorTreecd Psi(gen, tree);
    TensorTreecd Chi(gen, tree);
    TreeFunctions::sum(Psi, tree, Chi, true, true);

    for (const Node& node : tree) {
        const TensorShape& shape = node.shape();
        if (node.isToplayer()) {
            ASSERT_EQ(3, shape.order());
            ASSERT_EQ(1, shape.lastDimension());
            ASSERT_EQ(4 * 4, shape.lastBefore());
        } else if (node.isBottomlayer()) {
            ASSERT_EQ(2, shape.order());
            ASSERT_EQ(4, shape.lastDimension());
            ASSERT_EQ(2, shape.lastBefore());
        } else {
            ASSERT_EQ(3, shape.order());
            ASSERT_EQ(4, shape.lastDimension());
            ASSERT_EQ(4 * 4, shape.lastBefore());
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

