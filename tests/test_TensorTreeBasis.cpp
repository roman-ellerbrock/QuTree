//
// Created by Roman Ellerbrock on 2/3/20.
//
#include <gtest/gtest.h>
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"

TEST (TensorTreeBasis, TensorTreeBasis_Generator) {
    size_t n_leaf = 4;
    size_t n_node = 2;
    TensorShape tdim_top({n_node, n_node, 1});
    TensorShape tdim_upper({n_node, n_node, n_node});
    TensorShape tdim_bottom({n_leaf, n_node});

    for (size_t n_modes = 2; n_modes < 18; ++n_modes) {
        Tree tree = TreeFactory::balancedTree(n_modes, n_leaf, n_node);

        for (const Node& node : tree) {
            const TensorShape& tdim = node.shape();
            if (node.isToplayer()) {
                ASSERT_EQ(tdim_top, tdim);
            } else if (node.isBottomlayer()) {
                ASSERT_EQ(tdim_bottom, tdim);
            } else {
                ASSERT_EQ(tdim_upper, tdim);
            }
        }
    }
}

TEST (TensorTreeBasis, TensorTreeBasis_FileIO) {
    size_t n_leaf = 4;
    size_t n_node = 2;
    TensorShape tdim_top({n_node, n_node, 1});
    TensorShape tdim_upper({n_node, n_node, n_node});
    TensorShape tdim_bottom({n_leaf, n_node});
    size_t n_modes = 13;

    Tree tree = TreeFactory::balancedTree(n_modes, n_leaf, n_node);
    {
        ofstream os("TTBasis.IO.tmp.dat");
        tree.write(os);
        os.close();
    }
    Tree tree2("TTBasis.IO.tmp.dat");
    ASSERT_EQ(tree2.nNodes(), tree.nNodes());
}

TEST (TensorTreeBasis, TensorTreeBasis_Reindexing) {
    size_t n_modes = 9;
    Tree tree = TreeFactory::balancedTree(n_modes, 2, 4);

    map<size_t, size_t> Map;
    for (size_t k = 0; k < n_modes; ++k) {
        Map[k] = n_modes - 1 - k;
    }
    tree.reindexLeafModes(Map);

    size_t k = 0;
    for (const Node& node : tree) {
        if (node.isBottomlayer()) {
            const Leaf& leaf = node.getLeaf();
            ASSERT_EQ(k++, leaf.mode());
        }
    }
}

TEST (TensorTreeBasis, TensorTreeBasis_Train) {
    size_t nLeaves = 12;
    auto tree = TreeFactory::unbalancedTree(nLeaves, 4, 2, 6);
    ASSERT_EQ(2 * nLeaves - 1, tree.nNodes());
}

TEST (TensorTreeBasis, TensorTreeBasis_Copy) {
    /// Construct a tree and check that it works
    Tree tree = TreeFactory::balancedTree(12, 4, 3);
    ASSERT_EQ(true, tree.isWorking());

    {
        /// Copy-constructor test
        Tree tree_copy_construct(tree);
        ASSERT_EQ(true, tree_copy_construct.isWorking());

        /// Move constructor
        Tree tree_move_construct(move(tree_copy_construct));
        ASSERT_EQ(true, tree_move_construct.isWorking());
    }

    {
        /// Copy-asignment test
        Tree tree_copy_asign = tree;
        ASSERT_EQ(true, tree_copy_asign.isWorking());

        /// Move asignment operator
        Tree tree_move_asign = move(tree_copy_asign);
        ASSERT_EQ(true, tree_move_asign.isWorking());
    }
}