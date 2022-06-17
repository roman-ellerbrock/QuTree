//
// Created by Hannes Hoppe on 16.06.22.
//

#include <gtest/gtest.h>
#include "TreeClasses/SymTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Extension.h"
#include "TreeClasses/SpectralDecompositionTree.h"

class TTFactory : public ::testing::Test {
protected:
    TTFactory() {
        tree_ = TreeFactory::balancedTree(10, 4, 3);
        mt19937 gen(34676949);
        psi_ = SymTensorTree(gen, tree_, false);
        chi_ = SymTensorTree(gen, tree_, true);

        /// Operator initialization
        auto I = &LeafInterface::identity;
        for (size_t l = 0; l < tree_.nLeaves(); ++l) { I_.push_back(I, l); }
        stree_ = make_shared<SparseTree>(SparseTree(I_, tree_, false));
        SparseMatrixTreecd x1(stree_, tree_);
        SparseMatrixTreecd x2(stree_, tree_);
    }

    ~TTFactory() = default;

    Tree tree_;
    SymTensorTree psi_;
    SymTensorTree chi_;

    MLOcd I_;
    shared_ptr<SparseTree> stree_;
};

// will the swipe do the correct thong among the nodes?
TEST(Node, nextSCFNodeSwipe){
    Tree tree = TreeFactory::balancedTree(6, 2, 3);
    tree.info(std::cout);
}
