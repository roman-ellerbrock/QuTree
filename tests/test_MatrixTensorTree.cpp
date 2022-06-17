//
// Created by Roman Ellerbrock on 4/3/20.
//

#include <gtest/gtest.h>
#include "TreeClasses//MatrixTensorTree.h"
#include "TreeShape/TreeFactory.h"

double eps = 1e-7;

TEST (ExplicitEdgeWavefunction, Init) {
    Tree tree = TreeFactory::balancedTree(10, 2, 3);
    mt19937 gen(34676949);
    TensorTreecd Psi(gen, tree, false);

    /// transform to symmetric representation
    MatrixTensorTree Chi(Psi, tree, true);

    /// Check re-obtaining wavefunction
    TensorTreecd Psi2 = Chi.bottomUpNormalized(tree);
    auto S = TreeFunctions::dotProduct(Psi, Psi2, tree);
    for (const Node& node : tree) {
        const auto& s = S[node];
        ASSERT_NEAR(0., residual(s * s, identityMatrixcd(s.dim1())), eps);
    }

    /// Check top-down
    const MatrixTreecd rho = TreeFunctions::contraction(Psi, tree, true);
    const MatrixTreecd T = TreeFunctions::contraction(Psi2, Psi, S, tree);
    const TensorTreecd& Atilde = Chi.nodes();
    for (const Node& node : tree) {
        if (!node.isToplayer()) {
            const Node& parent = node.parent();
            auto x = contraction(Atilde[parent], Atilde[parent], node.childIdx());
            auto v1 = diagonalize(x).second;
            auto v2 = diagonalize(rho[node]).second;
            auto r = residual(v1, v2);
            ASSERT_NEAR(0., r, eps);
        }
    }
    /// @TODO: Add strong test for top-down.
    ASSERT_EQ(true, IsWorking(Chi, tree, eps));
}

