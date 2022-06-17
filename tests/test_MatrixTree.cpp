//
// Created by Roman Ellerbrock on 2/12/20.
//

#include <gtest/gtest.h>
#include "TreeClasses/MatrixTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "Util/RandomMatrices.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/TreeTransformation.h"
#include "TreeClasses/TensorTreeFunctions.h"


double mteps = 1e-8;

TEST (MatrixTree, Constructor) {
    Tree tree = TreeFactory::balancedTree(7, 5, 4);
    MatrixTreecd M(tree);
    ASSERT_EQ(tree.nNodes(), M.size());
}

TEST (MatrixTree, IO) {
    Tree tree = TreeFactory::balancedTree(7, 5, 4);
    MatrixTreecd M(tree);
    mt19937 gen(1988);
    for (const Node& node : tree) {
        Matrixcd& m = M[node];
        m = RandomMatrices::gue(m.dim1(), gen);
    }

    string filename("MatrixTree.dat");
    M.write(filename);
    MatrixTreecd Mread(filename);
    for (const Node& node : tree) {
        double delta = residual(M[node], Mread[node]);
        ASSERT_NEAR(0., delta, mteps);
    }
}


using namespace TreeFunctions;

TEST (MatrixTreeFunctions, dotProduct) {
    mt19937 gen(1923);
    Tree tree = TreeFactory::balancedTree(7, 5, 4);
    TensorTreecd Psi(gen, tree);
    orthogonal(Psi, tree);
    MatrixTreecd S = dotProduct(Psi, Psi, tree);
    for (const Node& node : tree) {
        const Matrixcd& s = S[node];
        Matrixcd Identity = identityMatrix<complex<double>>(s.dim1());
        double r = residual(Identity, s);
        ASSERT_NEAR(0., r, mteps);
    }
}

TEST (MatrixTreeFunctions, contraction) {
    mt19937 gen(1923);
    Tree tree = TreeFactory::balancedTree(7, 5, 4);
    TensorTreecd Psi(gen, tree, false);
    TensorTreecd Chi(gen, tree, false);
    MatrixTreecd S = dotProduct(Psi, Chi, tree);
    MatrixTreecd Rho = contraction(Psi, Chi, S, tree);
    ASSERT_EQ(tree.nNodes(), Rho.size());
    ASSERT_EQ(tree.nNodes(), S.size());
}

TEST (MatrixTreeFunctions, Density) {
    mt19937 gen(1923);
    Tree tree = TreeFactory::balancedTree(7, 5, 4);
    TensorTreecd Psi(gen, tree, true);
    MatrixTreecd Rho = contraction(Psi, tree, true);
    for (const Node& node : tree) {
        if (!node.isToplayer()) {
            Matrixcd& rho = Rho[node];
            ASSERT_EQ(rho.dim2(), rho.dim1());
            for (size_t j = 0; j < rho.dim2(); ++j) {
                for (size_t i = 0; i < rho.dim1(); ++i) {
                    if (j == 0 && i == 0) {
                        ASSERT_NEAR(1., abs(rho(j, i)), mteps);
                    } else {
                        ASSERT_NEAR(0., abs(rho(j, i)), mteps);
                    }
                }
            }
        }
    }
}

TEST (MatrixTreeFunctions, SpectralDecompositionTree_Calc) {
    mt19937 gen(1993);
    Tree tree = TreeFactory::balancedTree(12, 2, 2);
    TensorTreecd Psi(gen, tree);
    MatrixTreecd Rho = TreeFunctions::contraction(Psi, tree, true);
    SpectralDecompositionTreecd X(Rho, tree);
    ASSERT_EQ(Rho.size(), X.size());
    for (const Node& node : tree) {
        if (!node.isToplayer()) {
            ASSERT_NEAR(1., X[node].second(1), mteps);
        }
    }
}

TEST (MatrixTreeFunctions, SpectralDecompositionTree_Inverse) {
    Tree tree = TreeFactory::balancedTree(12, 4, 2);
    mt19937 gen(1993);
    MatrixTreecd H(tree);
    for (const Node& node : tree) {
        const TensorShape& dim = node.shape();
        auto mat = RandomMatrices::gue(dim.lastDimension(), gen);
        auto mat_dagger = mat.adjoint();
        H[node] = mat * mat_dagger;
    }

    SpectralDecompositionTreecd X(H, tree);
    auto H_inv = X.invert(tree, 1e-10);

    MatrixTreecd Identity(tree);
    for (const Node& node : tree) {
        Identity[node] = H_inv[node] * H[node];
    }
    for (const Node& node : tree) {
        const Matrixcd& I_test = Identity[node];
        auto r = residual(I_test, identityMatrix<complex<double>>(I_test.dim1()));
        ASSERT_NEAR(0., r, mteps);
    }
}

TEST (MatrixTreeFunctions, canonicalTransformation) {

    Tree tree = TreeFactory::balancedTree(12, 4, 2);
    mt19937 gen(1993);
    TensorTreecd Psi(gen, tree);

    canonicalTransformation(Psi, tree, true);
    auto rho = TreeFunctions::contraction(Psi, tree, true);
    double off = 0.;
    for (const auto& mat : rho) {
        for (size_t j = 0; j < mat.dim2(); ++j) {
            for (size_t i = 0; i < mat.dim1(); ++i) {
                if (i != j) { off += abs(mat(i, j)); }
            }
        }
    }
    ASSERT_NEAR(0., off, mteps);
}

TEST (MatrixTreeFunctions, CanonicalTransformation2) {

    Tree tree = TreeFactory::balancedTree(12, 4, 2);
    mt19937 gen(1993);
    TensorTreecd Psi(gen, tree, false);
    auto rho = TreeFunctions::contraction(Psi, tree, true);
    canonicalTransformation(Psi, tree, true);
    rho = TreeFunctions::contraction(Psi, tree, true);

    auto S = TreeFunctions::dotProduct(Psi, Psi, tree);

    rho = TreeFunctions::contraction(Psi, tree, true);
    double off = 0.;
    for (const auto& mat : rho) {
        for (size_t j = 0; j < mat.dim2(); j++) {
            for (size_t i = 0; i < mat.dim1(); ++i) {
                if (i != j) { off += abs(mat(i, j)); }
            }
        }
    }
    ASSERT_NEAR(0., off, mteps);
}


TEST(TreeTransformations, ContractionNormalized) {
    Tree tree = TreeFactory::balancedTree(12, 4, 2);
    // Increase number of states
    Node& top = tree.topNode();
    TensorShape& shape = top.shape();
    shape[shape.lastIdx()] += 1;
    tree.update();
    mt19937 gen(1993);
    TensorTreecd Psi(gen, tree);

    auto edgmtepsi = TreeFunctions::contractionNormalization(Psi, tree, true);

    for (const Edge& e : tree.edges()) {
        auto phi = edgmtepsi[e];
        Matrixcd deltaij = contraction(phi, phi, e.upIdx());
        Matrixcd I = identityMatrix<complex<double>>(deltaij.dim1());
        auto r = residual(deltaij, I);
        ASSERT_NEAR(0., r, 1e-7);
    }
}

TEST(TreeTransformations, CompressTensorTree) {
    Tree tree = TreeFactory::balancedTree(12, 4, 2);
    mt19937 gen(1993);
    TensorTreecd Psi(gen, tree);
    canonicalTransformation(Psi, tree, true);
    MatrixTreecd Rho = TreeFunctions::contraction(Psi, tree, true);
    SpectralDecompositionTreecd X(Rho, tree);

    TreeFunctions::adjust(Psi, tree, X, 1e-7);
    for (const Node& node : tree) {
        const TensorShape& shape = node.shape();
        /// Check Node TensorShape
        if (!node.isBottomlayer()) {
            ASSERT_EQ(1, shape.totalDimension());
        } else {
            ASSERT_EQ(4, shape.totalDimension());
        }
        /// Check Tensor
        Tensorcd Phiacc(shape);
        Phiacc(0) = 1.;
        ASSERT_NEAR(0., residual(Psi[node], Phiacc), mteps);
    }
}