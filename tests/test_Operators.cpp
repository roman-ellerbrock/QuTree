//
// Created by Roman Ellerbrock on 2020-01-24.
//
#include <gtest/gtest.h>
#include "TreeOperators/LeafOperator.h"
#include "TreeOperators/LeafMatrix.h"
#include "TreeOperators/MultiLeafOperator.h"
#include "TreeShape/LeafTypes/HO_Basis.h"
#include "TreeShape/Tree.h"
#include "TreeShape/TreeFactory.h"
#include "TreeOperators/SumOfProductsOperator_Implementation.h"

class HelperFactory : public ::testing::Test {
protected:
    HelperFactory() {
        Initialize();
    }

    ~HelperFactory() = default;
    Matrixcd X;
    LeafMatrixcd x;
    HO_Basis ho;

    void Initialize() {
        // Generate an bit-flip operator and Fmatrix
        X = Matrixcd(2, 2);
        X(0, 1) = 1.;
        X(1, 0) = 1.;
        x = LeafMatrixcd(X);
        ho = HO_Basis(10);
    }
};

TEST_F (HelperFactory, SPO_HO) {
    HO_Basis ho(10);
    ho.initialize(1., 0., 0., 1.);
    TensorShape tdim(vector<size_t>({10, 1}));
    Tensorcd A(tdim);
    ho.initSPF(A);
    Tensorcd xA(A.shape());
    ho.applyX(xA, A);
    string file("HO_Applyx.dat");
    xA.write(file);
    Tensorcd B(file);
    ASSERT_EQ(xA, B);
}

TEST_F (HelperFactory, SPO_SPOM) {
    // Check that applying a Matrix to a tensor
    // or the corresponding SPOM does the same

    TensorShape tdim(vector<size_t>({2, 1}));
    Tensorcd A(tdim);
    A(0) = 1.;
    Tensorcd XA = matrixTensor(X, A, 0);
    Tensorcd xA(tdim);
    x.apply(ho, xA, A);
    ASSERT_EQ(XA, xA);
}

TEST_F (HelperFactory, MPO_1) {
    MLOcd M(x, 1);
    mt19937 gen(time(nullptr));
    Tree tree = TreeFactory::balancedTree(4, 2, 2);
    TensorTreecd Chi(tree);
    Chi.fillRandom(gen, tree, false);
    auto Psi = M.apply(Chi, tree);

    /// Checking
/*		string filename("MLO.apply.dat");
		Psi.write(filename, false);
		TensorTreecd Xi(filename);
			ASSERT_EQ(Xi.size(), Psi.size());
		for (const Node& node : tree) {
				ASSERT_EQ(Xi[node], Psi[node]);
		}*/
}

TEST_F (HelperFactory, SOP_1) {
    MLOcd M(x, 1);
    MLOcd M2(x, 2);
    SOPcd S(M, 1.);
    S.push_back(M2, 0.5);
    auto M_sq = M * M;
    ASSERT_EQ(2, M_sq.size());
    complex<double> c(0.5);
    SOPcd cS = c * S;
    ASSERT_EQ(2, cS.size());
    SOPcd SS = S * S;
    ASSERT_EQ(4, SS.size());
}
