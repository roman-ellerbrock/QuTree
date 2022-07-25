#include <gtest/gtest.h>
#include "Tensor/mxpTensor.h"


TEST (mxpTensor, Constructor) {
    TensorShape shape({100, 100});
    mxpTensord A(shape);
    EXPECT_EQ(A.diag_.shape_[0], 100);

    EXPECT_EQ(A.off_.shape_[0], 100);
    EXPECT_EQ(A.off_.shape_[1], 100);

    EXPECT_EQ(A.offLP_.shape_[0], 100);
    EXPECT_EQ(A.offLP_.shape_[1], 100);
}

TEST (mxpTensor, ConstructFrom) {
    /// Convert Tensor into mixed prec. Tensor and back
    Tensord A = aranged({100, 100});

    mxpTensord mA(A);
    Tensord B = mA.convert();

    EXPECT_NEAR(0., residual(A, B), 1e-12);
}

TEST (mxpTensor, mdiagm) {
    size_t n = 10;
    Tensord A = aranged({n, n});
    Tensord B({n, n});
    Tensord diag = aranged({n});
    for (size_t i = 0; i < n; ++i) {
        B(i, i) = diag(i);
    }

    Tensord C({n, n});
    mdiagm(C, A, diag);

    auto C2 = A * B;
    EXPECT_NEAR(0., residual(C, C2), 1e-12);
}

TEST (mxpTensor, diagmm) {
    size_t n = 10;
    Tensord A = aranged({n, n});
    Tensord B({n, n});
    Tensord diag = aranged({n});
    for (size_t i = 0; i < n; ++i) {
        B(i, i) = diag(i);
    }

    Tensord C({n, n});
    diagmm(C, diag, A);

    auto C2 = B * A;
    EXPECT_NEAR(0., residual(C, C2), 1e-12);
}

TEST (mxpTensor, gemm) {
    /// perform A^2 in hp and mxp
    Tensord A = aranged({3, 3});
    Tensord A2r = A * A;

    mxpTensord B(A);
    mxpTensord C(A);
    mxpTensord B2 = gemm(B, C);

    Tensord A2 = B2.convert();
    cout << "ref:\n";
    A2r.print();
    cout << "A2:\n";
    A2.print();
    EXPECT_NEAR(0., residual(A2, A2r), 1e-12);
    getchar();
}
