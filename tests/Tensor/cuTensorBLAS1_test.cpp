#include "Tensor/cuTensorBLAS1.h"
#include <gtest/gtest.h>

#define n 3

class cuTensorBLAS1Fix : public ::testing::Test {
    public:
    cuTensorBLAS1Fix() : B_({n, n}), hRef_({n, n}), queue_(0, 1000) {
        hA_ = aranged({n, n});
        A_ = transferToGPUd(hA_);
    }

    cuTensord A_;
    cuTensord B_;
    Tensord hA_;
    Tensord hRef_;
    blas::Queue queue_;
};

TEST_F(cuTensorBLAS1Fix, axpy) {
    //template void axpy<f, Tensor<f>>(const Tensor<f>& A, Tensor<f>& B, f alpha, size_t inc_a, size_t inc_b);
    double alpha = 1.;
    size_t inc = 1;
    axpy(A_, B_, alpha, inc, inc, queue_);
    axpy(hA_, hRef_, alpha, inc, inc);
    auto hB = transferFromGPUd(B_);
    EXPECT_NEAR(0., residual(hB, hRef_), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, castdf) {
    /// Test cast
    cuTensorf Af(A_.shape_);
    cast(Af, A_);
    Tensorf hRes = transferFromGPUf(Af);

    Tensorf hRef(hA_.shape_);
    cast(hRef, hA_);
    EXPECT_NEAR(0., residual(hRes, hRef), 1e-6);
}

TEST_F(cuTensorBLAS1Fix, addEqual) {
    Tensord hA2 = aranged(A_.shape_);
    cuTensord A2 = transferToGPUd(hA2);
    A2 += A_;
    hA2 += hA_;
    cuTensord res = transferToGPUd(hA2);
    EXPECT_NEAR(0., residual(hRes, hRef), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, diagonal) {
    cuTensord diag = diagonal(A_, queue_);
    queue_.sync();
    Tensord ref = diagonal(hA_);
    Tensord res = transferFromGPUd(diag);
    EXPECT_NEAR(0., residual(ref, res), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, addDiagonal) {
    cuTensord diag = diagonal(A_, queue_);
    addDiagonal(A_, diag, 1., queue_);
    Tensord hdiag = diagonal(hA_);
    addDiagonal(hA_, hdiag);
    Tensord res = transferFromGPUd(A_);
    EXPECT_NEAR(0., residual(res, hA_), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, offDiagonal) {
    cuTensord off(A_.shape_);
    offDiagonal(off, A_, queue_);
    Tensord ref(hA_.shape_);
    offDiagonal(ref, hA_);
    Tensord res = transferFromGPUd(off);
    EXPECT_NEAR(0., residual(res, ref), 1e-14);
}

TEST_F(cuTensorBLAS1Fix, diagmm) {
    Tensord hdiag = aranged({A_.shape_[0]});
    cuTensord diag = transferToGPUd(hdiag);
    Tensord hdiag2 = transferFromGPUd(diag);

    diagmm(hRef_, hdiag, hA_, 1.);
    diagmm(B_, diag, A_, 1.);

    Tensord res = transferFromGPUd(B_);
    EXPECT_NEAR(0., residual(res, hRef_), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, mdiagm) {
    Tensord hdiag = aranged({A_.shape_[0]});
    cuTensord diag = transferToGPUd(hdiag);

    mdiagm(hRef_, hA_, hdiag, 1.);
    mdiagm(B_, A_, diag, 1.);

    Tensord res = transferFromGPUd(B_);
    EXPECT_NEAR(0., residual(res, hRef_), 1e-12);
}

