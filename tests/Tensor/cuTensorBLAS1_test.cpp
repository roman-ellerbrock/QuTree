#include "Tensor/cuTensorBLAS1.h"
#include <gtest/gtest.h>


class cuTensorBLAS1Fix : public ::testing::Test {
    public:
    cuTensorBLAS1Fix() : B_({100, 100}), hRef_({100, 100}), queue_(0, 1000) {
        hA_ = aranged({100, 100});
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
    axpy(A_, B_, alpha, 1, 1, queue_);
    axpy(hA_, hRef_, alpha, 1, 1);
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

TEST_F(cuTensorBLAS1Fix, diagmm) {
    Tensord hdiag = aranged({A_.shape_[0]});
    cuTensord diag = transferToGPUd(hdiag);

    diagmm(hRef_, hdiag, hA_, 1.);
    diagmm(B_, diag, A_, 1.);

    Tensord res = transferFromGPUd(B_);
    res.print();getchar();
    EXPECT_NEAR(0., residual(res, hRef_), 1e-12);
}

