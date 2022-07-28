#include "Tensor/cuTensorBLAS1.h"
#include <gtest/gtest.h>


class cuTensorBLAS1Fix : public ::testing::Test {
    public:
    cuTensorBLAS1Fix() : B_({100, 100}), hBr_({100, 100}), queue_(0, 1000) {
        hA_ = aranged({100, 100});
        A_ = transferToGPUd(hA_);
    }

    cuTensord A_;
    cuTensord B_;
    Tensord hA_;
    Tensord hBr_;
    blas::Queue queue_;
};

/*TEST_F(cuTensorBLAS1Fix, nrm2) {
    double a = nrm2(A_, 1, queue_);
    double ha = nrm2(hA_, 1);
    EXPECT_NEAR(ha, a, 1e-12);
}*/

TEST_F(cuTensorBLAS1Fix, axpy) {
    //template void axpy<f, Tensor<f>>(const Tensor<f>& A, Tensor<f>& B, f alpha, size_t inc_a, size_t inc_b);
    double alpha = 1.;
    axpy(A_, B_, alpha, 1, 1, queue_);
    axpy(hA_, hBr_, alpha, 1, 1);
    auto hB = transferFromGPUd(B_);
    EXPECT_NEAR(0., residual(hB, hBr_), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, castdf) {
    //template void axpy<f, Tensor<f>>(const Tensor<f>& A, Tensor<f>& B, f alpha, size_t inc_a, size_t inc_b);
    cuTensorf Af(A_.shape_);
    cast(Af, A_);
    Tensorf hAf = transferFromGPUf(Af);

    Tensorf hAf2(A_.shape_);
    cast(hAf2, hA_);

    EXPECT_NEAR(0., residual(hAf, hAf2), 1e-6);
}

