#include "Tensor/cuTensorBLAS1.h"
#include "Tensor/TensorBLAS2.h"
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

TEST_F(cuTensorBLAS1Fix, nrm2) {
    size_t inc = 1;
    double res = nrm2(A_, inc, queue_);
    double ref = nrm2(hA_, inc);
    EXPECT_NEAR(ref, res, 1e-12);
}

TEST_F(cuTensorBLAS1Fix, axpy) {
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
    Tensord res = transferFromGPUd(A2);
    EXPECT_NEAR(0., residual(res, hA2), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, substractEqual) {
    Tensord hA2 = aranged(A_.shape_);
    cuTensord A2 = transferToGPUd(hA2);
    A2 -= A_;
    hA2 -= hA_;
    Tensord res = transferFromGPUd(A2);
    EXPECT_NEAR(0., residual(res, hA2), 1e-12);
}

TEST_F(cuTensorBLAS1Fix, residual) {
    Tensord hA2 = deltad(A_.shape_);
    cuTensord A2 = transferToGPUd(hA2);
    double res = residual(A_, A2);
    double ref = residual(hA_, hA2);
    
    EXPECT_NEAR(ref, res, 1e-12);
}

TEST_F(cuTensorBLAS1Fix, diagonal) {
    cuTensord diag = diagonal(A_, queue_);
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

TEST_F(cuTensorBLAS1Fix, diagmPlusmdiag) {
    Tensord hdiag = aranged({A_.shape_[0]});
    cuTensord diag = transferToGPUd(hdiag);

    mdiagm(hRef_, hA_, hdiag, 1.);
    diagmm(hRef_, hdiag, hA_, 1.);

    diagmPlusmdiag(B_, A_, diag, 1.);
    Tensord res = transferFromGPUd(B_);
    EXPECT_NEAR(0., residual(res, hRef_), 1e-12);
}

TEST (cuTensor, dgemm) {
	cuTensord A({100, 100});
	cuTensord B({100, 100});
	cuTensord C({100, 100});

	int device = 0;	
	int batch_size = 1000;
	blas::Queue queue(device, batch_size);
	blas::set_device(device);
	gemm(C, A, B, 1., 0., blas::Op::NoTrans, blas::Op::NoTrans, queue);
}

TEST_F (cuTensorBLAS1Fix, operatorTimes) {
	B_ = A_;
	Tensord hB = hA_;
    Tensord hC = hA_ * hB;

    cuTensord C = A_ * B_;

    Tensord res = transferFromGPUd(C);
    EXPECT_NEAR(0., residual(hC, res), 1e-13);
}
