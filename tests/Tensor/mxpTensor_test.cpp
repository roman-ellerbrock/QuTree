#include <gtest/gtest.h>
#include "Tensor/mxpTensor.h"
#include <chrono>
#include "Tensor/cuTensor.h"
#include "Tensor/TensorBLAS2.h"


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
    Tensord A = aranged({10, 10});
    Tensord A2r = A * A;

    mxpTensord B(A);
    mxpTensord B2 = B * B;

    Tensord A2 = B2.convert();

    EXPECT_NEAR(0., residual(A2, A2r), 1e-12);
}

/*
TEST (mxpTensor, mxpgemmBenchmark) {
    /// perform A^2 in hp and mxp
	using namespace chrono;
    time_point<steady_clock> start, end;

    size_t n = 1;

    for (size_t dim = 1000; dim <= 15000; dim+=3000) {
    {
        Tensord A = aranged({dim, dim});
        Tensord B(A.shape_);
        start = steady_clock::now();
        for (size_t i = 0; i < n; ++i) {
            gemm(B, A, A);
        }
        end = steady_clock::now();
    }
    double ms_hp = duration_cast<nanoseconds>(end - start).count() / ((double) n * 1000000);

    {
        mxpTensord mA(aranged({dim, dim}));
        mxpTensord mB(aranged({dim, dim}));
        start = steady_clock::now();
        for (size_t i = 0; i < n; ++i) {
            gemm(mB, mA, mA);
        }
        end = steady_clock::now();
    }
    double ms_mxp = duration_cast<nanoseconds>(end - start).count() / ((double) n * 1000000);
    cout << dim << " " << ms_hp << " " << ms_mxp << " \n";
    }
    getchar();

}*/

template <class Tensor, class ...Queue>
void matm(Tensor& c, const Tensor& a, const Tensor& b, Queue& ... queue) {
	gemm(c, a, b, 1., 0., blas::Op::NoTrans, blas::Op::NoTrans, queue...);
}

TEST(Tensor, queueforward) {
    Tensord A = aranged({100, 100});
    Tensord B = A;
    Tensord C(A.shape_);
    Tensord C2 = A * B;
    matm(C, B, A);
    EXPECT_NEAR(0., residual(C, C2), 1e-12);

    using namespace polymorphic;
	cuTensord cuA = transfer<double, cuMemory, hostMemory>(A);
    cuTensord cuB = cuA;
    cuTensord cuC(cuA.shape_);
	int device = 0;	
	int batch_size = 1000;
	blas::Queue queue(device, batch_size);
	blas::set_device(device);
//	gemm(cuC, cuA, cuB, 1., 0., blas::Op::NoTrans, blas::Op::NoTrans, queue);
    matm(cuC, cuB, cuA, queue);
	C = transfer<double, hostMemory, cuMemory>(cuC);
    EXPECT_NEAR(0., residual(C, C2), 1e-12);
}
