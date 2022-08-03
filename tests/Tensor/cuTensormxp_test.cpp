#include <gtest/gtest.h>
#include "Tensor/cuTensormxp.h"
#include "Tensor/cuTensorBLAS1.h"
#include <chrono>

TEST(cuTensormxpF, constructor) {
    TensorShape shape({100, 100});
    cuTensordf A(shape);
    EXPECT_EQ(A.diag_.shape_[0], 100);

    EXPECT_EQ(A.off_.shape_[0], 100);
    EXPECT_EQ(A.off_.shape_[1], 100);

    EXPECT_EQ(A.offLP_.shape_[0], 100);
    EXPECT_EQ(A.offLP_.shape_[1], 100);
}

TEST (cuTensormxpF, ConstructFrom) {
    /// Convert Tensor into mixed prec. Tensor and back
    Tensord hA = aranged({100, 100});
    cuTensord A = transferToGPUd(hA);

    blas::Queue queue(0, 1000);
    cuTensordf mA(A, queue);
    cuTensord B = mA.convert(queue);
    Tensord hB = transferFromGPUd(B);

    EXPECT_NEAR(0., residual(hA, hB), 1e-12);
}

TEST (cuTensormxpF, gemm) {
    /// perform A^2 in hp and mxp
    Tensord hA = aranged({10, 10});
    Tensord ref = hA * hA;

    cuTensord A = transferToGPUd(hA);
    blas::Queue& queue = qutree::queue;
    cuTensordf B(A, queue);
    cuTensordf B2(A.shape_);
    gemm(B2, B, B, queue);

    cuTensord A2 = B2.convert(qutree::queue);
    Tensord hA2 = transferFromGPUd(A2);

    EXPECT_NEAR(0., residual(hA2, ref), 1e-12);
}

TEST (cuTensormxpF, gemm_return) {

    Tensord hA = aranged({10, 10});
    Tensord ref = hA * hA;

    cuTensord A = transferToGPUd(hA);
    blas::Queue& queue = qutree::queue;
    cuTensordf B(A, queue);
    cuTensordf B2 = gemm(B, B, queue);

    cuTensord A2 = B2.convert(qutree::queue);
    Tensord hA2 = transferFromGPUd(A2);

    EXPECT_NEAR(0., residual(hA2, ref), 1e-12);
}

TEST (cuTensormxpF, operatormult) {
    Tensord hA = aranged({10, 10});
    cuTensord Ad = transferToGPUd(hA);
    cuTensordf A(Ad, qutree::queue);
    cuTensordf A2 = A * A;
    cuTensord res = A2.convert(qutree::queue);
    cuTensord ref = Ad * Ad;

    EXPECT_NEAR(0., residual(ref, res), 1e-12);
}

/*
TEST (cuTensormxpBenchmark, gemmBenchmark) {
    /// perform A^2 in hp and mxp
	using namespace chrono;
    time_point<steady_clock> start, end;

    size_t n = 3;

    for (size_t dim = 1000; dim <= 7500; dim+=500) {
    {
        Tensord hA = aranged({dim, dim});
        cuTensord A = transferToGPUd(hA);
        cuTensord B(A.shape_);
        start = steady_clock::now();
        for (size_t i = 0; i < n; ++i) {
            B = A * A;
        }
        end = steady_clock::now();
    }

    double ms_hp = duration_cast<nanoseconds>(end - start).count() / ((double) n * 1000000);

    {
        Tensord hA = aranged({dim, dim});
        cuTensord A = transferToGPUd(hA);
        cuTensordf mA(A, qutree::queue);
        cuTensordf mB(A.shape_);
        start = steady_clock::now();
        for (size_t i = 0; i < n; ++i) {
            mB = mA * mA;
        }
        end = steady_clock::now();
    }

    double ms_mxp = duration_cast<nanoseconds>(end - start).count() / ((double) n * 1000000);

    {
        Tensorf hA = arangef({dim, dim});
        cuTensorf A = transferToGPUf(hA);
        cuTensorf B(A.shape_);
        start = steady_clock::now();
        for (size_t i = 0; i < n; ++i) {
            B = A * A;
        }
        end = steady_clock::now();
    }
    double ms_lp = duration_cast<nanoseconds>(end - start).count() / ((double) n * 1000000);
    cout << dim << "\t" << ms_lp << "\t" << ms_mxp << "\t" << ms_hp << "\n";
    }
    getchar();

}
*/
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
        Tensord mA(aranged({dim, dim}));
        Tensord mB(aranged({dim, dim}));
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

}
*/

/*TEST (mxpTensor, gemmAccuracy) {
    size_t dim = 10;
    Tensord R = randomd({dim, dim});
    Tensord diag = diagonal(randomd({dim, dim}));
    for (size_t n = 0; n < 11; ++n) {
        double eps = 1e-10 * pow(10, n);
        Tensord A = eps * R;
        for (size_t i = 0; i < dim; ++i) {
            A(i, i) += (1. - eps) * diag(i);
        }
        double norm = nrm2(A);
        A /= norm;
//        A.print();
        Tensord ref = A * A;
    
        Tensordf B(A);
        B = B * B;
        Tensord res = B.convert();
        Tensord Delta = ref - res;
        double delta = nrm2(Delta);

        Tensorf C(A.shape_);
        cast(C, A);
        C = C * C;
        Tensord resC(C.shape_);
        cast(resC, C);
        Tensord DeltaC = ref - resC;
        double deltaC = nrm2(DeltaC);

        cout << eps << " " << delta << " " << deltaC << endl;
    }
    getchar();
}*/
