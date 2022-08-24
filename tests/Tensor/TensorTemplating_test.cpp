#include <gtest/gtest.h>
#include "Tensor/Tensor.h"
#include "Tensor/cuTensorBLAS1.h"
#include "Tensor/Tensormxp.h"
#include "Tensor/cuTensormxp.h"
#include <chrono>

template <class Tensor>
void purification(Tensor& D, size_t nIter) {

    Tensor D2(D.shape_);
    for (size_t i = 0; i < nIter; ++i) {
        D2 = D * D;
        std::swap(D2, D);
    }
}

Tensord initialDensity(size_t n) {
    Tensord D = deltad({n, n});
    D += 1e-3 * randomd({n, n});
    return D;
}

template <class Tensor>
void benchmark(function<void(Tensor&, size_t)>& f, Tensor& A, size_t niter, string name) {
	using namespace chrono;
    time_point<steady_clock> start = steady_clock::now();
    f(A, niter);
    time_point<steady_clock> end = steady_clock::now();
    double ms = duration_cast<nanoseconds>(end - start).count() / ((double) niter * 1000000);
    cout << "function: " << name << " - time: " << ms << " ms\n";
}

TEST(TensorIntegration, purification) {
    constexpr int n = 5000;
    size_t niter = 3;

    /// Tensor<double, host>
    Tensord hD = initialDensity(n);
    function<void(Tensord&, size_t)> purificationd = purification<Tensord>;
    benchmark(purificationd, hD, niter, "Tensord");

    /// A = A_diag + A_off
    /// A^2 = A_diag^2 + A_diag * A_off + A_off * A_diag + A_off^2 [N^3]
    Tensordf mxpD(hD);
    function<void(Tensordf&, size_t)> purificationdf = purification<Tensordf>;
    benchmark(purificationdf, mxpD, niter, "Tensordf");

    cuTensord D = transferToGPUd(hD);
    function<void(cuTensord&, size_t)> cupurificationd = purification<cuTensord>;
    benchmark(cupurificationd, D, niter, "cuTensord");

    cuTensordf cuDdf(D, qutree::queue);
    function<void(cuTensordf&, size_t)> cupurificationdf = purification<cuTensordf>;
    benchmark(cupurificationdf, cuDdf, niter, "cuTensordf");
    getchar();
}
