#include "cuTensormxp.h"
#include "cuTensorBLAS1.h"

using d = double;
using f = float;

template <typename T>
using cuda = polymorphic::cuMemory<T>;

template class mxpTensor<d, f, cuda, blas::Queue>;

template <typename T, typename U>
void gemm(cuTensorm<T, U>& C, const cuTensorm<T, U>& A, const cuTensorm<T, U>& B, blas::Queue& queue) {
    T one = 1.;
    /// n^1 diagonal * diagonal
    hadamardProduct(C.diag_ ,A.diag_, B.diag_);

    /// n^2 diagonal * off & off * diagonal
    diagmm(C.off_, A.diag_, B.off_);
    mdiagm(C.off_, A.off_, B.diag_);
//    diagmPlusmdiag(C.off_, A.off_, B.diag_);

    /// perform n^3 in low precision
    U alpha = 1.;
    U beta = 0.;
    gemm(C.offLP_, A.offLP_ , B.offLP_, alpha, beta, blas::Op::NoTrans, blas::Op::NoTrans, queue);
    C.addLP(); // off(hp) += off(lp)
    C.toLP();
}

template void gemm(cuTensordf& C, const cuTensordf& A, const cuTensordf& B, blas::Queue& queue);

template <typename T, typename U>
cuTensorm<T, U> gemm(const cuTensorm<T, U>& A, const cuTensorm<T, U>& B, blas::Queue& queue) {
    cuTensorm<T, U> C(A.off_.shape_);
    gemm(C, A, B, queue);
    return C;
}

template cuTensordf gemm(const cuTensordf& A, const cuTensordf& B, blas::Queue& queue);

template <typename T, typename U>
cuTensorm<T, U> operator*(const cuTensorm<T, U>& L, const cuTensorm<T, U>& R) {
        return gemm(L, R, qutree::queue);
}

template cuTensordf operator*(const cuTensordf& L, const cuTensordf& R);
