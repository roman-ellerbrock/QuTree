#include "Tensormxp.h"
#include "hostMemory.h"

using d = double;
using f = float;
template <typename T>
using host = polymorphic::hostMemory<T>;

template class mxpTensor<d, f, host>;

template <typename T, typename U>
void gemm(Tensorm<T, U>& C, const Tensorm<T, U>& A,
const Tensorm<T, U>& B) {
    T one = 1.;
    /// n^1 diagonal * diagonal
    hadamardProduct(C.diag_ ,A.diag_, B.diag_);

    /// n^2 diagonal * off & off * diagonal
    diagmm(C.off_, A.diag_, B.off_);
    mdiagm(C.off_, A.off_, B.diag_);

    /// perform n^3 in low precision
    U alpha = 1.;
    U beta = 0.;
    gemm(C.offLP_, A.offLP_ , B.offLP_, alpha, beta, blas::Op::NoTrans, blas::Op::NoTrans);
    C.addLP(); // off(hp) += off(lp)
    C.toLP();
}

template void gemm(Tensorm<d, f>& C, const Tensorm<d, f>& A, const Tensorm<d, f>& B);

template <typename T, typename U>
Tensorm<T, U> gemm(const Tensorm<T, U>& A, const Tensorm<T, U>& B) {
    Tensorm<T, U> C(A.off_.shape_);
    gemm(C, A, B);
    return C;
}

template Tensordf gemm(const Tensordf& A, const Tensordf& B);

template <typename T, typename U>
Tensorm<T, U> operator*(const Tensorm<T, U>& L, const Tensorm<T, U>& R) {
        return gemm(L, R);
}

template Tensordf operator*(const Tensordf& L, const Tensordf& R);
