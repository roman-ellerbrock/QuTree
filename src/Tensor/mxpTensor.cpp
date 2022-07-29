#include "mxpTensor.h"
#include "hostMemory.h"

using d = double;
using f = float;
template <typename T>
using host = polymorphic::hostMemory<T>;

template class mxpTensor<d, f, host>;

template <typename T, typename U, template <typename> class Dev>
void gemm(mxpTensor<T, U, Dev>& C, const mxpTensor<T, U, Dev>& A,
const mxpTensor<T, U, Dev>& B) {
    /// n^1 diagonal * diagonal
    hadamardProduct(C.diag_ ,A.diag_, B.diag_);

    /// n^2 diagonal * off & off * diagonal
    diagmm(C.off_, A.diag_, B.off_);
    mdiagm(C.off_, A.off_, B.diag_);

    /// perform n^3 in low precision
    gemm(C.offLP_, A.offLP_ , B.offLP_);
    C.addLP(); // off(hp) += off(lp)
    C.toLP();
}

template void gemm(mxpTensor<d, f, host>& C, const mxpTensor<d, f, host>& A, const mxpTensor<d, f, host>& B);

template <typename T, typename U, template <typename> class Dev>
mxpTensor<T, U, Dev> gemm(const mxpTensor<T, U, Dev>& A, const mxpTensor<T, U, Dev>& B) {
    mxpTensor<T, U, Dev> C(A.off_.shape_);
    gemm(C, A, B);
    return C;
}

template mxpTensor<d, f, host> gemm(const mxpTensor<d, f, host>& A, const mxpTensor<d, f, host>& B);
