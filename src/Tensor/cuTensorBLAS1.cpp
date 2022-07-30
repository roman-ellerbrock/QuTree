#include "Tensor/TensorBLAS1.hpp"
#include "Tensor/cuTensorBLAS1.h"
#include "Tensor/cuTensorBLAS1.cuh"
#include "blas_queue.h"

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;
using namespace polymorphic;
//template double nrm2<d, cuTensor, blas::Queue>(const cuTensor<d>& A, size_t inc_a, blas::Queue&);

template void axpy<f, cuTensor<f>, blas::Queue>(const cuTensor<f>& A, cuTensor<f>& B, f alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<d, cuTensor<d>, blas::Queue>(const cuTensor<d>& A, cuTensor<d>& B, d alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cf, cuTensor<cf>, blas::Queue>(const cuTensor<cf>& A, cuTensor<cf>& B, cf alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cd, cuTensor<cd>, blas::Queue>(const cuTensor<cd>& A, cuTensor<cd>& B, cd alpha, size_t inc_a, size_t inc_b, blas::Queue&);

/// A += B
/*template<typename T>
void operator+=(cuTensor<T>& A, const cuTensor<T>& B) {
	T alpha = 1.;
	axpy(B, A, alpha, 1, 1, qutree::queue);
}

template void operator+=(cuTensor<f>& A, const cuTensor<f>& B);
template void operator+=(cuTensor<d>& A, const cuTensor<d>& B);
template void operator+=(cuTensor<cf>& A, const cuTensor<cf>& B);
template void operator+=(cuTensor<cd>& A, const cuTensor<cd>& B);
*/
/// cast vector from double to float
template <>
void cast<float, double, cuMemory>(cuTensorf& L, const cuTensord& R) {
    size_t nL = L.shape_.totalDimension();
    size_t nR = R.shape_.totalDimension();
//    errorif(nL != nR);
    cudaCast(L.data(), R.data(), nL);
}

template cuTensorf diagonal(const cuTensorf& A, blas::Queue& queue);
template cuTensord diagonal(const cuTensord& A, blas::Queue& queue);
template cuTensorcf diagonal(const cuTensorcf& A, blas::Queue& queue);
template cuTensorcd diagonal(const cuTensorcd& A, blas::Queue& queue);

template void addDiagonal(cuTensorf& B, const cuTensorf& diag, f alpha, blas::Queue& queue);
template void addDiagonal(cuTensord& B, const cuTensord& diag, d alpha, blas::Queue& queue);
template void addDiagonal(cuTensorcf& B, const cuTensorcf& diag, cf alpha, blas::Queue& queue);
template void addDiagonal(cuTensorcd& B, const cuTensorcd& diag, cd alpha, blas::Queue& queue);

template void offDiagonal(cuTensor<f>& off, const cuTensor<f>& full, blas::Queue& queue);
template void offDiagonal(cuTensor<d>& off, const cuTensor<d>& full, blas::Queue& queue);
template void offDiagonal(cuTensor<cf>& off, const cuTensor<cf>& full, blas::Queue& queue);
template void offDiagonal(cuTensor<cd>& off, const cuTensor<cd>& full, blas::Queue& queue);

template <>
void hadamardProduct<d, cuMemory>(cuTensord& C, const cuTensord& A, const cuTensord& B) {
    size_t nC = C.shape_.totalDimension();
    size_t nA = A.shape_.totalDimension();
    size_t nB = B.shape_.totalDimension();
//    errorif(nC != nA);
//    errorif(nC != nB);
    cudaHadamardProduct(C.data(), A.data(), B.data(), nC);
}

template <>
void diagmm<d, cuMemory>(cuTensord& C, const cuTensord& dA, const cuTensord& B, d factor) {
    size_t nrow = B.shape_.lastBefore();
    size_t ncol = B.shape_.lastDimension();
    /// add errors
    cudaDiagMatrixMatrix(C.data(), dA.data(), B.data(), factor, nrow, ncol);
}

template <>
void mdiagm<d, cuMemory>(cuTensord& C, const cuTensord& A, const cuTensord& dB, d factor) {
    size_t nrow = A.shape_.lastBefore();
    size_t ncol = A.shape_.lastDimension();
    /// add errors
    cudaMatrixDiagMatrix(C.data(), A.data(), dB.data(), factor, nrow, ncol);
}


