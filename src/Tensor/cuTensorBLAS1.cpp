#include "Tensor/TensorBLAS1.hpp"
#include "Tensor/cuTensorBLAS1.h"
#include "Tensor/cuTensorBLAS1.cuh"
#include "blas_queue.h"

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;
using namespace polymorphic;

template<typename T>
blas::real_type<T> nrm2(const cuTensor<T>& A, size_t incr, blas::Queue& queue) {
	return abs(blas::nrm2<T>(A.shape_.totalDimension() / incr, (const T*)&(A[0]), incr, queue));
}

template f nrm2<f>(const cuTensor<f>& A, size_t inc_a, blas::Queue&);
template d nrm2<d>(const cuTensor<d>& A, size_t inc_a, blas::Queue&);
template f nrm2<cf>(const cuTensor<cf>& A, size_t inc_a, blas::Queue&);
template d nrm2<cd>(const cuTensor<cd>& A, size_t inc_a, blas::Queue&);

template void axpy<f, cuTensor<f>, blas::Queue>(const cuTensor<f>& A, cuTensor<f>& B, f alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<d, cuTensor<d>, blas::Queue>(const cuTensor<d>& A, cuTensor<d>& B, d alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cf, cuTensor<cf>, blas::Queue>(const cuTensor<cf>& A, cuTensor<cf>& B, cf alpha, size_t inc_a, size_t inc_b, blas::Queue&);
template void axpy<cd, cuTensor<cd>, blas::Queue>(const cuTensor<cd>& A, cuTensor<cd>& B, cd alpha, size_t inc_a, size_t inc_b, blas::Queue&);

/// A += B
template<typename T>
void operator+=(cuTensor<T>& A, const cuTensor<T>& B) {
	T alpha = 1.;
	axpy(B, A, alpha, 1, 1, qutree::queue);
}

template void operator+=(cuTensor<f>& A, const cuTensor<f>& B);
template void operator+=(cuTensor<d>& A, const cuTensor<d>& B);
template void operator+=(cuTensor<cf>& A, const cuTensor<cf>& B);
template void operator+=(cuTensor<cd>& A, const cuTensor<cd>& B);

template<typename T>
void operator-=(cuTensor<T>& A, const cuTensor<T>& B) {
	T alpha = -1.;
	axpy<T>(B, A, alpha, 1, 1, qutree::queue);
}

template void operator-=(cuTensor<f>& A, const cuTensor<f>& B);
template void operator-=(cuTensor<d>& A, const cuTensor<d>& B);
template void operator-=(cuTensor<cf>& A, const cuTensor<cf>& B);
template void operator-=(cuTensor<cd>& A, const cuTensor<cd>& B);

template<typename T>
double residual(cuTensor<T> A, const cuTensor<T>& B, blas::Queue& queue) {
	A -= B;
	return abs(nrm2<T>(A, 1, queue));
}

template double residual<f>(cuTensorf A, const cuTensorf& B, blas::Queue& queue);
template double residual<d>(cuTensord A, const cuTensord& B, blas::Queue& queue);
template double residual<cf>(cuTensorcf A, const cuTensorcf& B, blas::Queue& queue);
template double residual<cd>(cuTensorcd A, const cuTensorcd& B, blas::Queue& queue);

/// cast vector from double to float
template <>
void cast<float, double, cuMemory>(cuTensorf& L, const cuTensord& R) {
    size_t nL = L.shape_.totalDimension();
    size_t nR = R.shape_.totalDimension();
//    errorif(nL != nR);
    cudaCast(L.data(), R.data(), nL);
}

template <>
void cast<double, float, cuMemory>(cuTensord& L, const cuTensorf& R) {
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

template<typename T>
void gemm(cuTensor<T>& c, const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b, blas::Queue& queue) {

	/// m: Number of rows of the matrix C and \(op(A)\). m >= 0.
	size_t m = nrows(a.shape_, op_a);
	assert(nrows(c.shape_) == m);

	/// n: Number of cols of the matrix C and \(op(B)\). n >= 0.
	size_t n = ncols(b.shape_, op_b);
	assert(ncols(c.shape_) == n);

	size_t k = ncols(a.shape_, op_a);
	assert(nrows(b.shape_, op_b) == k);

	size_t lda = (op_a == blas::Op::NoTrans) ? m : k;
	size_t ldb = (op_b == blas::Op::NoTrans) ? k : n;
	size_t ldc = m;

	blas::gemm(blas::Layout::ColMajor, op_a, op_b, m, n, k,
		alpha, a.data(), lda, b.data(), ldb, beta, c.data(), ldc,
		queue);
}

template void gemm(cuTensor<f>& c, const cuTensor<f>& a, const cuTensor<f>& b,
	f alpha, f beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template void gemm(cuTensor<d>& c, const cuTensor<d>& a, const cuTensor<d>& b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template void gemm(cuTensor<cf>& c, const cuTensor<cf>& a, const cuTensor<cf>& b,
	cf alpha, cf beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template void gemm(cuTensor<cd>& c, const cuTensor<cd>& a, const cuTensor<cd>& b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);

template<typename T>
cuTensor<T> gemm(const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b, blas::Queue& queue) {
    TensorShape shape = a.shape_;
    shape.back() = b.shape_.back();
    cuTensor<T> c(shape);
    gemm(c, a, b, alpha, beta, op_a, op_b, queue);
    return c;
}

template cuTensor<f> gemm(const cuTensor<f>& a, const cuTensor<f>& b,
	f alpha, f beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template cuTensor<d> gemm(const cuTensor<d>& a, const cuTensor<d>& b,
	d alpha, d beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template cuTensor<cf> gemm(const cuTensor<cf>& a, const cuTensor<cf>& b,
	cf alpha, cf beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);
template cuTensor<cd> gemm(const cuTensor<cd>& a, const cuTensor<cd>& b,
	cd alpha, cd beta, blas::Op op_a, blas::Op op_b, blas::Queue& queue);


template<typename T>
cuTensor<T> operator*(const cuTensor<T>& a, const cuTensor<T>& b) {
    T alpha = 1.;
    T beta = 0.;
    blas::Op no = blas::Op::NoTrans;
    return gemm(a, b, alpha, beta, no, no, qutree::queue);
}

template cuTensor<f> operator*(const cuTensor<f>& a, const cuTensor<f>& b);
template cuTensor<d> operator*(const cuTensor<d>& a, const cuTensor<d>& b);
template cuTensor<cf> operator*(const cuTensor<cf>& a, const cuTensor<cf>& b);
template cuTensor<cd> operator*(const cuTensor<cd>& a, const cuTensor<cd>& b);

template <typename T>
void diagmPlusmdiag(cuTensor<T>& C, const cuTensor<T>& A, const cuTensor<T>& diag, T factor) {
	size_t nrow = diag.shape_.front();
	cudaDiagmPlusmdiag(C.data(), A.data(), diag.data(), factor, nrow, nrow);
}

template void diagmPlusmdiag(cuTensord& C, const cuTensord& A, const cuTensord& diag, d factor);

