#ifndef TENSORBLAS1_HPP
#define TENSORBLAS1_HPP
#include "TensorBLAS1.h"
#include "Util/blaspp_device_extension.h"

/// = ||A||_2
template<typename T>
blas::real_type<T> nrm2(const Tensor<T>& A, size_t incr) {
	return abs(blas::nrm2<T>(A.shape_.totalDimension() / incr, (const T*)&(A[0]), incr));
}

/// b -> alpha * a + b
template<typename T, class Tensor, class ...Queue>
void axpy(const Tensor& A, Tensor& B, T alpha, size_t inc_a, size_t inc_b, Queue& ...queue) {
	size_t n = A.shape_.totalDimension() / inc_a;
	blas::axpy<T>(n, alpha, A.data(), inc_a, B.data(), inc_b, queue...);
}

template<typename T, template <typename> class Dev, class ...Queue>
void diagonal(Tensor<T, Dev>& diag, const Tensor<T, Dev>& A, Queue& ...queue) {
	size_t nrow = nrows(A.shape_);
	size_t ncol = ncols(A.shape_);
	size_t n = (nrow < ncol) ? nrow : ncol;

	blas::copy(n, A.data(), n + 1, diag.data(), 1, queue...);
}

template<typename T, template <typename> class Dev, class ...Queue>
Tensor<T, Dev> diagonal(const Tensor<T, Dev>& A, Queue& ...queue) {
	size_t nrow = nrows(A.shape_);
	size_t ncol = ncols(A.shape_);
	size_t n = (nrow < ncol) ? nrow : ncol;

	Tensor<T, Dev> diag({n});
	diagonal(diag, A, queue...);
	return diag;
}

template<typename T, template <typename> class Dev, class ...Queue>
void addDiagonal(Tensor<T, Dev>& B, const Tensor<T, Dev>& diag, T alpha, Queue& ...queue) {
	size_t incb = B.shape_.lastBefore() + 1;
	size_t n = diag.shape_.totalDimension();
	blas::axpy(n, alpha, diag.data(), 1, B.data(), incb, queue...);
}

template <typename T, template <typename> class Dev, class ... Queue>
void offDiagonal(Tensor<T, Dev>& off, const Tensor<T, Dev>& full, Queue& ...queue) {
    off = full;
	Tensor<T, Dev> diag = diagonal(full, queue...);
	T alpha = -1.;
	addDiagonal(off, diag, alpha, queue...);
}


#endif // TENSORBLAS1_HPP