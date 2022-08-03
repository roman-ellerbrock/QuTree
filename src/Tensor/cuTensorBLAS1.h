#ifndef CUTENSORBLAS1_H
#define CUTENSORBLAS1_H
#include "Tensor/TensorBLAS1.h"
#include "Tensor/cuTensor.h"
#include "blas_queue.h"

/**
 * \brief Calculate 2-norm of a Tensor
 * @param incr increment in loop
 * @return returns the 2-norm of a Tensor, i.e. ||A||_2
 */
template<typename T>
blas::real_type<T> nrm2(const cuTensor<T>& A, size_t incr = 1, blas::Queue& queue = qutree::queue);

template<typename T>
double residual(cuTensor<T> A, const cuTensor<T>& B, blas::Queue& queue = qutree::queue);

template<typename T>
void operator+=(cuTensor<T>& A, const cuTensor<T>& B);

template<typename T>
void operator-=(cuTensor<T>& A, const cuTensor<T>& B);

template<typename T>
void gemm(cuTensor<T>& c, const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b,
	blas::Queue& queue = qutree::queue);

template<typename T>
cuTensor<T> gemm(const cuTensor<T>& a, const cuTensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b,
	blas::Queue& queue = qutree::queue);

template <typename T>
cuTensor<T> operator*(const cuTensor<T>& L, const cuTensor<T>& R);

template <typename T>
void diagmPlusmdiag(cuTensor<T>& C, const cuTensor<T>& A, const cuTensor<T>& diag, T factor = (T) 1.);

#endif // CUTENSORBLAS1_H