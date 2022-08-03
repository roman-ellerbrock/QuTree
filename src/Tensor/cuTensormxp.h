#ifndef CUTENSORMXP_H
#define CUTENSORMXP_H
#include "Tensormxp.h"
#include "cuTensor.h"
#include "blas_queue.h"

template <typename T, typename U>
using cuTensorm = mxpTensor<T, U, polymorphic::cuMemory, blas::Queue>;

using cuTensordf = cuTensorm<double, float>;

template <typename T, typename U>
void gemm(cuTensorm<T, U>& C, const cuTensorm<T, U>& A, const cuTensorm<T, U>& B, blas::Queue& queue = qutree::queue);

template <typename T, typename U>
cuTensorm<T, U> gemm(const cuTensorm<T, U>& A, const cuTensorm<T, U>& B, blas::Queue& queue = qutree::queue);

template <typename T, typename U>
cuTensorm<T, U> operator*(const cuTensorm<T, U>& L, const cuTensorm<T, U>& R);

#endif // CUTENSORMXP_H
