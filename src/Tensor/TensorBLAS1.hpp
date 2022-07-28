#ifndef TENSORBLAS1_HPP
#define TENSORBLAS1_HPP
#include "TensorBLAS1.h"
#include "Util/blaspp_device_extension.h"

/// = ||A||_2
template<typename T, template <typename> class Tensor, class ...Queue>
double nrm2(const Tensor<T>& A, size_t incr, Queue& ... queue) {
	return abs(blas::nrm2<T>(A.shape_.totalDimension() / incr, (const T*)&(A[0]), incr, queue...));
}

/// b -> alpha * a + b
template<typename T, class Tensor, class ...Queue>
void axpy(const Tensor& A, Tensor& B, T alpha, size_t inc_a, size_t inc_b, Queue& ...queue) {
	size_t n = A.shape_.totalDimension() / inc_a;
	blas::axpy<T>(n, alpha, A.data(), inc_a, B.data(), inc_b, queue...);
}


#endif // TENSORBLAS1_HPP