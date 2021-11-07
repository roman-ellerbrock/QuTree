//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Core/TensorBLAS1.h"

typedef complex<double> cd;
typedef double d;

/// = ||A||_2
template<typename T>
T nrm2(const Tensor<T>& A, size_t incr) {
	return abs(blas::nrm2(A.shape_.totalDimension() / incr, &(A[0]), incr));
}

template cd nrm2(const Tensor<cd>&, size_t);
template d nrm2(const Tensor<d>&, size_t);


/// b -> alpha * a + b
template<typename T>
void axpy(const Tensor<T>& A, Tensor<T>& B, T alpha, size_t inc_a, size_t inc_b) {
	size_t n = A.shape_.totalDimension() / inc_a;
	blas::axpy(n, alpha, &(A[0]), inc_a, &(B[0]), inc_b);
}

template void axpy(const Tensor<cd>&, Tensor<cd>&, cd, size_t, size_t);
template void axpy(const Tensor<d>&, Tensor<d>&, d, size_t, size_t);


/// A += B
template<typename T>
void operator+=(Tensor<T>& A, const Tensor<T>& B) {
	T alpha = 1.;
	axpy(B, A, alpha);
}

template void operator+=(Tensor<cd>& A, const Tensor<cd>& B);
template void operator+=(Tensor<d>& A, const Tensor<d>& B);


/// A -= B
template<typename T>
void operator-=(Tensor<T>& A, const Tensor<T>& B) {
	T alpha = -1.;
	axpy(B, A, alpha);
}

template void operator-=(Tensor<cd>& A, const Tensor<cd>& B);
template void operator-=(Tensor<d>& A, const Tensor<d>& B);


/// || A - B ||_2
template<typename T>
T residual(Tensor<T>& A, const Tensor<T>& B) {
	A -= B;
	return nrm2(A);
}

template d residual(Tensor<d>& A, const Tensor<d>& B);
template cd residual(Tensor<cd>& A, const Tensor<cd>& B);


/// = A + B
template<typename T>
Tensor<T> operator+(Tensor<T>& A, const Tensor<T>& B) {
	A += B;
	return A;
}

template Tensor<d> operator+(Tensor<d>& A, const Tensor<d>& B);
template Tensor<cd> operator+(Tensor<cd>& A, const Tensor<cd>& B);

/// = A + B
template<typename T>
Tensor<T> operator-(Tensor<T>& A, const Tensor<T>& B) {
	A -= B;
	return A;
}

template Tensor<d> operator-(Tensor<d>& A, const Tensor<d>& B);
template Tensor<cd> operator-(Tensor<cd>& A, const Tensor<cd>& B);





