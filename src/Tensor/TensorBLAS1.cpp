//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Tensor/TensorBLAS1.h"
#include "Tensor/Tensor.hpp"
#include "stdafx.h"
#include <math.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

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
	blas::axpy(n, alpha, A.coeffs_, inc_a, B.coeffs_, inc_b);
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
double residual(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	return abs(nrm2(A));
}

template d residual(Tensor<d> A, const Tensor<d>& B);
template d residual(Tensor<cd> A, const Tensor<cd>& B);

/// = A + B
template<typename T>
Tensor<T> operator+(Tensor<T> A, const Tensor<T>& B) {
	A += B;
	return A;
}

template Tensor<d> operator+(Tensor<d> A, const Tensor<d>& B);
template Tensor<cd> operator+(Tensor<cd> A, const Tensor<cd>& B);

/// = A + B
template<typename T>
Tensor<T> operator-(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	return A;
}

template Tensor<d> operator-(Tensor<d> A, const Tensor<d>& B);
template Tensor<cd> operator-(Tensor<cd> A, const Tensor<cd>& B);

/// A *= alpha
template<typename T, typename U>
Tensor<T>& operator*=(Tensor<T>& A, const U alpha) {
	size_t n = A.shape_.totalDimension();
	blas::scal(n, (T) alpha, A.coeffs_, 1);
	return A;
}

template Tensor<d>& operator*=(Tensor<d>& A, const d alpha);
template Tensor<cd>& operator*=(Tensor<cd>& A, const cd alpha);
template Tensor<cd>& operator*=(Tensor<cd>& A, const d alpha);

///  = alpha * A
template<typename T, typename U>
Tensor<T> operator*(const U alpha, Tensor<T> A) {
	return A *= alpha;
}

template Tensor<d> operator*(const d alpha, Tensor<d> A);
template Tensor<cd> operator*(const cd alpha, Tensor<cd> A);
template Tensor<cd> operator*(const d alpha, Tensor<cd> A);

///  = A * alpha
template<typename T, typename U>
Tensor<T> operator*(Tensor<T> A, const U alpha) {
	return A *= alpha;
}

template Tensor<d> operator*(Tensor<d> A, const d alpha);
template Tensor<cd> operator*(Tensor<cd> A, const cd alpha);
template Tensor<cd> operator*(Tensor<cd> A, const d alpha);

///  A /= alpha
template<typename T, typename U>
Tensor<T>& operator/=(Tensor<T>& A, const U alpha) {
	return A *= (1. / alpha);
}

template Tensor<d>& operator/=(Tensor<d>& A, const d alpha);
template Tensor<cd>& operator/=(Tensor<cd>& A, const cd alpha);
template Tensor<cd>& operator/=(Tensor<cd>& A, const d alpha);

///  = A / alpha
template<typename T, typename U>
Tensor<T> operator/(Tensor<T> A, const U alpha) {
	return A /= alpha;
}

template Tensor<d> operator/(Tensor<d> A, const d alpha);
template Tensor<cd> operator/(Tensor<cd> A, const cd alpha);
template Tensor<cd> operator/(Tensor<cd> A, const d alpha);

template<typename T>
Tensor<T> productElementwise(const Tensor<T>& A, const Tensor<T>& B) {
	assert(A.shape_.totalDimension() == B.shape_.totalDimension());
	Tensor<T> C(A.shape_);
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		C(i) = A(i) * B(i);
	}
	return C;
}

template Tensor<cd> productElementwise(const Tensor<cd>& A, const Tensor<cd>& B);
template Tensor<d> productElementwise(const Tensor<d>& A, const Tensor<d>& B);

template<typename T>
Tensor<T> conj(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = conj(A(i));
	}
	return A;
}

template Tensor<cd> conj(Tensor<cd> A);
template Tensor<d> conj(Tensor<d> A);

template<typename T>
Tensor<T> diagonal(const Tensor<T>& A) {
	size_t nrow = nrows(A.shape_);
	size_t ncol = ncols(A.shape_);
	size_t n = (nrow < ncol) ? nrow : ncol;
	Tensor<T> diag({n});
	for (size_t i = 0; i < n; ++i) {
		diag(i) = A(i, i);
	}
	return diag;
}

template Tensor<cd> diagonal(const Tensor<cd>& A);
template Tensor<d> diagonal(const Tensor<d>& A);

template<typename T>
T trace(const Tensor<T>& A) {
	Tensor<T> diag = diagonal(A);
	T tr = 0.;
	for (size_t i = 0; i < diag.shape_.totalDimension(); ++i) {
		tr += diag(i);
	}
	return tr;
}

template cd trace(const Tensor<cd>& A);
template d trace(const Tensor<d>& A);

/// dest[dim1, dim2] = beta * dest[dim1,dim2] + A[dim2, dim1]
template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta) {
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
			dest[j + dim2 * i] = beta * dest[j + dim2 * i] + src[i + dim1 * j];
		}
	}
}

template void transpose(d *dest, const d *src, size_t dim1, size_t dim2, d beta);
template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2, cd beta);

/// A[dim1, dim2] --> A[dim2, dim1]
template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2) {
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
			dest[j + dim2 * i] = src[i + dim1 * j];
		}
	}
}

template void transpose(d *dest, const d *src, size_t dim1, size_t dim2);
template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2);


/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
void transpose(Tensor<T>& dest, const Tensor<T>& src) {
	size_t dim1 = src.shape_.lastBefore();
	size_t dim2 = src.shape_.lastDimension();
	transpose(dest.coeffs_, src.coeffs_, dim1, dim2);
	dest.shape_ = transposeToFront(src.shape_, src.shape_.lastIdx());
}

template void transpose(Tensor<d>& dest, const Tensor<d>& src);
template void transpose(Tensor<cd>& dest, const Tensor<cd>& src);


/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
Tensor<T> transpose(const Tensor<T>& src) {
	Tensor<T> dest(src.shape_);
	transpose(dest, src);
	return dest;
}

template Tensor<d> transpose(const Tensor<d>& src);
template Tensor<cd> transpose(const Tensor<cd>& src);

/// \brief perform matrix adjoint A[bef, last] --> A*[last, bef]
template<typename T>
void adjoint(Tensor<T>& dest, const Tensor<T>& src) {
	transpose(dest, src);
	dest = conj(dest);
}

template void adjoint(Tensor<d>& dest, const Tensor<d>& src);
template void adjoint(Tensor<cd>& dest, const Tensor<cd>& src);

/// \brief perform matrix adjoint A[bef, last] --> A*[last, bef]
template<typename T>
Tensor<T> adjoint(const Tensor<T>& src) {
	Tensor<T> dest(src.shape_);
	adjoint(dest, src);
	return dest;
}

template Tensor<d> adjoint(const Tensor<d>& src);
template Tensor<cd> adjoint(const Tensor<cd>& src);

/// dest[b, a, c] = src[a, b, c]
template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C) {
	for (size_t c = 0; c < C; ++c) {
		transpose(&dest[c * A * B], &src[c * A * B], A, B);
	}
}

template void transposeAB(cd *dest, const cd *src, size_t A, size_t B, size_t C);
template void transposeAB(d *dest, const d *src, size_t A, size_t B, size_t C);

/// dest[a, c, b] = src[a, b, c]
template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C) {
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			memcpy(&dest[c * A + b * A * C], &src[b * A + c * A * B], A * sizeof(T));
		}
	}
}

template void transposeBC(d *dest, const d *src, size_t A, size_t B, size_t C);
template void transposeBC(cd *dest, const cd *src, size_t A, size_t B, size_t C);

/// dest[...,...,k] = src[...,k,...] if back == false
template<typename T>
void transpose(Tensor<T>& dest, const Tensor<T>& src, size_t k, bool back) {
	assert(dest.shape_.totalDimension() == src.shape_.totalDimension());
	size_t a = src.shape_.before(k);
	size_t b = src.shape_[k];
	size_t c = src.shape_.after(k);
	if (back) {
		b = src.shape_.lastBefore() / src.shape_.before(k);
		c = src.shape_.lastDimension();
	}
	transposeBC(dest.coeffs_, src.coeffs_, a, b, c);

	/// transpose shape
	dest.shape_ = transpose(src.shape_, k, back);
}

template void transpose(Tensor<cd>& dest, const Tensor<cd>& src, size_t k, bool back);
template void transpose(Tensor<d>& dest, const Tensor<d>& src, size_t k, bool back);

template<typename T>
Tensor<T> transpose(const Tensor<T>& src, size_t k, bool back) {
	TensorShape shapeT = transpose(src.shape_, k, back);
	Tensor<T> dest(shapeT);
	transpose(dest, src, k, back);
	return dest;
}

template Tensor<cd> transpose(const Tensor<cd>& src, size_t k, bool back);
template Tensor<d> transpose(const Tensor<d>& src, size_t k, bool back);

template<typename T>
double isCloseToIdentity(const Tensor<T>& A) {
	double eps = 0.;
	for(size_t j = 0; j < A.shape_.lastDimension(); ++j) {
		for (size_t i = 0; i < A.shape_.lastBefore(); ++i) {
			if (i == j) {
				eps += pow(abs(A(i, j) - 1.), 2);
			} else {
				eps += pow(abs(A(i, j)), 2);
			}
		}
	}
	return sqrt(eps);
}

template double isCloseToIdentity(const Tensor<d>& A);
template double isCloseToIdentity(const Tensor<cd>& A);

size_t nrows(const TensorShape& shape, blas::Op op) {
	size_t n = shape.lastBefore();
	if (op != blas::Op::NoTrans) { n = shape.lastDimension(); }
	return n;
}

size_t ncols(const TensorShape& shape, blas::Op op) {
	size_t n = shape.lastDimension();
	if (op != blas::Op::NoTrans) { n = shape.lastBefore(); }
	return n;
}

template<typename T, typename U>
void vectorTensor(Tensor<T>& B, const Tensor<U>& vec, size_t k) {
	for (size_t aft = 0; aft < B.shape_.after(k); ++aft) {
		for (size_t act = 0; act < B.shape_[k]; ++act) {
			for (size_t bef = 0; bef < B.shape_.before(k); ++bef) {
				B(bef, act, aft, k) *= vec(act);
			}
		}
	}
}

template void vectorTensor(Tensor<cd>& B, const Tensor<cd>& a, size_t k);
template void vectorTensor(Tensor<cd>& B, const Tensor<d>& a, size_t k);
template void vectorTensor(Tensor<d>& B, const Tensor<d>& a, size_t k);

