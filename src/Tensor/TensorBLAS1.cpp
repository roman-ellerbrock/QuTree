//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Tensor/TensorBLAS1.h"
#include "Tensor/Tensor.hpp"
#include "stdafx.h"
#include <math.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

/// = ||A||_2
template<class Tensor, class ...Queue>
double nrm2(const Tensor& A, size_t incr, Queue& ... queue) {
	return abs(blas::nrm2(A.shape_.totalDimension() / incr, &(A[0]), incr, queue...));
}

template double nrm2<Tensorf>(const Tensorf& A, size_t incr);
template double nrm2<Tensord>(const Tensord& A, size_t incr);
template double nrm2<Tensorcf>(const Tensorcf& A, size_t incr);
template double nrm2<Tensorcd>(const Tensorcd& A, size_t incr);

/// b -> alpha * a + b
template<typename T, class Tensor, class ...Queue>
void axpy(const Tensor& A, Tensor& B, T alpha, size_t inc_a, size_t inc_b, Queue& ...queue) {
	size_t n = A.shape_.totalDimension() / inc_a;
	blas::axpy(n, alpha, A.data(), inc_a, B.data(), inc_b, queue...);
}

template void axpy<f, Tensor<f>>(const Tensor<f>& A, Tensor<f>& B, f alpha, size_t inc_a, size_t inc_b);
template void axpy<d, Tensor<d>>(const Tensor<d>& A, Tensor<d>& B, d alpha, size_t inc_a, size_t inc_b);
template void axpy<cf, Tensor<cf>>(const Tensor<cf>& A, Tensor<cf>& B, cf alpha, size_t inc_a, size_t inc_b);
template void axpy<cd, Tensor<cd>>(const Tensor<cd>& A, Tensor<cd>& B, cd alpha, size_t inc_a, size_t inc_b);

/// A += B
template<typename T>
void operator+=(Tensor<T>& A, const Tensor<T>& B) {
	T alpha = 1.;
	axpy(B, A, alpha);
}

template void operator+=(Tensor<f>& A, const Tensor<f>& B);
template void operator+=(Tensor<d>& A, const Tensor<d>& B);
template void operator+=(Tensor<cf>& A, const Tensor<cf>& B);
template void operator+=(Tensor<cd>& A, const Tensor<cd>& B);

/// A -= B
template<typename T>
void operator-=(Tensor<T>& A, const Tensor<T>& B) {
	T alpha = -1.;
	axpy(B, A, alpha);
}

template void operator-=(Tensor<f>& A, const Tensor<f>& B);
template void operator-=(Tensor<d>& A, const Tensor<d>& B);
template void operator-=(Tensor<cf>& A, const Tensor<cf>& B);
template void operator-=(Tensor<cd>& A, const Tensor<cd>& B);

/// || A - B ||_2
template<class Tensor, class ...Queue>
double residual(Tensor A, const Tensor& B, Queue& ...queue) {
	A -= B;
	return abs(nrm2<Tensor, Queue...>(A, queue...));
}

template double residual<Tensorf>(Tensorf A, const Tensorf& B);
template double residual<Tensord>(Tensord A, const Tensord& B);
template double residual<Tensorcf>(Tensorcf A, const Tensorcf& B);
template double residual<Tensorcd>(Tensorcd A, const Tensorcd& B);

/// = A + B
template<typename T>
Tensor<T> operator+(Tensor<T> A, const Tensor<T>& B) {
	A += B;
	return A;
}

template Tensor<f> operator+(Tensor<f> A, const Tensor<f>& B);
template Tensor<d> operator+(Tensor<d> A, const Tensor<d>& B);
template Tensor<cf> operator+(Tensor<cf> A, const Tensor<cf>& B);
template Tensor<cd> operator+(Tensor<cd> A, const Tensor<cd>& B);

/// = A + B
template<typename T>
Tensor<T> operator-(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	return A;
}

template Tensor<f> operator-(Tensor<f> A, const Tensor<f>& B);
template Tensor<d> operator-(Tensor<d> A, const Tensor<d>& B);
template Tensor<cf> operator-(Tensor<cf> A, const Tensor<cf>& B);
template Tensor<cd> operator-(Tensor<cd> A, const Tensor<cd>& B);

/// A *= alpha
template<typename T, typename U>
Tensor<T>& operator*=(Tensor<T>& A, const U alpha) {
	size_t n = A.shape_.totalDimension();
	blas::scal(n, (T) alpha, A.data(), 1);
	return A;
}

template Tensor<f>& operator*=(Tensor<f>& A, const f alpha);
template Tensor<d>& operator*=(Tensor<d>& A, const d alpha);
template Tensor<cf>& operator*=(Tensor<cf>& A, const cf alpha);
template Tensor<cf>& operator*=(Tensor<cf>& A, const f alpha);
template Tensor<cd>& operator*=(Tensor<cd>& A, const cd alpha);
template Tensor<cd>& operator*=(Tensor<cd>& A, const d alpha);

///  = alpha * A
template<typename T, typename U>
Tensor<T> operator*(const U alpha, Tensor<T> A) {
	return A *= alpha;
}

template Tensor<f> operator*(const f alpha, Tensor<f> A);
template Tensor<d> operator*(const d alpha, Tensor<d> A);
template Tensor<cf> operator*(const cf alpha, Tensor<cf> A);
template Tensor<cd> operator*(const cd alpha, Tensor<cd> A);
template Tensor<cf> operator*(const f alpha, Tensor<cf> A);
template Tensor<cd> operator*(const d alpha, Tensor<cd> A);


///  = A * alpha
template<typename T, typename U>
Tensor<T> operator*(Tensor<T> A, const U alpha) {
	return A *= alpha;
}

template Tensor<f> operator*(Tensor<f> A, const f alpha);
template Tensor<d> operator*(Tensor<d> A, const d alpha);
template Tensor<cf> operator*(Tensor<cf> A, const cf alpha);
template Tensor<cf> operator*(Tensor<cf> A, const f alpha);
template Tensor<cd> operator*(Tensor<cd> A, const cd alpha);
template Tensor<cd> operator*(Tensor<cd> A, const d alpha);

///  A /= alpha
template<typename T, typename U>
Tensor<T>& operator/=(Tensor<T>& A, const U alpha) {
	return A *= (T) ((U)1. / alpha);
}

template Tensor<f>& operator/=(Tensor<f>& A, const f alpha);
template Tensor<d>& operator/=(Tensor<d>& A, const d alpha);
template Tensor<cf>& operator/=(Tensor<cf>& A, const f alpha);
template Tensor<cd>& operator/=(Tensor<cd>& A, const d alpha);
template Tensor<cf>& operator/=(Tensor<cf>& A, const cf alpha);
template Tensor<cd>& operator/=(Tensor<cd>& A, const cd alpha);

///  = A / alpha
template<typename T, typename U>
Tensor<T> operator/(Tensor<T> A, const U alpha) {
	return A /= (T) alpha;
}

template Tensor<f> operator/(Tensor<f> A, const f alpha);
template Tensor<d> operator/(Tensor<d> A, const d alpha);
template Tensor<cf> operator/(Tensor<cf> A, const f alpha);
template Tensor<cf> operator/(Tensor<cf> A, const cf alpha);
template Tensor<cd> operator/(Tensor<cd> A, const d alpha);
template Tensor<cd> operator/(Tensor<cd> A, const cd alpha);

template <typename T, typename U, template <typename> class Dev>
void cast(Tensor<T, Dev>& L, const Tensor<U, Dev>& R) {
    for (size_t I = 0; I < L.shape_.totalDimension(); ++I) {
        L(I) = (T) R(I);
    }
}

template void cast(Tensor<d>& L, const Tensor<f>& R);
template void cast(Tensor<f>& L, const Tensor<d>& R);
template void cast(Tensor<cd>& L, const Tensor<cf>& R);
template void cast(Tensor<cf>& L, const Tensor<cd>& R);

template<typename T>
Tensor<T> productElementwise(const Tensor<T>& A, const Tensor<T>& B) {
	assert(A.shape_.totalDimension() == B.shape_.totalDimension());
	Tensor<T> C(A.shape_);
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		C(i) = A(i) * B(i);
	}
	return C;
}

template <typename T, template <typename> class Dev>
void mdiagm(Tensor<T, Dev>& C, const Tensor<T, Dev>& B, const Tensor<T, Dev>& diag) {
    /**
     * @brief Multiply a dense matrix with a diagonal matrix
     * @param C output Matrix that is written on
     * @param diag diagonal Matrix
     * @param B dense input Matrix
     * 
     * C_ij = C_ij + B_ij * A_ii
     */
    if (B.shape_.order() != 2) {
        cerr << "Error: expected order of B to be two but received " << B.shape_.order() << "\n";
        exit(3);
    }
    if (diag.shape_.order() != 1) {
        cerr << "Error: expected order of diag to be one but received " << diag.shape_.order() << "\n";
        exit(3);
    }
    size_t n = diag.shape_[0];
    size_t inc_a = 1;
    size_t inc_b = 1;
    for (size_t i = 0; i < n; ++i) {
        T alpha = diag(i);
        const T* Bstart = B.data() + i * n;
        T* Cstart = C.data() + i * n;
	    blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b);
    }
}

template void mdiagm(Tensor<f>& C, const Tensor<f>&, const Tensor<f>&);
template void mdiagm(Tensor<d>& C, const Tensor<d>&, const Tensor<d>&);
template void mdiagm(Tensor<cf>& C, const Tensor<cf>&, const Tensor<cf>&);
template void mdiagm(Tensor<cd>& C, const Tensor<cd>&, const Tensor<cd>&);

template <typename T, template <typename> class Dev>
void diagmm(Tensor<T, Dev>& C, const Tensor<T, Dev>& diag, const Tensor<T, Dev>& B) {
    /**
     * @brief Multiply a dense matrix with a diagonal matrix
     * @param C output Matrix that is written on
     * @param B dense input Matrix
     * @param diag diagonal Matrix
     * 
     * C_ij = C_ij + A_ii * B_ij
     */
    if (B.shape_.order() != 2) {
        cerr << "Error: expected order of B to be two but received " << B.shape_.order() << "\n";
        exit(3);
    }
    if (diag.shape_.order() != 1) {
        cerr << "Error: expected order of diag to be one but received " << diag.shape_.order() << "\n";
        exit(3);
    }

    size_t n = diag.shape_[0];
    size_t inc_a = n;
    size_t inc_b = n;
    for (size_t i = 0; i < n; ++i) {
        T alpha = diag(i);
        const T* Bstart = B.data() + i;
        T* Cstart = C.data() + i;
	    blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b);
    }
}

template void diagmm(Tensor<f>& C, const Tensor<f>&, const Tensor<f>&);
template void diagmm(Tensor<d>& C, const Tensor<d>&, const Tensor<d>&);
template void diagmm(Tensor<cf>& C, const Tensor<cf>&, const Tensor<cf>&);
template void diagmm(Tensor<cd>& C, const Tensor<cd>&, const Tensor<cd>&);

template<typename T>
Tensor<T> conj(Tensor<T> A) {
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = conj(A(i));
	}
	return A;
}

template< >
Tensor<double> conj<double>(Tensor<double> A) {
	return A;
}

template< >
Tensor<float> conj<float>(Tensor<float> A) {
	return A;
}

template Tensor<cf> conj(Tensor<cf> A);
template Tensor<cd> conj(Tensor<cd> A);

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

template <typename T, typename U, template <typename> class Dev>
void offDiagonal(Tensor<T, Dev>& off, const Tensor<U, Dev>& full) {
    off = full;
    for (size_t i = 0; i < off.shape_[0]; ++i) {
        off(i, i) = 0;
    }
}

template void offDiagonal(Tensor<f>& off, const Tensor<f>& full);
template void offDiagonal(Tensor<d>& off, const Tensor<d>& full);
template void offDiagonal(Tensor<cf>& off, const Tensor<cf>& full);
template void offDiagonal(Tensor<cd>& off, const Tensor<cd>& full);

template<typename T>
T trace(const Tensor<T>& A) {
	Tensor<T> diag = diagonal(A);
	T tr = 0.;
	for (size_t i = 0; i < diag.shape_.totalDimension(); ++i) {
		tr += diag(i);
	}
	return tr;
}

/// dest[dim1, dim2] = beta * dest[dim1,dim2] + A[dim2, dim1]
template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta) {
	blas::scal(dim1 * dim2, beta, dest, 1);
	for (size_t j = 0; j < dim2; ++j) {
		const T* x = src + dim1 * j;
		size_t incx = 1;
		T* y = dest + j;
		size_t incy = dim2;
		blas::axpy(dim1, 1., x, incx, y, incy);
	}
}

template void transpose(f *dest, const f *src, size_t dim1, size_t dim2, f beta);
template void transpose(d *dest, const d *src, size_t dim1, size_t dim2, d beta);
template void transpose(cf *dest, const cf *src, size_t dim1, size_t dim2, cf beta);
template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2, cd beta);

/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
void transpose(Tensor<T>& dest, const Tensor<T>& src) {
	size_t dim1 = src.shape_.lastBefore();
	size_t dim2 = src.shape_.lastDimension();
	transpose(dest.data(), src.data(), dim1, dim2);
	dest.shape_ = transposeToFront(src.shape_, src.shape_.lastIdx());
}

template void transpose(Tensor<f>& dest, const Tensor<f>& src);
template void transpose(Tensor<d>& dest, const Tensor<d>& src);
template void transpose(Tensor<cf>& dest, const Tensor<cf>& src);
template void transpose(Tensor<cd>& dest, const Tensor<cd>& src);


/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
Tensor<T> transpose(const Tensor<T>& src) {
	Tensor<T> dest(src.shape_);
	transpose(dest, src);
	return dest;
}

template Tensor<f> transpose(const Tensor<f>& src);
template Tensor<d> transpose(const Tensor<d>& src);
template Tensor<cf> transpose(const Tensor<cf>& src);
template Tensor<cd> transpose(const Tensor<cd>& src);

/// \brief perform matrix adjoint A[bef, last] --> A*[last, bef]
template<typename T>
void adjoint(Tensor<T>& dest, const Tensor<T>& src) {
	transpose(dest, src);
	dest = conj(dest);
}

template void adjoint(Tensor<f>& dest, const Tensor<f>& src);
template void adjoint(Tensor<d>& dest, const Tensor<d>& src);
template void adjoint(Tensor<cf>& dest, const Tensor<cf>& src);

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

/// dest[a, c, b] = src[a, b, c]
template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C) {
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			const T* x = src + b * A + c * A * B;
			      T* y = dest + c * A + b * A * C;
			blas::copy(A, x, 1, y, 1);
		}
	}
}

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
	transposeBC(dest.data(), src.data(), a, b, c);

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

// https://stackoverflow.com/questions/50338955/how-can-i-concisely-write-a-lot-of-explicit-function-template-instantiations
template<typename... Ts>
auto instantiateTensorBLAS() {
    static auto funcs = std::tuple_cat(std::make_tuple(
//        nrm2<Ts>,
//		axpy<Ts>,
//		residual<Ts>,
		productElementwise<Ts>,
		trace<Ts>,
		diagonal<Ts>,
		transposeAB<Ts>,
		transposeBC<Ts>
    )...);

    return &funcs;
}

template auto instantiateTensorBLAS<f, d, cf, cd>();
