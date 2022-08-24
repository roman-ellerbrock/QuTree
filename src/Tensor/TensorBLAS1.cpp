//
// Created by Roman Ellerbrock on 11/7/21.
//

#include "Tensor/TensorBLAS1.hpp"
#include "Tensor/Tensor.hpp"
#include "stdafx.h"
#include <math.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

using f = float;
using d = double;
using cf = complex<f>;
using cd = complex<d>;

template f nrm2<f>(const Tensorf& A, size_t incr);
template d nrm2<d>(const Tensord& A, size_t incr);
template f nrm2<cf>(const Tensorcf& A, size_t incr);
template d nrm2<cd>(const Tensorcd& A, size_t incr);

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
	axpy<T>(B, A, alpha);
}

template void operator-=(Tensor<f>& A, const Tensor<f>& B);
template void operator-=(Tensor<d>& A, const Tensor<d>& B);
template void operator-=(Tensor<cf>& A, const Tensor<cf>& B);
template void operator-=(Tensor<cd>& A, const Tensor<cd>& B);

template<typename T>
double residual(Tensor<T> A, const Tensor<T>& B) {
	A -= B;
	return (double) abs(nrm2<T>(A));
}

/// || A - B ||_2
template double residual<f>(Tensorf A, const Tensorf& B);
template double residual<d>(Tensord A, const Tensord& B);
template double residual<cf>(Tensorcf A, const Tensorcf& B);
template double residual<cd>(Tensorcd A, const Tensorcd& B);

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
	if (alpha == (U) 1.) { return A; }
	if (alpha == (U) 0.) { A.zero(); return A; }
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

template<typename T, template <typename> class Dev>
void hadamardProduct(Tensor<T, Dev>& C, const Tensor<T, Dev>& A, const Tensor<T, Dev>& B) {
	assert(A.shape_.totalDimension() == B.shape_.totalDimension());
	for (size_t i = 0; i < A.shape_.totalDimension(); i++) {
		C(i) = A(i) * B(i);
	}
}

template void hadamardProduct<f>(Tensorf&, const Tensorf& A, const Tensorf& B);
template void hadamardProduct<d>(Tensord&, const Tensord& A, const Tensord& B);
template void hadamardProduct<cf>(Tensorcf&, const Tensorcf& A, const Tensorcf& B);
template void hadamardProduct<cd>(Tensorcd&, const Tensorcd& A, const Tensorcd& B);

template <typename T, template <typename> class Dev, class ...Queue>
void mdiagm(Tensor<T, Dev>& C, const Tensor<T, Dev>& B, const Tensor<T, Dev>& diag, T factor, Queue& ...queue) {
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
		T alpha = diag(i) * factor;
		const T* Bstart = B.data() + i * n;
		T* Cstart = C.data() + i * n;
		blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b, queue...);
    }
}

template void mdiagm(Tensor<f>& C, const Tensor<f>&, const Tensor<f>&, f);
template void mdiagm(Tensor<d>& C, const Tensor<d>&, const Tensor<d>&, d);
template void mdiagm(Tensor<cf>& C, const Tensor<cf>&, const Tensor<cf>&, cf);
template void mdiagm(Tensor<cd>& C, const Tensor<cd>&, const Tensor<cd>&, cd);

template <typename T, template <typename> class Dev, class ...Queue>
void diagmm(Tensor<T, Dev>& C, const Tensor<T, Dev>& diag, const Tensor<T, Dev>& B, T factor, Queue& ...queue) {
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
        T alpha = diag(i) * factor;
        const T* Bstart = B.data() + i;
        T* Cstart = C.data() + i;
	    blas::axpy(n, alpha, Bstart, inc_a, Cstart, inc_b, queue...);
    }
}

template void diagmm(Tensor<f>& C, const Tensor<f>&, const Tensor<f>&, f);
template void diagmm(Tensor<d>& C, const Tensor<d>&, const Tensor<d>&, d);
template void diagmm(Tensor<cf>& C, const Tensor<cf>&, const Tensor<cf>&, cf);
template void diagmm(Tensor<cd>& C, const Tensor<cd>&, const Tensor<cd>&, cd);

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

template void diagonal(Tensorf&, const Tensorf&);
template void diagonal(Tensord&, const Tensord&);
template void diagonal(Tensorcf&, const Tensorcf&);
template void diagonal(Tensorcd&, const Tensorcd&);

template Tensorf diagonal(const Tensorf&);
template Tensord diagonal(const Tensord&);
template Tensorcf diagonal(const Tensorcf&);
template Tensorcd diagonal(const Tensorcd&);

template void addDiagonal(Tensorf& B, const Tensorf& diag, f alpha);
template void addDiagonal(Tensord& B, const Tensord& diag, d alpha);
template void addDiagonal(Tensorcf& B, const Tensorcf& diag, cf alpha);
template void addDiagonal(Tensorcd& B, const Tensorcd& diag, cd alpha);

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
	if (beta != (T)1) {blas::scal(dim1 * dim2, beta, dest, 1);}
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
			dest[i * dim2 + j] += src[j * dim1 + i];
		}
/*		const T* x = src + dim1 * j;
		size_t incx = 1;
		T* y = dest + j;
		size_t incy = dim2;
		blas::axpy(dim1, 1., x, incx, y, incy);
		*/
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
			memcpy(&dest[c * A + b * A * C], &src[b * A + c * A * B], A * sizeof(T));
//			const T* x = src + b * A + c * A * B;
//			      T* y = dest + c * A + b * A * C;
//			blas::copy(A, x, 1, y, 1);
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
//		hadamardProduct<Ts>,
		trace<Ts>,
//		diagonal<Ts>,
		transposeAB<Ts>,
		transposeBC<Ts>
    )...);

    return &funcs;
}

template auto instantiateTensorBLAS<f, d, cf, cd>();
