//
// Created by Roman Ellerbrock on 5/22/21.
//

#include "Core/TensorBLAS.h"
#include <cblas.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta) {
	// A[dim1, dim2] --> A[dim2, dim1]
	/// simple in-place transpose
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
//			dest[j + dim2 * i] = src[i + dim1 * j];
			dest[j + dim2 * i] = beta * dest[j + dim2 * i] + src[i + dim1 * j];
		}
	}
}

template<typename T, int blocksize>
void transpose2(T *dest, const T *src, size_t lda, size_t ldb) {
	/// Incorporate blocking to minimize cache misses
	// dest[b, a, c] = src[a, b, c]
	/// simple in-place transpose
	size_t nblockA = lda / blocksize;
	size_t nblockB = ldb / blocksize;
	for (size_t b = 0; b < nblockB; ++b) {
		for (size_t a = 0; a < nblockA; ++a) {
			for (size_t bb = 0; bb < blocksize; ++bb) {
				for (size_t aa = 0; aa < blocksize; ++aa) {
					// dest[a', b'] <-- src[b', a']
					// a' = a * blocksize + aa
					// b' = b * blocksize + bb
					// idxDest = b' * lda + a' = (b * blocksize + bb) * lda + a * blocksize + aa
					// idxSrc  = a' * ldb + b' = (a * blocksize + aa) * ldb + b * blocksize + bb
					size_t idxDest = (b * blocksize + bb) * lda + a * blocksize + aa;
					size_t idxSrc = (a * blocksize + aa) * ldb + b * blocksize + bb;
					dest[idxDest] = src[idxSrc];
				}
			}
		}
	}

	// right block and square block
	size_t restA = lda % blocksize;
	size_t offsetA = lda - restA;
	for (size_t a = 0; a < restA; ++a) {
		for (size_t b = 0; b < ldb; ++b) {
			// dest[a', b'] <-- src[a', b']
			size_t aa = offsetA + a;
			dest[aa + b * lda] = src[aa * ldb + b];
		}
	}
	// bottom block
	size_t restB = lda % blocksize;
	size_t offsetB = ldb - restB;
	for (size_t b = 0; b < restB; ++b) {
		for (size_t a = 0; a < nblockA * blocksize; ++a) {
			// dest[a', b'] <-- src[a', b']
			size_t bb = offsetB + b;
			dest[a + bb * lda] = src[a * ldb + bb];
		}
	}
}

template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[b, a, c] = src[a, b, c]
	for (size_t c = 0; c < C; ++c) {
		//	transpose2<T,4>(&dest[c * A * B], &src[c * A * B], A, B);
		transpose(&dest[c * A * B], &src[c * A * B], A, B);
	}
}

template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[a, c, b] = src[a, b, c]
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			memcpy(&dest[c * A + b * A * C], &src[b * A + c * A * B], A * sizeof(T));
/*			for (size_t a = 0; a < A; ++a) {
				dest[a + c * A + b * A * C] = src[a + b * A + c * A * B];
			}*/
		}
	}
}

template<typename T, typename U>
void matrixTensor1(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero) {

	if (zero) { C.zero(); }

	size_t dimafter = active * before;
	size_t dimafterC = activeC * before;
	for (size_t aft = 0; aft < after; ++aft) {
		for (size_t act = 0; act < active; ++act) {
			for (size_t actC = 0; actC < activeC; ++actC) {
				for (size_t bef = 0; bef < before; ++bef) {
					C[dimafterC * aft + actC * before + bef] +=
						h[act * activeC + actC] * B[dimafter * aft + act * before + bef];
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor2(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	Tensor<T>& D, size_t before, size_t active, size_t activeC, size_t after, bool zero) {
	/// D is work tensor with shape of C.
	typedef complex<double> cd;
	typedef double d;
	/// If not double/complex<double> type, fall back to straightforward implementation
	if constexpr(!((is_same<U, cd>::value && is_same<T, cd>::value))
		|| (is_same<U, d>::value && is_same<T, d>::value)) {
		matrixTensor1(C, h, B, before, active, activeC, after, zero);
		return;
	}

	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }

/*	if (before == 1) {
		size_t m = activeC;
		size_t k = active; //activeB
		size_t n = after;
		if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
				(void *) &z, (void *) &h[0], m, (void *) &B[0], k, (void *) &zz, (void *) &C[0], m);
		} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
				z, (double *) &h[0], m, (double *) &B[0], k, zz, (double *) &C[0], m);
		}
		return;
	}*/

	size_t m = activeC;
	size_t k = active; //activeB
	size_t n = before;

	size_t pref = before * active;
	size_t prefC = before * activeC;
	T zer = 0.;
	for (size_t aft = 0; aft < after; ++aft) {
		if constexpr(is_same<U, cd>::value && is_same<T, cd>::value) {
			cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
				(void *) &z, (void *) &h[0], m, (void *) &B[aft * pref], n, (void *) &zer, (void *) &D[aft * prefC], m);
		} else if constexpr(is_same<U, d>::value && is_same<T, d>::value) {
			cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
				z, (double *) &h[0], m, (double *) &B[aft * pref], n, zer, (double *) &D[aft * prefC], m);
		}
//		transpose2<T, 4>(&C[aft * prefC], &D[aft * prefC], activeC, before);
		transpose(&C[aft * prefC], &D[aft * prefC], activeC, before, zz);
	}
	if (!zero) { }
}

template<typename T, typename U>
void matrixTensor3(Tensor<T>& hKet, const Matrix<U>& h, const Tensor<T>& Ket,
	Tensor<T>& Ket_work, Tensor<T>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	/// hKet[a, b2, c] = h[b2, b] * Ket[a, b, c]
	/// A, B, B2, C are dimensions of indices a, b, b2, c
	/// Transpose tensor and map to BLAS ?geem
	typedef complex<double> cd;
	typedef double d;

	size_t AC = A * C;
	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	transposeAB(&Ket_work[0], &Ket[0], A, B, C);
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, B2, AC, B,
		(void *) &z, (void *) &h[0], B2, (void *) &Ket_work[0], B,
		(void *) &zz, (void *) &hKet_work[0], B2);

	transposeAB(&hKet[0], &hKet_work[0], B2, A, C);
}

//====================================================================

//====================================================================

template<typename T>
void contraction2(Matrix<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	Tensor<T>& bra_work, Tensor<T>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	typedef complex<double> cd;
	typedef double d;
	/// If not double/complex<double> type, fall back to straightforward implementation
	if constexpr(!(is_same<T, cd>::value || is_same<T, d>::value)) {
		contraction1(h, bra, ket, A, B, B2, C, zero);
		return;
	}

	transposeBC(&bra_work[0], &bra[0], A, B, C);
	transposeBC(&ket_work[0], &ket[0], A, B2, C);

	size_t AC = A * C;
	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	if constexpr(is_same<T, cd>::value) {
		cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
			(void *) &z, (void *) &bra_work[0], AC, (void *) &ket_work[0], AC, (void *) &zz, (void *) &h[0], B);
	} else if constexpr(is_same<T, d>::value) {
		cblas_dgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
			z, (double *) &bra_work[0], AC, (double *) &ket_work[0], AC, zz, (double *) &h[0], B);
	}
}

/// ========================================================================
/// Wrappers for matrixTensor product and Tensor contraction
/// ========================================================================

/// Wrapper for matrix-Tensor product
template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, Tensor<T>& workC, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {

	TensorShape tdim(B.shape());
	TensorShape tdimC(C.shape());

	if (mode >= tdim.order()) {
		cerr << "matrixTensor error: mode too large.\n";
		exit(1);
	}

	size_t after = tdim.after(mode);
	size_t before = tdim.before(mode);
	size_t active1 = A.dim1();
	size_t active2 = A.dim2();

	if (!(A.dim2() == tdim[mode])) {
		cerr << "matrix Tensor error: active dimension wrong.\n";
		exit(1);
	}
	if (!(A.dim1() == tdimC[mode])) {
		cerr << "matrix Tensor error: left active dimension wrong.\n";
		exit(1);
	}
	if (tdim.before(mode) != tdimC.before(mode)) {
		cerr << "matrix Tensor error: before dimension wrong.\n";
		exit(1);
	}
	if (tdim.after(mode) != tdimC.after(mode)) {
		cerr << "matrix Tensor error: after dimension wrong.\n";
		exit(1);
	}
	if (tdimC.totalDimension() != workC.shape().totalDimension()) {
		cerr << "matrix Tensor error: work array has wrong dimension.\n";
		exit(1);
	}

	matrixTensor2(C, A, B, workC, before, active1, active2, after, zero);
}

template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	Tensor<T> workC(C);
	matrixTensorBLAS(C, workC, A, B, mode, zero);
}

template<typename T>
void contractionBLAS(Matrix<T>& h, Tensor<T>& workA, Tensor<T>& workB, const Tensor<T>& A, const Tensor<T>& B,
	size_t mode, bool zero) {

	TensorShape tdimA(A.shape());
	TensorShape tdimB(B.shape());
	if (mode >= tdimA.order()) {
		cerr << "contraction error: mode too large for Bra.\n";
		exit(1);
	}
	if (mode >= tdimB.order()) {
		cerr << "contraction error: mode too large for Ket.\n";
		exit(1);
	}

	size_t after = tdimA.after(mode);
	size_t before = tdimA.before(mode);
	size_t activeA = tdimA[mode];
	size_t activeB = tdimB[mode];

	if (!(h.dim1() == activeA)) {
		cerr << "contraction error: ket active dimension wrong.\n";
		exit(1);
	}
	if (!(h.dim2() == activeB)) {
		cerr << "contraction error: bra active dimension wrong.\n";
		exit(1);
	}
	if (tdimA.before(mode) != tdimB.before(mode)) {
		cerr << "contraction error: before dimension wrong.\n";
		exit(1);
	}
	if (tdimA.after(mode) != tdimB.after(mode)) {
		cerr << "contraction error: after dimension wrong.\n";
		exit(1);
	}

	if (tdimA.totalDimension() != workA.shape().totalDimension()) {
		cerr << "contraction error: work array A has wrong dimension.\n";
		exit(1);
	}
	if (tdimB.totalDimension() != workB.shape().totalDimension()) {
		cerr << "contraction error: work array B has wrong dimension.\n";
		exit(1);
	}

	contraction2(h, A, B, workA, workB, before, activeA, activeB, after, zero);
}

/// Wrapper for matrix-Tensor product
template<typename T>
void contractionBLAS(Matrix<T>& h, const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero) {
	Tensor<T> Awork(A.shape());
	Tensor<T> Bwork(B.shape());
	contractionBLAS(h, Awork, Bwork, A, B, mode, zero);
}

template<typename T>
Matrix<T> contractionBLAS(const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape tdimA(A.shape());
	TensorShape tdimB(B.shape());
	if (mode >= tdimA.order()) {
		cerr << "contraction error: mode too large for Bra.\n";
		exit(1);
	}
	if (mode >= tdimB.order()) {
		cerr << "contraction error: mode too large for Ket.\n";
		exit(1);
	}

	size_t dim1 = tdimA[mode];
	size_t dim2 = tdimB[mode];
	Matrix<T> h(dim1, dim2);
	contractionBLAS(h, A, B, mode, zero);
	return h;
}

template<typename T, typename U>
Tensor<T> matrixTensorBLAS(const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero) {
	TensorShape Cshape = B.shape();
	Cshape.setDimension(A.dim1(), mode);
	Tensor<T> C(Cshape);
	matrixTensorBLAS(C, A, B, mode, zero);
	return C;
}

/// ========================================================================
/// Template Instantiations
/// ========================================================================

typedef complex<double> cd;

typedef double d;

template void matrixTensor1(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor2(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	Tensorcd& D, size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor3(Tensor<cd>& hKet, const Matrix<cd>& h, const Tensor<cd>& Ket,
	Tensor<cd>& Ket_work, Tensor<cd>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void contraction2(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	Tensor<cd>& bra_work, Tensor<cd>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2, cd beta);

template void transpose2<complex<double>, 4>(cd *dest, const cd *src, size_t lda, size_t ldb);

template void transposeAB(cd *dest, const cd *src, size_t A, size_t B, size_t C);

/// ==== Wrappers ====
// complex double
template void matrixTensorBLAS(Tensor<cd>& C, Tensor<cd>& workC, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void matrixTensorBLAS(Tensor<cd>& C, const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<cd>& h, Tensor<cd>& workA, Tensor<cd>& workB, const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<cd>& h, const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Tensor<cd> matrixTensorBLAS(const Matrix<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
template Matrix<cd> contractionBLAS(const Tensor<cd>& A, const Tensor<cd>& B, size_t mode, bool zero);
// double
template void matrixTensorBLAS(Tensor<d>& C, Tensor<d>& workC, const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void matrixTensorBLAS(Tensor<d>& C, const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<d>& h, Tensor<d>& workA, Tensor<d>& workB, const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template void contractionBLAS(Matrix<d>& h, const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template Tensor<d> matrixTensorBLAS(const Matrix<d>& A, const Tensor<d>& B, size_t mode, bool zero);
template Matrix<d> contractionBLAS(const Tensor<d>& A, const Tensor<d>& B, size_t mode, bool zero);
