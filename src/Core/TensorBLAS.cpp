//
// Created by Roman Ellerbrock on 5/22/21.
//

#include "Core/TensorBLAS.h"
#include <cblas.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2) {
	// A[dim1, dim2] --> A[dim2, dim1]
	/// simple in-place transpose
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
			dest[j + dim2 * i] = src[i + dim1 * j];
		}
	}
}

template<typename T, int blocksize>
void transpose2(T *dest, const T *src, size_t lda, size_t ldb) {
	/// Incorporate blocking to minimize cache misses
	// dest[b, a, c] = src[a, b, c]
	/// simple in-place transpose
	cerr << "check for rectangular shapes.\n";
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
						h[act * before + actC] * B[dimafter * aft + act * before + bef];
				}
			}
		}
	}
}

template<typename T, typename U>
void matrixTensor2(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	Tensorcd& D, size_t before, size_t active, size_t activeC, size_t after, bool zero) {

	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }

	if (before == 1) {
		size_t m = activeC;
		size_t k = active; //activeB
		size_t n = after;
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
			(void *) &z, (void *) &h[0], m, (void *) &B[0], k, (void *) &zz, (void *) &C[0], m);
		return;
	}

	size_t m = activeC;
	size_t k = active; //activeB
	size_t n = before;

	size_t pref = before * active;
	size_t prefC = before * activeC;
	for (size_t aft = 0; aft < after; ++aft) {
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
			(void *) &z, (void *) &h[0], m, (void *) &B[aft * pref], n, (void *) &zz, (void *) &D[aft * prefC], m);
//		transpose2<T, 4>(&C[aft * prefC], &D[aft * prefC], activeC, before);
		transpose(&C[aft * prefC], &D[aft * prefC], activeC, before);
	}
}

template<typename T, typename U>
void matrixTensor3(Tensor<T>& hKet, const Matrix<U>& h, const Tensor<T>& Ket,
	Tensor<T>& Ket_work, Tensor<T>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	/// hKet[a, b2, c] = h[b2, b] * Ket[a, b, c]
	/// A, B, B2, C are dimensions of indices a, b, b2, c
	/// Transpose tensor and map to BLAS ?geem

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
template<typename T, typename U>
void contraction1(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {

	if (zero) { h.zero(); }

	for (size_t a = 0; a < A; ++a) {
		for (size_t b = 0; b < B; ++b) {
			for (size_t b2 = 0; b2 < B2; ++b2) {
				for (size_t c = 0; c < C; ++c) {
//					h(b, b2) += conj(bra(a, b, c)) * ket(a, b2, c);
					h(b, b2) += conj(bra[a + b * A + c * A * B]) * ket(a + b2 * A + c * A * B2);
				}
			}
		}
	}
}

template<typename T, typename U>
void contraction2(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	Tensorcd& bra_work, Tensorcd& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {

//	Tensorcd bra_cpy(bra.shape());
//	Tensorcd ket_work(bra.shape());
	transposeBC(&bra_work[0], &bra[0], A, B, C);
	transposeBC(&ket_work[0], &ket[0], A, B2, C);

	size_t AC = A * C;
	T z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
		(void *) &z, (void *) &bra_work[0], AC, (void *) &ket_work[0], AC, (void *) &zz, (void *) &h[0], B);
}

/// ========================================================================
/// Template Instantiations
/// ========================================================================

typedef complex<double> cd;

typedef double d;

template void matrixTensor1(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor2(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	Tensorcd& D,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor3(Tensor<cd>& hKet, const Matrix<cd>& h, const Tensor<cd>& Ket,
	Tensor<cd>& Ket_work, Tensor<cd>& hKet_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void contraction1(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void contraction2(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	Tensorcd& bra_work, Tensorcd& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void transpose(cd *dest, const cd *src, size_t dim1, size_t dim2);

template void transpose2<complex<double>,4>(cd *dest, const cd *src, size_t lda, size_t ldb);

template void transposeAB(cd *dest, const cd *src, size_t A, size_t B, size_t C);
