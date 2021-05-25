//
// Created by Roman Ellerbrock on 5/22/21.
//

#include "Core/TensorMath.h"
#include <cblas.h>

#define swap(type, x, y) { type _tmp; _tmp = x; x = y; y = _tmp; }

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2) {
	// A[dim1, dim2] --> A[dim2, dim1]
	/// simple in-place transpose
	for (size_t j = 0; j < dim2; ++j) {
		for (size_t i = 0; i < dim1; ++i) {
/*			T tmp = A[dim1 * j + i];
			A[dim1 * j + i] = A[dim2 * i + j];
			A[dim2 * i + j] = tmp;*/
//			swap(T, A[j + dim2 * i], A[i + dim1 * j]);
			dest[j + dim2 * i] = src[i + dim1 * j];
		}
	}
}

template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[b, a, c] = src[a, b, c]
	/// simple in-place transpose
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			for (size_t a = 0; a < A; ++a) {
				dest[b + a * B + c * A * B] = src[a + b * A + c * A * B];
			}
		}
	}
}

template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C) {
	// dest[a, c, b] = src[a, b, c]
	/// simple in-place transpose
	for (size_t c = 0; c < C; ++c) {
		for (size_t b = 0; b < B; ++b) {
			for (size_t a = 0; a < A; ++a) {
				dest[a + c * A + b * A * C] = src[a + b * A + c * A * B];
			}
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
	size_t before, size_t active, size_t activeC, size_t after, bool zero) {

	if (zero) { C.zero(); }

	//  C[be, acC, af] += h[acC, acB] * B[be, acB, af]
	if (before == 1) {
		//  C[acC, af] += h[acC, acB] * B[acB, af]
		size_t m = activeC;
		size_t k = active; //activeB
		size_t n = after;
		complex<double> z = 1.0;
		complex<double> zz = 0.0;
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, k,
			(void *) &z, (void *) &h[0], m, (void *) &B[0], k, (void *) &zz, (void *) &C[0], m);
		return;
	}

	//  C[be, acC] += h[acC, acB] * B[be, acB]
	//  BLAS: C[acC, be] += h[acC, acB] * B[be, acB]
	//  transposed: C[be, acC]
	size_t m = activeC;
	size_t k = active; //activeB
	size_t n = before;
	size_t pref = before * active;
	size_t prefC = before * activeC;
	complex<double> z = 1.0;
	complex<double> zz = 0.0;
	Tensorcd D = C;
	for (size_t aft = 0; aft < after; ++aft) {
		cblas_zgemm(CblasColMajor, CblasNoTrans, CblasTrans, m, n, k,
			(void *) &z, (void *) &h[0], m, (void *) &B[aft * pref], n, (void *) &zz, (void *) &D[aft * prefC], m);
		transpose(&C[aft * prefC], &D[aft * prefC], activeC, before);
	}
}

template<typename T, typename U>
void matrixTensor3(Tensor<T>& hKet, const Matrix<U>& h, const Tensor<T>& Ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero) {
	/// hKet[a, b2, c] = h[b2, b] * Ket[a, b, c]
	/// A, B, B2, C are dimensions of indices a, b, b2, c
	/// Transpose tensor and map to BLAS ?geem

	Tensorcd Ket_cpy(Ket.shape());
	Tensorcd hKet_cpy = hKet;

	transposeAB(&Ket_cpy[0], &Ket[0], A, B, C);
	size_t AC = A * C;
	complex<double> z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, B2, AC, B,
		(void *) &z, (void *) &h[0], B2, (void *) &Ket_cpy[0], B, (void *) &zz, (void *) &hKet_cpy[0], B2);

	transposeAB(&hKet[0], &hKet_cpy[0], B2, A, C);
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
	size_t A, size_t B, size_t B2, size_t C, bool zero) {

	Tensorcd bra_cpy(bra.shape());
	transposeBC(&bra_cpy[0], &bra[0], A, B, C);
	Tensorcd ket_cpy(bra.shape());
	transposeBC(&ket_cpy[0], &ket[0], A, B2, C);

	size_t AC = A * C;
	complex<double> z = 1.0, zz = 1.;
	if (zero) { zz = 0.; }
	cblas_zgemm(CblasColMajor, CblasConjTrans, CblasNoTrans, B, B2, AC,
		(void *) &z, (void *) &bra_cpy[0], AC, (void *) &ket_cpy[0], AC, (void *) &zz, (void *) &h[0], B);
}

typedef complex<double> cd;

typedef double d;

template void matrixTensor1(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor2(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void matrixTensor3(Tensor<cd>& C, const Matrix<cd>& h, const Tensor<cd>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template void contraction1(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template void contraction2(Matrix<cd>& h, const Tensor<cd>& bra, const Tensor<cd>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero);
