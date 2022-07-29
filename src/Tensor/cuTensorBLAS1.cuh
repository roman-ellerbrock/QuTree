#ifndef CUTENSORBLAS1_CUH
#define CUTENSORBLAS1_CUH
#include <stdafx.h>

template <typename TL, typename TR>
void cudaCast(TL* dest, const TR* src, size_t n);

template <typename T>
void cudaHadamardProduct(T* C, const T* A, const T* B, size_t n);

template <typename T>
void cudaDiagMatrixMatrix(T* C, const T* dA, const T* B, T factor, size_t nrow, size_t ncol);

#endif // CUTENSORBLAS1_CUH