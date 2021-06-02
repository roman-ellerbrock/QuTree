//
// Created by Roman Ellerbrock on 5/22/21.
//

#ifndef TENSORBLAS_H
#define TENSORBLAS_H
#include "Tensor.h"

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2);

template<typename T, int blocksize>
void transpose2(T *dest, const T *src, size_t lda, size_t ldb);

template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C);

template <typename T, typename U>
void matrixTensor1(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);

template <typename T, typename U>
void matrixTensor2(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	Tensor<T>& D, size_t before, size_t active, size_t activeC, size_t after, bool zero);

template <typename T, typename U>
void matrixTensor3(Tensor<T>& C, const Matrix<U>& h, const Tensor<T>& B,
	Tensor<T>& Ket_work, Tensor<T>& hKet_work,
	size_t before, size_t active, size_t activeC, size_t after, bool zero);


template<typename T, typename U>
void contraction1(Matrix<U>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

template<typename T>
void contraction2(Matrix<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	Tensor<T>& bra_work, Tensor<T>& ket_work,
	size_t A, size_t B, size_t B2, size_t C, bool zero);

/// ==== Wrappers ====
template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T, typename U>
Tensor<T> matrixTensorBLAS(const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T>
void contractionBLAS(Matrix<T>& h, const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T>
Matrix<T> contractionBLAS(const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero = true);


#endif //TENSORBLAS_H
