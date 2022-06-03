//
// Created by Roman Ellerbrock on 5/22/21.
//

#ifndef TENSORBLAS_H
#define TENSORBLAS_H
#include "Tensor.h"

void dgeem(Matrixd& h, const Matrixd& bra, const Matrixd& ket);

template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta = 0.);

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

template<typename T>
void general_contraction(const Tensor<T>& A,
                         const Tensor<T>& B,
                         Tensor<T>& result,
                         const vector<size_t> &A_indices,
                         const vector<size_t> &B_indices);

template<typename T, typename Q>
void general_contraction(const Tensor<T>& A,
                         const Tensor<T>& B,
                         Tensor<T>& result,
                         const std::vector<std::pair<Q,Q>>& contraction_pairs);

template<typename T>
void general_transpose_bd(T* dst, const T* src, size_t a, size_t b, size_t c, size_t d, size_t e);

template<typename T>
bool is_contraction_legal(const Tensor<T> &TensorA,
                          const Tensor<T> &TensorB,
                          vector<size_t> Acontraction,
                          vector<size_t> Bcontraction);

template<typename T>
void general_transpose(Tensor<T>& dst, const Tensor<T>& src, size_t index_one, size_t index_two);

template<typename T>
void general_transpose_to_order(Tensor<T>& dst, const Tensor<T>& src, const vector<size_t> &form);

/// ==== Wrappers ====
template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, Tensor<T>& workC, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T, typename U>
void matrixTensorBLAS(Tensor<T>& C, const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T, typename U>
Tensor<T> matrixTensorBLAS(const Matrix<U>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T>
void contractionBLAS(Matrix<T>& h, Tensor<T>& workA, Tensor<T>& workB, const Tensor<T>& A,
	const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T>
void contractionBLAS(Matrix<T>& h, const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero = true);

template<typename T>
Matrix<T> contractionBLAS(const Tensor<T>& A, const Tensor<T>& B, size_t mode, bool zero = true);


#endif //TENSORBLAS_H
