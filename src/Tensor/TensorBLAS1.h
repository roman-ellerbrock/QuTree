//
// Created by Roman Ellerbrock on 11/7/21.
//

#ifndef TENSORBLAS1_H
#define TENSORBLAS1_H
#include "Tensor.h"
#include <blas.hh>

//////////////////////////////////////////////////////////
/// Level 1 Tensor BLAS (operations with O(n^d))
//////////////////////////////////////////////////////////
/**
 * Rationale:
 * Routines that operate on Tensors that have runtime
 * O(n^d) or less and where the TensorShape
 * is {n, ..., n} (d-th order).
 * Includes fundamental creation, copy, move, +,-,*,/, ||.||, etc.
 */

/**
 * \brief Calculate 2-norm of a Tensor
 * @param incr increment in loop
 * @return returns the 2-norm of a Tensor, i.e. ||A||_2
 */
template<typename T, template <typename> class Tensor, class ...Queue>
double nrm2(const Tensor<T>& A, size_t incr = 1, Queue& ...queue);

/**
 * \brief Perform vector addition
 * @param inc_a increment in loop for a
 * @param inc_b increment in loop for b
 * @return returns alpha*(a+b)
 */
template<typename T, class Tensor, class ...Queue>
void axpy(const Tensor& A, Tensor& B, T alpha = 1., size_t inc_a = 1, size_t inc_b = 1, Queue& ...queue);

template<typename T>
void operator+=(Tensor<T>& A, const Tensor<T>& B);

template<typename T>
void operator-=(Tensor<T>& A, const Tensor<T>& B);

template<typename T, template <typename> class Tensor, class ...Queue>
double residual(Tensor<T> A, const Tensor<T>& B, Queue& ...queue);

template<typename T>
Tensor<T> operator+(Tensor<T> A, const Tensor<T>& B);

template<typename T>
Tensor<T> operator-(Tensor<T> A, const Tensor<T>& B);

template<typename T, typename U>
Tensor<T>& operator*=(Tensor<T>& A, U a);

template<typename T, typename U>
Tensor<T> operator*(U alpha, Tensor<T> A);

template<typename T, typename U>
Tensor<T> operator*(Tensor<T> A, U alpha);

template<typename T, typename U>
Tensor<T>& operator/=(Tensor<T>& A, U alpha);

template<typename T, typename U>
Tensor<T> operator/(Tensor<T> A, U alpha);

template <typename T, typename U, template <typename> class Dev>
void cast(Tensor<T, Dev>& L, const Tensor<U, Dev>& R);

template<typename T, template <typename> class Dev = polymorphic::hostMemory>
void hadamardProduct(Tensor<T, Dev>& C, const Tensor<T, Dev>& A, const Tensor<T, Dev>& B);

template <typename T, template <typename> class Tensor, class ...Queue>
void mdiagm(Tensor<T>& C, const Tensor<T>& B, const Tensor<T>& diag, T factor = (T)1., Queue& ...queue);

template <typename T, template <typename> class Dev, class ...Queue>
void diagmm(Tensor<T, Dev>& C, const Tensor<T, Dev>& diag, const Tensor<T, Dev>& B, T factor = (T)1., Queue& ...queue);

template<typename T>
[[nodiscard]] Tensor<T> conj(Tensor<T> A);

template<typename T>
Tensor<T> diagonal(const Tensor<T>& A);

template <typename T, typename U, template <typename> class Dev>
void offDiagonal(Tensor<T, Dev>& off, const Tensor<U, Dev>& full);

template<typename T>
T trace(const Tensor<T>& A);

/// \brief perform matrix transpose A[bef, last] --> A[last, bef] (multiply with old beta
template<typename T>
void transpose(T *dest, const T *src, size_t dim1, size_t dim2, T beta = 0.);

/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
//template<typename T>
//void transpose(T *dest, const T *src, size_t dim1, size_t dim2);

/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
void transpose(Tensor<T>& dest, const Tensor<T>& src);

/// \brief perform matrix transpose A[bef, last] --> A[last, bef]
template<typename T>
Tensor<T> transpose(const Tensor<T>& src);

/// \brief perform matrix adjoint A[bef, last] --> A*[last, bef]
template<typename T>
void adjoint(Tensor<T>& dest, const Tensor<T>& src);

/// \brief perform matrix adjoint A[bef, last] --> A*[last, bef]
template<typename T>
Tensor<T> adjoint(const Tensor<T>& src);

/// \brief perform tensor transpose dest[b, a, c] <-- src[a, b, c]
template<typename T>
void transposeAB(T *dest, const T *src, size_t A, size_t B, size_t C);

/// \brief perform tensor transpose dest[a, c, b] <-- src[a, b, c]
template<typename T>
void transposeBC(T *dest, const T *src, size_t A, size_t B, size_t C);

/// \brief perform tensor transpose dest[...,...,k] = src[...,k,...] if back == false
template<typename T>
void transpose(Tensor<T>& dest, const Tensor<T>& src, size_t k, bool back = false);

template<typename T>
Tensor<T> transpose(const Tensor<T>& src, size_t k, bool back = false);

/// \brief checks whether Tensor A is an identity Tensor
template<typename T>
double isCloseToIdentity(const Tensor<T>& A);

size_t nrows(const TensorShape& shape, blas::Op op = blas::Op::NoTrans);
size_t ncols(const TensorShape& shape, blas::Op op = blas::Op::NoTrans);

// \brief B(I(k), i_k) *= a(i_k)
template<typename T, typename U>
void vectorTensor(Tensor<T>& B, const Tensor<U>& a, size_t k);

#endif //TENSORBLAS1_H
