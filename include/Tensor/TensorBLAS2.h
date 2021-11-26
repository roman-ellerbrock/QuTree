//
// Created by Roman Ellerbrock on 11/7/21.
//

#ifndef TENSORBLAS2_H
#define TENSORBLAS2_H
#include "Tensor/TensorBLAS1.h"
#include <lapack.hh>


/// ===================== Matrix/Tensor Products ====================

/**
 * \brief c = beta * c + alpha * op_a(a) * op_b(b)
 * default: c = a*b
 * @tparam T
 * @param c
 * @param a
 * @param b
 * @param alpha
 * @param beta
 * @param op_a
 * @param op_b
 */
template<typename T>
void gemm(Tensor<T>& c, const Tensor<T>& a, const Tensor<T>& b,
	T alpha = 1., T beta = 0., blas::Op op_a = blas::Op::NoTrans, blas::Op op_b = blas::Op::NoTrans);

template<typename T>
Tensor<T> gemm(const Tensor<T>& a, const Tensor<T>& b, T alpha = 1.,
	blas::Op op_a = blas::Op::NoTrans, blas::Op op_b = blas::Op::NoTrans);

/// perform u^cT a u
template<typename T>
[[nodiscard]] Tensor<T> unitarySimilarityTrafo(Tensor<T> a, const Tensor<T>& u);

/**
  * \brief Perform matrix-tensor product on k-th index.
  *
  * Performs the operation C = alpha * h(k)^op * A + beta * C.
  *
  * @param C Output tensor
  * @param h input matrix
  * @param B input tensor
  * @param k1 target index in tensor B (and C)
  * @param alpha factor for h*B
  * @param beta factor for previous C
  * @param op_h ConjTrans/Trans/NoTrans on h
  */
template<typename T>
void matrixTensor(Tensor<T>& C, const Tensor<T>& h, const Tensor<T>& B,
	size_t k, T alpha = 1., T beta = 0., blas::Op op_h = blas::Op::NoTrans);

/**
  * \brief Perform matrix-tensor product on k-th index.
  *
  * Performs the operation C = alpha * h(k)^op * A + beta * C.
  *
  * @param C Output tensor
  * @param h input matrix
  * @param B input tensor
  * @param k1 target index in tensor B (and C)
  * @param alpha factor for h*B
  * @param beta factor for previous C
  * @param op_h ConjTrans/Trans/NoTrans on h
  */
template <typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& B,
	size_t k, T alpha = 1., T beta = 0., blas::Op op_h = blas::Op::NoTrans);

/**
 * \brief Perform tensor-hole contraction leaving out k-th index.
 *
 * Performs the operations h = alpha * tr(bra, ket)_(k) + beta * h
 *
 * @param h resulting matrix
 * @param bra input tensor <Bra|
 * @param ket output tensor |Ket>
 * @param k hole-index left out in contraction
 * @param alpha multiply contracted result by alpha
 * @param beta multiply previous h by beta
 */
template<typename T>
void contraction(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha = 1., T beta = 0.);

/**
 * \brief Perform tensor-hole contraction leaving out k-th index.
 *
 * Performs the operations h = alpha * tr(bra, ket)_(k) + beta * h
 *
 * @param h resulting matrix
 * @param bra input tensor <Bra|
 * @param ket output tensor |Ket>
 * @param k hole-index left out in contraction
 * @param alpha multiply contracted result by alpha
 * @param beta multiply previous h by beta
 */
template<typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha = 1.);


#endif //TENSORBLAS2_H
