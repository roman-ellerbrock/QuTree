//
// Created by Roman Ellerbrock on 11/7/21.
//

#ifndef TENSORBLAS2_H
#define TENSORBLAS2_H
#include "Tensor/TensorBLAS1.h"
#include <lapack.hh>

template<typename T>
void gemm(Tensor<T>& c, const Tensor<T>& a, const Tensor<T>& b,
	T alpha, T beta,
	blas::Op op_a, blas::Op op_b);

template<typename T>
Tensor<T> gemm(const Tensor<T>& a, const Tensor<T>& b,
	T alpha, blas::Op op_a, blas::Op op_b);

template<typename T>
void gemmRef(Tensor<T>& c, Tensor<T> a, Tensor<T> b,
	T alpha, T beta, blas::Op op_a, blas::Op op_b);

template<typename T, typename U>
void matrixTensorRef(Tensor<T>& C, const Tensor<U>& h, const Tensor<T>& B,
	size_t k, bool zero = true);

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
	size_t k, T alpha, T beta, blas::Op op_h);

template<typename T>
void contractionRef(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
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



/**
 * \brief calculates the QR-decomposition of a Tensor A.
 *
 * Treats A (A(i_1,i_2,...,i_d)) as a flattened matrix A_(Ji_d) and calculates A=Q*R
 * @param Q
 * @param A
 */
template<typename T>
void qr(Tensor<T>& Q, const Tensor<T>& A);

template<typename T>
Tensor<T> qr(const Tensor<T>& A);

template<typename T>
Tensor<T> qr(Tensor<T> A, size_t k);

template<typename T>
void gramSchmidt(Tensor<T>& A);

template<typename T>
void gramSchmidt(Tensor<T>& A, size_t k);

/// SVD
template <typename T>
using SVD = tuple<Tensor<T>, Tensor<T>, Tensord>;
typedef SVD<complex<double>> SVDcd;
typedef SVD<double> SVDd;

template<typename T>
void svd(Tensor<T>& AthenU, Tensor<T>& V, Tensord& sigma);

template<typename T>
SVD<T> svd(Tensor<T> A);

template<typename T>
SVD<T> svd(Tensor<T> A, size_t k);

/// Eigenvectors/-values
template <typename T>
using SpectralDecomposition = pair<Tensor<T>, Tensord>;
typedef SpectralDecomposition<complex<double>> SpectralDecompositioncd;
typedef SpectralDecomposition<double> SpectralDecompositiond;

template <typename T>
void diagonalize(SpectralDecomposition<T>& x);


#endif //TENSORBLAS2_H
