//
// Created by Roman Ellerbrock on 11/18/21.
//

#ifndef TENSORLAPACK_H
#define TENSORLAPACK_H
#include "Tensor/TensorBLAS2.h"
#include "Tensor/SVD.h"

/// ===================== Orthonormalization ====================
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

/// ===================== SVD ====================
template<typename T>
Tensor<T> toTensor(const SVD<T>& x, size_t k);

template<typename T>
Tensor<T> toTensor(const SVD<T>& x);

template<typename T>
void svd(SVD<T>& x, Tensor<T> A);

template<typename T>
SVD<T> svd(Tensor<T> A);

template<typename T>
SVD<T> svd(const Tensor<T>& A, size_t k);

template<typename T>
void regularize(SVD<T>& x, size_t k, double eps, mt19937& gen = rng::gen);

template<typename T>
void normalize(SVD<T>& x, size_t k, double eps, mt19937& gen = rng::gen);

template<typename T>
Tensor<T> normalize(const Tensor<T>& A, size_t k, double eps, mt19937& gen = rng::gen);

/// ===================== Diagonalization ====================
template<typename T>
Tensor<T> toTensor(const SpectralDecomposition<T>& X);

template<typename T>
void phaseConvention(Tensor<T>& mat);

template <typename T>
void heev(SpectralDecomposition<T>& x);

template <typename T>
SpectralDecomposition<T> heev(const Tensor<T>& A);


#endif //TENSORLAPACK_H
