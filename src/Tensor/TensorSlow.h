//
// Created by Roman Ellerbrock on 11/17/21.
//

#ifndef TENSORSLOW_H
#define TENSORSLOW_H
#include "Tensor/Tensor.h"
#include <blas.hh>
#include "Tensor/SVD.h"

template<typename T>
void gemmRef(Tensor<T>& c, Tensor<T> a, Tensor<T> b,
	T alpha = 1., T beta = 0.,
	blas::Op op_a = blas::Op::NoTrans, blas::Op op_b = blas::Op::NoTrans);

template<typename T>
Tensor<T> gemmRef(Tensor<T> a, Tensor<T> b, T alpha = 1.,
	blas::Op op_a = blas::Op::NoTrans, blas::Op op_b = blas::Op::NoTrans);

template<typename T, typename U>
void matrixTensorRef(Tensor<T>& C, const Tensor<U>& h, const Tensor<T>& B,
	size_t k, bool zero = true);

template<typename T>
void contractionRef(Tensor<T>& h, const Tensor<T>& bra, const Tensor<T>& ket,
	size_t k, T alpha = 1., T beta = 0.);

template <typename T>
Tensor<T> toTensorRef(const SVD<T>& x);

template<typename T>
Tensor<T> toTensorRef(const SpectralDecomposition<T>& X);


#endif //TENSORSLOW_H
