//
// Created by Roman Ellerbrock on 11/7/21.
//

#ifndef TENSORBLAS1_H
#define TENSORBLAS1_H
#include "Tensor.h"
#include <blas.hh>

//////////////////////////////////////////////////////////
/// map Lvl 1 BLAS to Tensors
//////////////////////////////////////////////////////////

template<typename T>
T nrm2(const Tensor<T>& A, size_t incr = 1);

template<typename T>
void axpy(const Tensor<T>& A, Tensor<T>& B, T alpha = 1., size_t inc_a = 1, size_t inc_b = 1);

template<typename T>
void operator+=(Tensor<T>& A, const Tensor<T>& B);

template<typename T>
void operator-=(Tensor<T>& A, const Tensor<T>& B);

template<typename T>
T residual(Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Tensor<T> operator+(const Tensor<T>& A, const Tensor<T>& B);

template<typename T>
Tensor<T> operator-(const Tensor<T>& A, const Tensor<T>& B);



template<typename T>
Tensor<T>& operator*=(Tensor<T>& A, T a);

template<typename T>
Tensor<T>& operator/=(Tensor<T>& A, T a);


#endif //TENSORBLAS1_H
