//
// Created by Roman Ellerbrock on 12/4/21.
//

#ifndef DOTPRODUCT_H
#define DOTPRODUCT_H

#include "TensorTree.h"

template <typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T> Ket, const Edge& edge);

template <typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket, const Edge& edge);

template <typename T>
void dotProduct(TensorTree<T>& S, const TensorTree<T>& Bra, const TensorTree<T>& Ket);

#endif //DOTPRODUCT_H
