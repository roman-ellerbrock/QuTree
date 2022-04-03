//
// Created by Roman Ellerbrock on 12/4/21.
//

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "TensorTree.h"
#include "TensorTreeFactory.h"
#include "Operator/ProductOperator.h"
#include "Operator/SumOfProductsOperator.h"

template <typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& Ket, const Edge& edge);

template <typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket, const Edge& edge);

template <typename T>
void contraction(TensorTree<T>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const ProductOperator<T>& P = ProductOperator<T>());

template <typename T>
void contraction(vector<TensorTree<T>>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const SumOfProductsOperator<T>& H = SumOfProductsOperator<T>());

template <typename T>
TensorTree<T> product(const TensorTree<T>& S, TensorTree<T> Ket);

template <typename T>
TensorTree<T> fullContraction(const TensorTree<T>& Bra, const TensorTree<T>& S, TensorTree<T> Ket);

#endif //CONTRACTIONS_H
