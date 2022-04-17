//
// Created by Roman Ellerbrock on 12/4/21.
//

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#include "TensorTree.h"
#include "TensorTreeFactory.h"
#include "Operator/ProductOperator.h"
#include "Operator/SumOfProductsOperator.h"

/**
 * Node targeting routines
 */
template<typename T>
void apply(Tensor<T>& Ket, const TensorTree<T>& pmat,
	const ProductOperator<T>& P, const Node* node);

template<typename T>
void apply(Tensor<T>& Ket, const vector<TensorTree<T>>& S,
	const SOP<T>& H, const Node *node);

/**
 * Edge targeting routines
 */

template <typename T>
Tensor<T> matrixTensor(const Tensor<T>& h, const Tensor<T>& Ket, const Edge& edge);

template <typename T>
Tensor<T> contraction(const Tensor<T>& bra, const Tensor<T>& ket, const Edge& edge);

template<typename T>
void contraction(TensorTree<T>& S, const Tensor<T>& Bra, Tensor<T> Ket,
	const ProductOperator<T>& P, const Edge* edge);

template<typename T>
void contraction(vector<TensorTree<T>>& Svec, const Tensor<T>& Bra,
	Tensor<T> Ket, const SumOfProductsOperator<T>& H, const Edge* edge);

template<typename T>
void apply(Tensor<T>& Ket, const TensorTree<T>& S,
	const ProductOperator<T>& P, const Edge *edge);

template<typename T>
void apply(Tensor<T>& Ket, const vector<TensorTree<T>>& S,
	const SOP<T>& H, const Edge *edge);

/**
 * TensorTree targeting routines
 */


template <typename T>
void contraction(TensorTree<T>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const ProductOperator<T>& P = ProductOperator<T>());

template <typename T>
void contraction(vector<TensorTree<T>>& S, const TensorTree<T>& Bra, TensorTree<T> Ket,
	const SumOfProductsOperator<T>& H = SumOfProductsOperator<T>());

template<typename T>
void apply(TensorTree<T>& Ket, const TensorTree<T>& pmat,
	const ProductOperator<T>& P);

template <typename T>
double residual(const TensorTree<T>& Psi1, const TensorTree<T>& Psi2, const Tree& tree);

#endif //CONTRACTIONS_H
