//
// Created by Roman Ellerbrock on 4/2/22.
//

#ifndef TENSORTREEFACTORY_H
#define TENSORTREEFACTORY_H
#include "TensorTree.h"
#include "Operator/SumOfProductsOperator.h"

template<typename T>
TensorTree<T> matrixTree(const Tree& tree, const SubTreeParameters& par = {});

template <typename T>
TensorTree<T> matrixTree(const Tree& tree, const ProductOperator<T>& P);

template <typename T>
vector<TensorTree<T>> matrixTree(const Tree& tree, const SumOfProductsOperator<T>& H);

template <typename T>
TensorTree<T> matrixTree(const TensorTree<T>& Psi, const ProductOperator<T>& P = {});

/**
 * The following code allows to use matrixTreecd(...) for overloaded template functions
 * using perfect forwarding.
 * See this: https://stackoverflow.com/a/9864472/1407466
 *
 * This (simpler code) does not work:
 * constexpr auto matrixTreecd = matrixTree<complex<double>>;
 * constexpr auto matrixTreed = matrixTree<double>;
 */

/// matrixTree
template <typename... Args>
auto matrixTreecd(Args&&... args)
-> decltype(matrixTree<complex<double>>(std::forward<Args>(args)...)) {
	return matrixTree<complex<double>>(std::forward<Args>(args)...);
}

template <typename... Args>
auto matrixTreed(Args&&... args)
-> decltype(matrixTree<double>(std::forward<Args>(args)...)) {
	return matrixTree<double>(std::forward<Args>(args)...);
}

#endif //TENSORTREEFACTORY_H
