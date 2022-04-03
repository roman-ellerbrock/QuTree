//
// Created by Roman Ellerbrock on 8/7/21.
//

#ifndef CONTRACTSOP_H
#define CONTRACTSOP_H
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeOperators/TensorOperators/TTOMatrixTree.h"
#include "TreeOperators/TensorOperators/TTOHoleTree.h"

template <typename T>
void contractSOP(TensorTreeOperator<T>& A, const SOP<T>& S,
	size_t maxIter, const Tree& optree, ostream *os);

template <typename T>
void iterate(TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& optree);

template <typename T>
double error(const TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& optree);

#endif //CONTRACTSOP_H
