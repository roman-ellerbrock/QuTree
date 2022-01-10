//
// Created by Roman Ellerbrock on 8/24/21.
//
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeOperators/TensorOperators/TTOHoleTree.h"
#include "TreeOperators/TensorOperators/TTOMatrixTree.h"
#include "TreeOperators/TensorOperators/TTOcontraction_Implementation.h"
#include "TreeOperators/TensorOperators/TTOrepresentation.h"
#include "TreeOperators/TensorOperators/contractSOP.h"

typedef double d;
typedef complex<double> cd;

/// Check other .cpp's for more instantiations (e.g. TensorTreeOperator.cpp)

template class TTOHoleTree<d>;
template class TTOHoleTree<cd>;

template class TTOMatrixTree<d>;
template class TTOMatrixTree<cd>;

template class TTOrepresentation<d>;
template class TTOrepresentation<cd>;

template class TTOcontraction<d>;
template class TTOcontraction<cd>;

template Tensor<d> apply(const Tensor<d>& Phi, const TensorTreeOperator<d>& H,
	const TTOrepresentation<d>& rep, const TTOcontraction<d>& con, const Node& node);
template Tensor<cd> apply(const Tensor<cd>& Phi, const TensorTreeOperator<cd>& H,
	const TTOrepresentation<cd>& rep, const TTOcontraction<cd>& con, const Node& node);

/*
template TensorTreeOperator<d> contractSOP(TensorTreeOperator<T> A, const SOP<T>& S,
	size_t maxIter, const Tree& optree, ostream *os);

template void iterate(TensorTreeOperator<d>& A, const SOP<T>& S, const Tree& optree);

template double error(const TensorTreeOperator<d>& A, const SOP<T>& S, const Tree& optree);


template TensorTreeOperator<T> contractSOP(TensorTreeOperator<T> A, const SOP<T>& S,
	size_t maxIter, const Tree& optree, ostream *os);

template void iterate(TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& optree);

template double error(const TensorTreeOperator<T>& A, const SOP<T>& S, const Tree& optree);
*/
