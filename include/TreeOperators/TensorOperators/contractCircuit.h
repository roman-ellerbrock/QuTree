//
// Created by Roman Ellerbrock on 8/21/21.
//

#ifndef CONTRACTCIRCUIT_H
#define CONTRACTCIRCUIT_H
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeOperators/SOPVector.h"

template<typename T>
[[nodiscard]] TensorTreeOperator<T> contractCircuit(const SOP<T>& S, const TensorTreeOperator<T>& U,
	size_t maxIter, const Tree& optree, ostream *os);

template<typename T>
[[nodiscard]] TensorTreeOperator<T> contractCircuit(const SOPVector<T>& circuit,
	size_t maxIter, const Tree& optree, ostream *os);

#endif //CONTRACTCIRCUIT_H
