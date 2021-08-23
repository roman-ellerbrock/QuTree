//
// Created by Roman Ellerbrock on 8/7/21.
//

#ifndef CONTRACTSOP_H
#define CONTRACTSOP_H
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeOperators/TensorOperators/TTOMatrixTree.h"
#include "TreeOperators/TensorOperators/TTOHoleTree.h"

TensorTreeOperator contractSOP(TensorTreeOperator A, const SOPd& S, size_t maxIter, const Tree& optree, ostream *os);

void iterate(TensorTreeOperator& A, const SOPd& S, const Tree& optree);

double error(const TensorTreeOperator& A, const SOPd& S, const Tree& optree);

#endif //CONTRACTSOP_H
