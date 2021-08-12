//
// Created by Roman Ellerbrock on 8/7/21.
//

#ifndef CONTRACTSOP_H
#define CONTRACTSOP_H
#include "TreeOperators/TensorOperators/TensorOperatorTree.h"
#include "TreeOperators/TensorOperators/TTNOMatrixTree.h"
#include "TreeOperators/TensorOperators/TTNOHoleTree.h"

TensorOperatorTree contractSOP(TensorOperatorTree A, const SOPd& S, const Tree& optree);

void iterate(TensorOperatorTree& A, const SOPd& S, const Tree& optree);

double error(const TensorOperatorTree& A, const SOPd& S, const Tree& optree);

#endif //CONTRACTSOP_H
