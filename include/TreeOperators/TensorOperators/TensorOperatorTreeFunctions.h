//
// Created by Roman Ellerbrock on 5/23/20.
//

#ifndef TENSOROPERATORTREEFUNCTIONS_H
#define TENSOROPERATORTREEFUNCTIONS_H
#include "TreeOperators/TensorOperators/MatrixListTree.h"
#include "TreeOperators/TensorOperators/TensorOperatorTree.h"

namespace TreeFunctions {

	void DirectSum(TensorOperatorTree& H, Tree& tree, const TensorOperatorTree& G);

	MatrixListTree Represent(const TensorTreecd& Psi, const TensorOperatorTree& H,
		const Tree& tree);

	MatrixListTree Contraction(const TensorTreecd& Psi, const TensorOperatorTree& H,
		MatrixListTree& Hrep, const Tree& tree);

}

#endif //TENSOROPERATORTREEFUNCTIONS_H
