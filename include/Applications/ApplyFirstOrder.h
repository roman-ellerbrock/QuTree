//
// Created by Roman Ellerbrock on 12/16/20.
//

#ifndef APPLYFIRSTORDER_H
#define APPLYFIRSTORDER_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SOPVector.h"

typedef TensorTreecd Wavefunction;

namespace TreeFunctions {
	void applyOperatorSCF(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPVectorcd& sop, const Tree& tree);

	double applyOperatorIteration(TensorTreecd& Psi, MatrixTreecd& rho,
		const SOPcd& sop, const Tree& tree);

	void applyOperatorIteration(TensorTreecd& HPsi, const TensorTreecd& Psi, MatrixTreecd& rho,
		SparseMatrixTreescd& Hmats, SparseMatrixTreescd& HHoles,
		const SOPcd& sop, const SparseTree& stree, const Tree& tree);

	Tensorcd applyOperatorLocal(const Tensorcd& Phi, const MatrixTreecd& rho,
		const SparseMatrixTreescd& mats, const SparseMatrixTreescd& holes,
		const SOPcd& sop, const Node& node);

	double fidelity(const TensorTreecd& Psi, const TensorTreecd& Psilast,
		const SparseMatrixTreescd& umat, const MatrixTreecd& rho,
		const SOPcd& u, const SparseTree& stree, const Tree& tree);
}

#endif //APPLYFIRSTORDER_H
