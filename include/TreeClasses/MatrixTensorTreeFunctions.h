//
// Created by Roman Ellerbrock on 1/3/21.
//

#ifndef MATRIXTENSORTREEFUNCTIONS_H
#define MATRIXTENSORTREEFUNCTIONS_H
#include "MatrixTensorTree.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"

typedef pair<SparseMatrixTreecd, SparseMatrixTreecd> SparseMatrixTreePaircd;
typedef vector<SparseMatrixTreePaircd> SparseMatrixTreePairscd;

namespace TreeFunctions {

	void Represent(SparseMatrixTreecd& mat, const MatrixTensorTree& Psi,
		const MLOcd& M, const Tree& tree);

	void Contraction(SparseMatrixTreecd& hole, const MatrixTensorTree& Psi,
		const SparseMatrixTreecd& mat, const SparseTree& marker, const Tree& tree);

	void Represent(SparseMatrixTreePaircd& mats, const MatrixTensorTree& Psi,
		const MLOcd& M, const Tree& tree);

	void Represent(SparseMatrixTreePairscd& matset, const MatrixTensorTree& Psi,
		const SOPcd& H, const Tree& tree);

	Tensorcd symApplyDown(const Tensorcd& Phi, const SparseMatrixTreecd& hHole,
		const Node& node);

	Tensorcd symApply(const Tensorcd& Phi,
		const SparseMatrixTreePaircd& mats, const MLOcd& M, const Node& node);

	Tensorcd symApply(Tensorcd Phi,
		const SparseMatrixTreePairscd& hMatSet,
		const SOPcd& H, const Node& node);
}

#endif //MATRIXTENSORTREEFUNCTIONS_H
