//
// Created by Roman Ellerbrock on 1/27/21.
//

#ifndef SYMTENSORTREE_H
#define SYMTENSORTREE_H
#include "MatrixTensorTreeFunctions.h"

typedef SparseMatrixTreePaircd SymMatrixTree;
typedef SparseMatrixTreePairscd SymMatrixTrees;

class SymTensorTree {
public:
	SymTensorTree(const Tree& tree) : weighted_(tree), up_(tree), down_(tree) {}
	SymTensorTree(const TensorTreecd& Psi, const Tree& tree) {
		initializeFromTT(Psi, tree);
	}
	~SymTensorTree() {}

	void initializeFromTT(const TensorTreecd& Psi, const Tree& tree);

	void rebuild(const Tree& tree);

	TensorTreecd weighted_;
	TensorTreecd up_;
	TensorTreecd down_;
private:
};

Tensorcd solveSLE(const Matrixcd& B, const Tensorcd& A, size_t idx);

Tensorcd normalizedTensor(const Tensorcd& weighted, size_t k);

bool isWorking(const SymTensorTree& Psi, const Tree& tree);

template<typename T>
void QROrthogonalDown(TensorTree<T>& weighted, TensorTree<T>& down, const TensorTreecd& up, const Tree& tree);

namespace TreeFunctions {
	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Psi,
		const MLOcd& M, const Tree& tree);

	void symRepresent(SymMatrixTrees& mat, const SymTensorTree& Psi,
		const SOPcd& M, const Tree& tree);

/*	Tensorcd symApplyDown(const Tensorcd& Phi, const SparseMatrixTreecd& hHole,
		const Node& node);

	Tensorcd symApply(const Tensorcd& Phi,
		const SparseMatrixTreePaircd& mats, const MLOcd& M, const Node& node);

	Tensorcd symApply(Tensorcd Phi,
	const SparseMatrixTreePairscd& hMatSet,
	const SOPcd& H, const Node& node);
*/

	Tensorcd symApply(const Tensorcd& Phi,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node);

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree);
}

#endif //SYMTENSORTREE_H
