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

	TensorTreecd bottomUpNormalized(const Tree& tree);
	void canonicalRepresentation(const Tree& tree);
	void normalizeCanonical(const Node& node);
	void normalizeCanonical(const Tree& tree);

	void print(const Tree& tree);

	TensorTreecd weighted_;
	TensorTreecd up_;
	TensorTreecd down_;
private:
};

Tensorcd canonicalTensor(Tensorcd w, bool up, bool down);
Matrixcd calculateB(const Tensorcd& weighted, size_t k);

bool isWorking(const SymTensorTree& Psi, const Tree& tree);

template<typename T>
void createWeighted(TensorTree<T>& weighted, TensorTree<T>& down, const TensorTreecd& up, const Tree& tree);

namespace TreeFunctions {
	void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Bra,
		const Tensorcd& Ket, const SparseMatrixTreecd& mat, const Node& hchild);

	void symContraction(SparseMatrixTreecd& hole, const TensorTreecd& Bra,
		const TensorTreecd& Ket, const SparseMatrixTreecd& mat, const Tree& tree);

	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Bra,
		const SymTensorTree& Ket,
		const MLOcd& M, const Tree& tree);

	void symRepresent(SymMatrixTrees& mats, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const SOPcd& S, const Tree& tree);

	Tensorcd symApply(const Tensorcd& Ket,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node);

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree);
}

#endif //SYMTENSORTREE_H
