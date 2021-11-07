//
// Created by Roman Ellerbrock on 1/27/21.
//

#ifndef SYMTENSORTREE_H
#define SYMTENSORTREE_H
#include "MatrixTensorTreeFunctions.h"


class SymMatrixTree : public SparseMatrixTreePaircd {
public:
	~SymMatrixTree() = default;

	SymMatrixTree(const MLOcd& M, const Tree& tree)
	: SparseMatrixTreePaircd({M, tree}, {M, tree}) {
		initialize(M, tree);
	}

	void initialize(const MLOcd& M, const Tree& tree) {
		first = SparseMatrixTreecd(M, tree);
		second = SparseMatrixTreecd(M, tree);
	}

	void print(const Tree& tree) const {
		cout << "Up:\n";
		first.print();
		cout << "Down: " << endl;
		second.print();
	}
};

class SymMatrixTrees : public vector<SymMatrixTree> {
public:
	SymMatrixTrees() = default;
	~SymMatrixTrees() = default;
	SymMatrixTrees(const SOPcd& S, const Tree& tree) {
		initialize(S, tree);
	}

	void initialize(const SOPcd& S, const Tree& tree) {
		clear();
		for (const MLOcd& M : S) {
			push_back(SymMatrixTree(M, tree));
		}
	}

	void print(const Tree& tree) const {
		for (size_t i = 0; i < size(); ++i) {
			cout << "part = " << i << endl;
			const SymMatrixTree& mat = operator[](i);
			mat.print(tree);
		}
	}

};

class SymTensorTree {
public:
	explicit SymTensorTree() : eps_(1e-7) {}
	~SymTensorTree() = default;

	void initialize(const Tree& tree);
	SymTensorTree(mt19937& gen, const Tree& tree, bool delta_lowest = true);

	void normalizeUp(const Tree& tree);
	void normalizeDown(const Tree& tree);
	void normalize(const Tree& tree);

	TensorTreecd toConventional(const Tree& tree) const;

	void print(const Tree& tree) const {
		cout << "weighted:\n";
		weighted_.print(tree);
		cout << "up:\n";
		up_.print(tree);
		cout << "down:\n";
		down_.print(tree);

	}

	TensorTreecd weighted_;
	TensorTreecd up_;
	TensorTreecd down_;
	double eps_;
};

Tensorcd canonicalTensor(Tensorcd w, bool up, bool down);
Matrixcd calculateB(const Tensorcd& weighted, size_t k);

bool isWorking(const SymTensorTree& Psi, const Tree& tree);

void symDotProduct(const SymTensorTree& Bra, const SymTensorTree& Ket, const Tree& tree);

template<typename T>
void createWeighted(TensorTree<T>& weighted, TensorTree<T>& down, const TensorTreecd& up, const Tree& tree);

namespace TreeFunctions {

	void contractionUp(MatrixTreecd& S, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const Tree& tree);

	void contractionDown(MatrixTreecd& Sdown, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const MatrixTreecd& S, const Tree& tree);

	vector<double> dotProduct(const SymTensorTree& Bra, SymTensorTree Ket,
		const Tree& tree);

		void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Bra,
		const Tensorcd& Ket, const SparseMatrixTreecd& mat, const Node& hchild);

	void symContraction(SparseMatrixTreecd& hole, const TensorTreecd& Bra,
		const TensorTreecd& Ket, const SparseMatrixTreecd& mat, const Tree& tree);

	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Bra,
		const SymTensorTree& Ket,
		const MLOcd& M, const Tree& tree);

	void symRepresent(SymMatrixTrees& mats, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const SOPcd& S, const Tree& tree);

	void applySCF(SymTensorTree& HPsi, SymMatrixTrees& mats, const SymTensorTree& Psi,
		const SOPcd& H, const Tree& tree, double eps, size_t max_iter, ostream* os = nullptr);

	Tensorcd symApply(const Tensorcd& Ket,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node);

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree);

	void symApply(Tensorcd& HPhi, const Tensorcd& Phi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Node& node);

	MatrixTreecd symDotProduct(const SymTensorTree& Bra, const SymTensorTree& Ket,
		const Tree& tree);

	void generate(SymTensorTree& psi, mt19937& gen, const SparseTree& stree, const Tree& tree);

}

#endif //SYMTENSORTREE_H
