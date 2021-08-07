//
// Created by Roman Ellerbrock on 8/4/21.
//

#ifndef TTNOMATRIXTREE_H
#define TTNOMATRIXTREE_H
#include "TensorOperatorTree.h"
#include "TreeClasses/MatrixTree.h"
#include "TreeOperators/SumOfProductsOperator.h"

void toTensor(Tensord& A, const MLOd& M, size_t part, const TensorShape& shape, const Leaf& leaf) {
	/**
	 * Rationale:
	 * Build Tensor representation of the LeafOperator that acts on 'leaf' and store it in A(:,part).
	 * Note: currently only works if there is a single operator acting on 'leaf' in M
	 */
	if (!M.isActive(leaf.mode())) {
		auto dim = (size_t) sqrt((double) shape.lastBefore() + 1e-12);
		auto ml = identityMatrixd(dim);
		for (size_t i = 0; i < shape.lastBefore(); ++i) {
			A(i, part) = ml[i];
		}
	}
	for (size_t k = 0; k < M.size(); ++k) {
		if (M.isActive(k, leaf.mode())) {
			const auto& L = M[k];
			const Matrixd ml = toMatrix(*L, leaf);
			for (size_t i = 0; i < shape.lastBefore(); ++i) {
				A(i, part) = ml[i];
			}
		}
	}
}

class TTNOMatrixTree: public MatrixTreed {
public:
	using MatrixTreed::NodeAttribute<Matrix<double>>::attributes_;

	TTNOMatrixTree(const Tree& tree, const SOPd& S) {
		size_t npart = S.size();
		for (const Node& node : tree) {
			const TensorShape& shape = node.shape();
			size_t ntensor = shape.lastDimension();
			attributes_.emplace_back(Matrixd(ntensor, npart));
		}
	}

	void calculate(const TensorOperatorTree& A, const SOPd& S, const Tree& tree) {
		for (const Node& node : tree) {
			Tensord& B();
			if (node.isBottomlayer()) {
				for (size_t l = 0; l < S.size(); ++l) {
					const MLOd& M = S[l];
				}
			} else {
			}
		}
	}

	~TTNOMatrixTree() = default;
};


#endif //TTNOMATRIXTREE_H
