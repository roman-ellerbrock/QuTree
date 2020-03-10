//
// Created by Roman Ellerbrock on 3/10/20.
//

#ifndef SOPMATRIXTREES_H
#define SOPMATRIXTREES_H
#include "TreeClasses/SparseMatrixTreeFunctions.h"

template <typename T>
class SOPMatrixTrees {
public:
	SOPMatrixTrees(const SOP<T>& H, const TensorTree<T>& Psi, const Tree& tree) {
		for (const auto& M : H) {
			matrices_.push_back(SparseMatrixTree<T>(M, tree));
			contractions_.push_back(SparseMatrixTree<T>(M, tree));
		}
	}

	void Update(const SOP<T>& H, const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		SparseMatrixTreeFunctions::Represent(matrices_, H, Bra, Ket, tree);
		SparseMatrixTreeFunctions::Contraction(contractions_, matrices_, H, Bra, Ket, tree);
	}

private:
	vector<SparseMatrixTree<T>> matrices_;
	vector<SparseMatrixTree<T>> contractions_;
};


#endif //SOPMATRIXTREES_H
