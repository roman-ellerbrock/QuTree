//
// Created by Roman Ellerbrock on 1/4/21.
//

#ifndef SYMXMATRIXTREES_H
#define SYMXMATRIXTREES_H
#include "MatrixTensorTree.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "MatrixTensorTreeFunctions.h"
#include "XMatrixTrees.h"
#include "TreeClasses/SymTensorTree.h"
#include "TreeClasses/SymMatrixTreeFunctions.h"

SOPcd symXsop(const Tree& tree);

class SymXMatrixTrees {
public:

	SymXMatrixTrees(const Tree& tree)
	: xops_(symXsop(tree)) {
		xmat_.clear();
		for (const auto& x : xops_) {
			SparseMatrixTreecd mat(x, tree);
			SparseMatrixTreecd hole(x, tree, false, true);
			xmat_.emplace_back(SymMatrixTree(mat, hole));
		}
	}

	~SymXMatrixTrees() = default;

	void update(const SymTensorTree& Psi, const Tree& tree) {
		TreeFunctions::symRepresent(xmat_, Psi, Psi, xops_, tree);
	}

	SOPcd xops_;

	SymMatrixTrees xmat_;
};


#endif //SYMXMATRIXTREES_H
