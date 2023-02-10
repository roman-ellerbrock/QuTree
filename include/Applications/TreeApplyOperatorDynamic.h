//
// Created by Roman Ellerbrock on 2/20/20.
//

#ifndef TREEAPPLYOPERATORDYNAMIC_H
#define TREEAPPLYOPERATORDYNAMIC_H
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeOperators/SumOfProductsOperator.h"
#include "TreeOperators/SOPVector.h"
#include "Applications/TreeApplyOperator.h"

namespace TreeFunctions {
	using namespace TreeFunctions;

	void applyOperator(TensorTreecd& Psi, Tree& tree,
		const SOPVectorcd& sop, const SOPVectorcd& sop_adj, double eps,
		size_t max_spf, size_t n_plus);

	void applyOperator(TensorTreecd& Psi, Tree& tree,
		const SOPcd& sop, const SOPcd& sop_adj,
		double eps, size_t max_spf, size_t n_plus);

	Tensorcd applyLowerDynamic(Tensorcd Phi, const SOPcd& sop,
		vector<SparseMatrixTreecd>& Hmats, const vector<SparseMatrixTreecd>& Holes,
		Node& node, Node& parent, Tensorcd& upPhi, double eps, size_t max_spf, size_t n_plus);

	size_t nOccupied_local(const Vectord& ev, double rho);

}


#endif //TREEAPPLYOPERATORDYNAMIC_H
