//
// Created by Roman Ellerbrock on 8/6/21.
//

#include "TreeOperators/TensorOperators/TTNOMatrixTree.h"

void toTensor(Tensord& A, const MLOd& M, size_t part, const Leaf& leaf) {
	/**
	 * Rationale:
	 * Build Tensor representation of the LeafOperator that acts on 'leaf' and store it in A(:,part).
	 * Note: currently only works if there is a single operator acting on 'leaf' in M
	 */
	const TensorShape& shape = A.shape();
	auto dim = (size_t) sqrt((double) shape.lastBefore() + 1e-12);
	auto sigma = identityMatrixd(dim);
	for (size_t k = 0; k < M.size(); ++k) {
		if (M.isActive(k, leaf.mode())) {
			const auto& L = M[k];
			const Matrixd mk = toMatrix(*L, leaf);
			sigma = mk * sigma;
		}
	}
	for (size_t i = 0; i < shape.lastBefore(); ++i) {
		A(i, part) = sigma[i];
	}
}

Tensord toTensor(const SOPd& S, const Leaf& leaf) {
	const auto& node = (const Node&) leaf.parent();
	const TensorShape& shape = node.shape();
	TensorShape shape_sop(shape);
	size_t idx = shape_sop.lastIdx();
	shape_sop.setDimension(S.size(), idx);
	Tensord C(shape_sop);
	for (size_t l = 0; l < S.size(); ++l) {
		const MLOd& M = S[l];
		toTensor(C, M, l, leaf);
	}
	return C;
}

double prodMk(const vector<size_t>& idx, const vector<Matrixd>& Mk, size_t l, int skip) {
	double factor = 1.;
	for (size_t k = 0; k < Mk.size(); ++k) {
		if (k != skip) {
			factor *= Mk[k](idx[k], l);
		}
	}
	return factor;
}

