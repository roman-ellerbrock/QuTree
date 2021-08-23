//
// Created by Roman Ellerbrock on 8/6/21.
//

#include "TreeOperators/TensorOperators/TTOMatrixTree.h"

Tensord toTensor(const SOPd& S, const Leaf& leaf) {
	const auto& node = (const Node&) leaf.parent();
	const TensorShape& shape = node.shape();
	TensorShape shape_sop(shape);
	shape_sop.setDimension(S.size(), shape.lastIdx());
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

