//
// Created by Roman Ellerbrock on 8/6/21.
//

#include "TreeOperators/TensorOperators/TTOMatrixTree.h"


template <typename T>
Tensor<T> toTensor(const SOP<T>& S, const Leaf& leaf) {
	const auto& node = (const Node&) leaf.parent();
	const TensorShape& shape = node.shape();
	TensorShape shape_sop(shape);
	shape_sop.setDimension(S.size(), shape.lastIdx());
	Tensor<T> C(shape_sop);
	for (size_t l = 0; l < S.size(); ++l) {
		const MLO<T>& M = S[l];
		toTensor(C, M, l, leaf);
	}
	return C;
}

template <typename T>
T prodMk(const vector<size_t>& idx, const vector<Matrix<T>>& Mk, size_t l, int skip) {
	T factor = 1.;
	for (size_t k = 0; k < Mk.size(); ++k) {
		if (k != skip) {
			factor *= Mk[k](idx[k], l);
		}
	}
	return factor;
}

typedef double d;
typedef complex<double> cd;

template Tensor<cd> toTensor(const SOP<cd>& S, const Leaf& leaf);
template Tensor<d> toTensor(const SOP<d>& S, const Leaf& leaf);

template d prodMk(const vector<size_t>& idx, const vector<Matrix<d>>& Mk, size_t l, int skip = -1);
template cd prodMk(const vector<size_t>& idx, const vector<Matrix<cd>>& Mk, size_t l, int skip = -1);

