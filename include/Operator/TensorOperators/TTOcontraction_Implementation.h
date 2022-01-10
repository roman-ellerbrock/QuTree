//
// Created by Roman Ellerbrock on 9/5/21.
//
#include "TreeOperators/TensorOperators/TTOcontraction.h"

template <typename T>
Tensor<T> apply(const Tensor<T>& Phi, const TensorTreeOperator<T>& H,
	const TTOrepresentation<T>& rep, const TTOcontraction<T>& con,
	const Node& node) {

	Tensor<T> HPhi(Phi.shape());
	const TensorShape& opshape = H[node].shape();
	double eps = 1e-12;
	if (node.isBottomlayer()) {
		const Leaf& leaf = node.getLeaf();
		for (size_t l = 0; l < opshape.lastDimension(); ++l) {
			vector<size_t> ls = {0, l};
			HPhi += con.applyMatrices(Phi, ls, H, rep, node, -1);
		}
	} else {
		auto ls = indexMapping(0, opshape);
		for (size_t L = 0; L < opshape.totalDimension(); ++L) {
			indexMapping(ls, L, opshape);
			if (abs(H[node](ls)) < eps) {continue;}
			HPhi += con.applyMatrices(Phi, ls, H, rep, node, -1);
		}
	}
	return HPhi;
}

