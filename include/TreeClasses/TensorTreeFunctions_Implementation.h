//
// Created by Roman Ellerbrock on 5/24/20.
//

#ifndef TENSORTREEFUNCTIONS_IMPLEMENTATION_H
#define TENSORTREEFUNCTIONS_IMPLEMENTATION_H
#include "TensorTreeFunctions.h"

namespace TreeFunctions {
	size_t nOccupied(const Vectord& p, double eps) {
		size_t n = 0;
		for (size_t i = 0; i < p.Dim(); ++i) {
			if (p(i) > eps) { n++; }
		}
		if (n == 0) { n = 1; }
		return n;
	}

	template <typename T>
	void AdjustNode(Tensor<T>& Phi, Tensor<T>& A, Node& node, Node& parent,
		const SpectralDecomposition<T>& x, double eps) {
		const Vectord& p = x.second;
		size_t n_occ = nOccupied(p, eps);

		/// Adjust node content
		TensorShape& shape = node.shape();
		size_t node_idx = node.parentIdx();
		shape.setDimension(n_occ, node_idx);
		Phi = Phi.AdjustDimensions(shape);

		/// Adjust parent content
		TensorShape& pshape = parent.shape();
		size_t k = node.childIdx();
		Phi = Phi.AdjustDimensions(shape);
		pshape.setDimension(n_occ, k);
		A = A.AdjustDimensions(pshape);

	}

	template <typename T>
	void Adjust(TensorTree<T>& Psi, Tree& tree,
		const SpectralDecompositionTree<T>& X, double eps) {
		for (Node& node : tree) {
			if (!node.isToplayer()) {
				Node& parent = node.parent();
				AdjustNode(Psi[node], Psi[parent], node, parent, X[node], eps);
			}
		}
	}
}


#endif //TENSORTREEFUNCTIONS_IMPLEMENTATION_H
