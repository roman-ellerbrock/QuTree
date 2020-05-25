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
		/// do not switch to node.ChildIdx here:
		/// This does not work with TensorOperators
		shape.setDimension(n_occ, shape.lastIdx());
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

	template <typename T>
	void ToplayerSum(Tensor<T>& A, const Tensor<T>& B) {
		const TensorShape& Ashape = A.shape();
		const TensorShape& Bshape = B.shape();
		auto dims = Ashape.dimensions();
		for (size_t i = 0; i < Bshape.order() - 1; ++i) {

		}
	}

	template <typename T>
	void DirectSum(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi) {
		for (Node& node : tree) {
			Psi[node] = Tensor_Extension::DirectSum(Psi[node], Chi[node],
				!node.isBottomlayer(), !node.isToplayer());
			node.shape() = Psi[node].shape();
		}
	}
}


#endif //TENSORTREEFUNCTIONS_IMPLEMENTATION_H
