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

	template<typename T>
	void AdjustNode(Tensor<T>& Phi, Tensor<T>& A, Node& node, Node& parent,
		const SpectralDecomposition<T>& x, double eps) {
		const Vectord& p = x.second;
		size_t n_occ = nOccupied(p, eps);

		/// Adjust node content
		TensorShape& shape = node.shape();
		/// do not switch to node.ChildIdx here:
		/// This does not work with TensorOperators
		shape.setDimension(n_occ, shape.lastIdx());
		Phi = Phi.adjustDimensions(shape);

		/// Adjust parent content
		TensorShape& pshape = parent.shape();
		size_t k = node.childIdx();
		Phi = Phi.adjustDimensions(shape);
		pshape.setDimension(n_occ, k);
		A = A.adjustDimensions(pshape);
	}

	template<typename T>
	void Adjust(TensorTree<T>& Psi, Tree& tree,
		const SpectralDecompositionTree<T>& X, double eps) {
		for (Node& node : tree) {
			if (!node.isToplayer()) {
				Node& parent = node.parent();
				AdjustNode(Psi[node], Psi[parent], node, parent, X[node], eps);
			}
		}
	}

	template<typename T>
	void Adjust(TensorTree<T>& Psi, const Tree& newtree) {
		for (const Node& node : newtree) {
			Tensor<T>& Phi = Psi[node];
			Phi = Phi.adjustDimensions(node.shape());
		}
	}

	template<typename T>
	void Sum(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi,
		bool sharedLeafs, bool sumToplayer) {
		/**
		 * \brief Perform sum on TensorTrees.
		 *
		 * This is the generalized sum for Tensor Trees. It performs a direct sum of the basis
		 * and can perform either a regular sum (Psi + Chi) at the toplayer or a direct sum
		 * at toplayer (Psi [+] Chi) depending on the boolean sumToplayer. The Boolean
		 * sharedLeafs tells whether the LeafInterface is the same or not.
		 *
		 * @param Psi left element of sum
		 * @param Chi left element of sum
		 * @param tree Treestructure of Psi & Chi, will get changed during sum.
		 * @param sharedLeafs Tells whether leaf basis is same for Psi & Chi.
		 * @param sumToplayer If true, the toplayer will be a regular sum, otherwise a direct sum.
		 */

		for (Node& node : tree) {
			bool before = !(node.isBottomlayer() && sharedLeafs);
			bool last = !(node.isToplayer() && sumToplayer);
			Psi[node] = Tensor_Extension::directSum(Psi[node], Chi[node],
				before, last);
			node.shape() = Psi[node].shape();
		}
	}

	template <typename T>
	void Product(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi) {
		for (Node& node : tree) {
			Psi[node] = Tensor_Extension::directProduct(Psi[node], Chi[node]);
			node.shape() = Psi[node].shape();
		}
	}
}


#endif //TENSORTREEFUNCTIONS_IMPLEMENTATION_H
