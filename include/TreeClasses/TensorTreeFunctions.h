//
// Created by Roman Ellerbrock on 5/24/20.
//

#ifndef TENSORTREEFUNCTIONS_H
#define TENSORTREEFUNCTIONS_H
#include "TreeClasses/SpectralDecompositionTree.h"

namespace TreeFunctions {
	/// Compress the TensorTree for a given accuracy.
	template <typename T>
	void adjust(TensorTree<T>& Psi, Tree& tree,
		const SpectralDecompositionTree<T>& X, double eps);

	template <typename T>
	void adjust(TensorTree<T>& Psi, const Tree& newtree);

	/// Perform a (generalized) sum of tensor trees. By default, a regular sum will be performed.
	template <typename T>
	void sum(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi,
		bool sharedLeafs = true, bool sumToplayer = true);

	/// Perform a product of tensor trees.
	template <typename T>
	void product(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi);

}

#endif //TENSORTREEFUNCTIONS_H
