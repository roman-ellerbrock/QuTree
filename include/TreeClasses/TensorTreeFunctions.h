//
// Created by Roman Ellerbrock on 5/24/20.
//

#ifndef TENSORTREEFUNCTIONS_H
#define TENSORTREEFUNCTIONS_H
#include "TreeClasses/SpectralDecompositionTree.h"

namespace TreeFunctions {
	template <typename T>
	void Adjust(TensorTree<T>& Psi, Tree& tree,
		const SpectralDecompositionTree<T>& X, double eps);


	template <typename T>
	void ToplayerSum(Tensor<T>& A, const Tensor<T>& B);

	template <typename T>
	void DirectSum(TensorTree<T>& Psi, Tree& tree, const TensorTree<T>& Chi);

}

#endif //TENSORTREEFUNCTIONS_H
