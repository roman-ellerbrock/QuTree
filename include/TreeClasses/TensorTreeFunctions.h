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
}

#endif //TENSORTREEFUNCTIONS_H
