//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef TREETRANSFORMATION_H
#define TREETRANSFORMATION_H
#include "TreeClasses/MatrixTreeFunctions.h"

namespace TreeFunctions {

	template <typename T>
	TensorTree<T> DensityWeighting(TensorTree<T> Psi, const Tree& tree);

	template <typename T>
	TensorTree<T> ContractionNormalization(TensorTree<T> Psi, const Tree& tree, bool orthogonal = true);

	template <typename T>
	void Transform(TensorTree<T>& Chi, const TensorTree<T>& Psi, const MatrixTree<T>& M, const MatrixTree<T>& M_inv,
		const Tree& tree);
}

#endif //TREETRANSFORMATION_H
