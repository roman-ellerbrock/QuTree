//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef MATRIXTREETRANSFORMATION_H
#define MATRIXTREETRANSFORMATION_H
#include "TreeClasses/MatrixTreeFunctions.h"

namespace TreeFunctions {

	template <typename T>
	TensorTree<T> DensityWeighting(TensorTree<T> Psi, const Tree& tree);

	template <typename T>
	TensorTree<T> ContractionNormalization(TensorTree<T> Psi, const Tree& tree);

}

#endif //MATRIXTREETRANSFORMATION_H
