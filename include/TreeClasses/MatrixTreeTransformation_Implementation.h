//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef MATRIXTREETRANSFORMATION_IMPLEMENTATION_H
#define MATRIXTREETRANSFORMATION_IMPLEMENTATION_H

#include "TreeClasses/MatrixTreeTransformation.h"
#include "TreeClasses/SpectralDecompositionTree.h"


namespace TreeFunctions {

	template <typename T>
	TensorTree<T> ContractionNormalization(TensorTree<T> Psi, const Tree& tree) {
		auto S = DotProduct(Psi, tree);
		auto rho = Contraction(Psi, tree, S);
		auto rho_x = SpectralDecompositionTree<T>(rho);
		auto sqrt_rho = sqrt(rho_x);
		auto inv_sqrt_rho = inverse(sqrt_rho);
	}

}

#endif //MATRIXTREETRANSFORMATION_IMPLEMENTATION_H
