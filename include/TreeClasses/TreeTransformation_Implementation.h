//
// Created by Roman Ellerbrock on 4/9/20.
//

#ifndef TREETRANSFORMATION_IMPLEMENTATION_H
#define TREETRANSFORMATION_IMPLEMENTATION_H

#include "TreeClasses/TreeTransformation.h"
#include "TreeClasses/SpectralDecompositionTree.h"


namespace TreeFunctions {


	template <typename T>
	TensorTree<T> DotProductNormalization(TensorTree<T> Psi, const Tree& tree) {

	}

	template <typename T>
	TensorTree<T> ContractionNormalization(TensorTree<T> Psi, const Tree& tree, bool orthogonal) {
		/// Calculate contraction matrixtree
		MatrixTree<T> rho = Contraction(Psi, tree, orthogonal);
		auto rho_x = SpectralDecompositionTree<T>(rho);

		/// sqrt(rho)
		auto sqrt_rho_x = sqrt(rho_x);
		auto sqrt_rho = to_matrixtree(sqrt_rho_x);

		/// Inverse of sqrt(rho)
		auto inv_sqrt_rho_x = inverse(sqrt_rho_x);
		auto inv_sqrt_rho = to_matrixtree(inv_sqrt_rho_x);
	}

	template <typename T>
	void Transform(TensorTree<T>& Psi, const MatrixTree<T>& M, const MatrixTree<T>& M_inv,
		const Tree& tree) {
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
//				multStateAB(M[node], Psi[node], node.parent)
				Psi[node] = MatrixTensor(M[node], Psi[node], node.parentIdx());

				const Node& parent = node.parent();
			}
		}
	}

}

#endif //TREETRANSFORMATION_IMPLEMENTATION_H
