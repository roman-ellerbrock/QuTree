//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREEMANAGER_H
#define MATRIXTREEMANAGER_H
#include "Tree/MatrixTree.h"
#include "Tree/TensorTree.h"

namespace MatrixTreeManager {

	template<typename T>
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Phi, Tensor<T> AChi, const Node& node);

	template<typename T>
	Matrix<T> DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const TTBasis& basis);

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const MatrixTree<T>& S, const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const TTBasis& basis);

	template<typename T>
	void ContractionLocal(MatrixTree<T>& Rho, const MatrixTree<T>& S, const Tensor<T>& Bra,
		Tensor<T> Ket, const Node& node);
}

#endif //MATRIXTREEMANAGER_H
