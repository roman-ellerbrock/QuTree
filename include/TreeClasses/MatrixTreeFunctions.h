//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREEFUNCTIONS_H
#define MATRIXTREEFUNCTIONS_H
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/MatrixTree.h"

namespace MatrixTreeFunctions {

	template<typename T>
	void DotProductLocal(MatrixTree<T>& S, const Tensor<T>& Phi, Tensor<T> AChi, const Node& node);

	template<typename T>
	void DotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree);

	template<typename T>
	MatrixTree<T> DotProduct(const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree);

	template<typename T>
	void ContractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
		const MatrixTree<T>& S, const Node& node);

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const MatrixTree<T>& S, const Tree& tree);

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree);

	template<typename T>
	void Contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const Tree& tree, bool orthogonal);

	template<typename T>
	MatrixTree<T> Contraction(const TensorTree<T>& Psi, const Tree& tree, bool orthogonal);

}

#endif //MATRIXTREEFUNCTIONS_H
