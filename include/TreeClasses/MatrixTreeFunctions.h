//
// Created by Roman Ellerbrock on 2/11/20.
//

#ifndef MATRIXTREEFUNCTIONS_H
#define MATRIXTREEFUNCTIONS_H
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/MatrixTree.h"

namespace TreeFunctions {

	template<typename T>
	void dotProductLocal(MatrixTree<T>& S, const Tensor<T>& Phi, Tensor<T> AChi, const Node& node);

	template<typename T>
	void dotProduct(MatrixTree<T>& S, const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree);

	template<typename T>
	MatrixTree<T> dotProduct(const TensorTree<T>& Psi, const TensorTree<T>& Chi, const Tree& tree);

	template<typename T>
//	void ContractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
//		const MatrixTree<T>& S, const Node& node);
	void contractionLocal(MatrixTree<T>& Rho, const Tensor<T>& Bra, Tensor<T> Ket,
		const Node& node, const MatrixTree<T> *S_opt);

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi,
		const TensorTree<T>& Chi, const MatrixTree<T>& S, const Tree& tree);

	template<typename T>
	MatrixTree<T> contraction(const TensorTree<T>& Psi, const TensorTree<T>& Chi,
		const MatrixTree<T>& S, const Tree& tree);

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi, const Tree& tree, bool orthogonal);

	template<typename T>
	MatrixTree<T> contraction(const TensorTree<T>& Psi, const Tree& tree, bool orthogonal);

}

#endif //MATRIXTREEFUNCTIONS_H
