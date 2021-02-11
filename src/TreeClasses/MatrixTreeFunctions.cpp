//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "TreeClasses/MatrixTreeFunctions_Implementation.h"
#include "TreeClasses/TreeTransformation_Implementation.h"

namespace TreeFunctions {
	typedef complex<double> cd;

	template void dotProductLocal(MatrixTree<cd>& S, const Tensor<cd>& Bra, Tensor<cd> Ket, const Node& node);
	template void dotProduct(MatrixTree<cd>& S, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const Tree& tree);
	template MatrixTree<cd> dotProduct(const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void contractionLocal(MatrixTree<cd>& Rho, const Tensor<cd>& Bra, Tensor<cd> Ket, const Node& node,
		const MatrixTree<cd> *S);
	template void contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const Tree& tree, const MatrixTree<cd> *S);
	template void contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const MatrixTree<cd>& S, const Tree& tree);
	template void contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const Tree& tree, bool orthogonal);
	template MatrixTree<cd> contraction(const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const MatrixTree<cd>& S, const Tree& tree);
	template MatrixTree<cd> contraction(const TensorTree<cd>& Psi, const Tree& tree, bool orthogonal);

	typedef double d;

	template void dotProductLocal<d>(MatrixTree<d>& S, const Tensor<d>& Bra, Tensor<d> Ket, const Node& node);
	template void dotProduct<d>(MatrixTree<d>& S, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const Tree& tree);
	template MatrixTree<d> dotProduct(const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void contractionLocal(MatrixTree<d>& Rho, const Tensor<d>& Bra, Tensor<d> Ket, const Node& node,
		const MatrixTree<d> *S);
	template void contraction(MatrixTree<d>& Rho, const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree,
		const MatrixTree<d> *S);
	template void contraction(MatrixTree<d>& Rho, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const MatrixTree<d>& S, const Tree& tree);
	template void contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const Tree& tree, bool orthogonal);
	template MatrixTree<d> contraction(const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const MatrixTree<d>& S, const Tree& tree);
	template MatrixTree<d> contraction(const TensorTree<d>& Psi, const Tree& tree, bool orthogonal);
}
