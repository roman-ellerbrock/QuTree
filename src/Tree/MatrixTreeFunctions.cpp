//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "MatrixTreeFunctions_Implementation.h"

namespace MatrixTreeFunctions {
	typedef complex<double> cd;

	template void DotProductLocal(MatrixTree<cd>& S, const Tensor<cd>& Phi, Tensor<cd> AChi, const Node& node);
	template void DotProduct(MatrixTree<cd>& S, const TensorTree<cd>& Psi, const TensorTree<cd>& Chi,
		const TTBasis& basis);
	template MatrixTree<cd> DotProduct(const TensorTree<cd>& Psi, const TensorTree<cd>& Chi, const TTBasis& basis);

	template void ContractionLocal(MatrixTree<cd>& Rho, const Tensor<cd>& Bra, Tensor<cd> Ket, const Node& node,
		const MatrixTree<cd>* S);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const TensorTree<cd>& Chi,
		const TTBasis& basis, const MatrixTree<cd>* S);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const TensorTree<cd>& Chi,
		const MatrixTree<cd>& S, const TTBasis& basis);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const TTBasis& basis, bool orthogonal);
	template MatrixTree<cd> Contraction(const TensorTree<cd>& Psi, const TensorTree<cd>& Chi,
		const MatrixTree<cd>& S, const TTBasis& basis);
	template MatrixTree<cd> Contraction(const TensorTree<cd>& Psi, const TTBasis& basis, bool orthogonal);

	typedef double d;

	template void DotProductLocal<d>(MatrixTree<d>& S, const Tensor<d>& Phi, Tensor<d> AChi, const Node& node);
	template void DotProduct<d>(MatrixTree<d>& S, const TensorTree<d>& Psi, const TensorTree<d>& Chi,
		const TTBasis& basis);
	template MatrixTree<d> DotProduct(const TensorTree<d>& Psi, const TensorTree<d>& Chi, const TTBasis& basis);

	template void ContractionLocal(MatrixTree<d>& Rho, const Tensor<d>& Bra, Tensor<d> Ket, const Node& node,
		const MatrixTree<d>* S);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const TensorTree<d>& Chi, const TTBasis& basis,
		const MatrixTree<d>* S);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const TensorTree<d>& Chi,
		const MatrixTree<d>& S, const TTBasis& basis);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const TTBasis& basis, bool orthogonal);
	template MatrixTree<d> Contraction(const TensorTree<d>& Psi, const TensorTree<d>& Chi,
		const MatrixTree<d>& S, const TTBasis& basis);
	template MatrixTree<d> Contraction(const TensorTree<d>& Psi, const TTBasis& basis, bool orthogonal);
}
