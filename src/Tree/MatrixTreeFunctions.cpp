//
// Created by Roman Ellerbrock on 2/11/20.
//

#include "MatrixTreeFunctions_Implementation.h"

namespace MatrixTreeFunctions {
	typedef complex<double> cd;

	template void DotProductLocal(MatrixTree<cd>& S, const Tensor<cd>& Bra, Tensor<cd> Ket, const Node& node);
	template void DotProduct(MatrixTree<cd>& S, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const TTBasis& tree);
	template MatrixTree<cd> DotProduct(const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const TTBasis& tree);

	template void ContractionLocal(MatrixTree<cd>& Rho, const Tensor<cd>& Bra, Tensor<cd> Ket, const Node& node,
		const MatrixTree<cd>* S);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const TTBasis& tree, const MatrixTree<cd>* S);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const MatrixTree<cd>& S, const TTBasis& tree);
	template void Contraction(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi, const TTBasis& tree, bool orthogonal);
	template MatrixTree<cd> Contraction(const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const MatrixTree<cd>& S, const TTBasis& tree);
	template MatrixTree<cd> Contraction(const TensorTree<cd>& Psi, const TTBasis& tree, bool orthogonal);


	typedef double d;

	template void DotProductLocal<d>(MatrixTree<d>& S, const Tensor<d>& Bra, Tensor<d> Ket, const Node& node);
	template void DotProduct<d>(MatrixTree<d>& S, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const TTBasis& tree);
	template MatrixTree<d> DotProduct(const TensorTree<d>& Bra, const TensorTree<d>& Ket, const TTBasis& tree);

	template void ContractionLocal(MatrixTree<d>& Rho, const Tensor<d>& Bra, Tensor<d> Ket, const Node& node,
		const MatrixTree<d>* S);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Bra, const TensorTree<d>& Ket, const TTBasis& tree,
		const MatrixTree<d>* S);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const MatrixTree<d>& S, const TTBasis& tree);
	template void Contraction(MatrixTree<d>& Rho, const TensorTree<d>& Psi, const TTBasis& tree, bool orthogonal);
	template MatrixTree<d> Contraction(const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const MatrixTree<d>& S, const TTBasis& tree);
	template MatrixTree<d> Contraction(const TensorTree<d>& Psi, const TTBasis& tree, bool orthogonal);
}
