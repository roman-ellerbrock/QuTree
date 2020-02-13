//
// Created by Roman Ellerbrock on 2/12/20.
//
#include "Tree/SparseMatrixTreeFunctions_Implementation.h"

namespace SparseMatrixTreeFunctions {
	typedef complex<double> cd;

	template void Represent(SparseMatrixTree<cd>& hmat,
		const MLO<cd>& M, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const TTBasis& basis);

	template void Represent(SparseMatrixTree<cd>& hmat, const MLO<cd>& M,
		const TensorTree<cd>& Psi, const TTBasis& basis);

	template SparseMatrixTree<cd> Represent(const MLO<cd>& M, const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const TTBasis& basis);

	template SparseMatrixTree<cd> Represent(const MLO<cd>& M, const TensorTree<cd>& Psi,
		const TTBasis& basis);

	typedef double d;
	template void Represent(SparseMatrixTree<d>& hmat,
		const MLO<d>& M, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const TTBasis& basis);

	template void Represent(SparseMatrixTree<d>& hmat, const MLO<d>& M,
		const TensorTree<d>& Psi, const TTBasis& basis);

	template SparseMatrixTree<d> Represent(const MLO<d>& M, const TensorTree<d>& Bra,
		const TensorTree<d>& Ket, const TTBasis& basis);

	template SparseMatrixTree<d> Represent(const MLO<d>& M, const TensorTree<d>& Psi,
		const TTBasis& basis);

}
