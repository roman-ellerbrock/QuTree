//
// Created by Roman Ellerbrock on 2/12/20.
//
#include "TreeClasses/SparseMatrixTreeFunctions_Implementation.h"

namespace SparseMatrixTreeFunctions {
	typedef complex<double> cd;

	template void Represent(SparseMatrixTree<cd>& hmat,
		const MLO<cd>& M, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const Tree& tree);

	template void Represent(SparseMatrixTree<cd>& hmat, const MLO<cd>& M,
		const TensorTree<cd>& Psi, const Tree& tree);

	template SparseMatrixTree<cd> Represent(const MLO<cd>& M, const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const Tree& tree);

	template SparseMatrixTree<cd> Represent(const MLO<cd>& M, const TensorTree<cd>& Psi,
		const Tree& tree);

	template Tensor<cd> Apply(const SparseMatrixTree<cd>& mats, const Tensor<cd>& Phi,
		const MLO<cd>& M, const Node& node);

	template void Contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void Contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Psi,
		const SparseMatrixTree<cd>& mats, const Tree& tree);


	typedef double d;
	template void Represent(SparseMatrixTree<d>& hmat,
		const MLO<d>& M, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const Tree& tree);

	template void Represent(SparseMatrixTree<d>& hmat, const MLO<d>& M,
		const TensorTree<d>& Psi, const Tree& tree);

	template SparseMatrixTree<d> Represent(const MLO<d>& M, const TensorTree<d>& Bra,
		const TensorTree<d>& Ket, const Tree& tree);

	template SparseMatrixTree<d> Represent(const MLO<d>& M, const TensorTree<d>& Psi,
		const Tree& tree);

	template Tensor<d> Apply(const SparseMatrixTree<d>& mats, const Tensor<d>& Phi,
		const MLO<d>& M, const Node& node);

	template void Contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template void Contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra,
		const SparseMatrixTree<d>& mats, const Tree& tree);
}