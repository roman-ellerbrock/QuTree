//
// Created by Roman Ellerbrock on 2/12/20.
//
#include "TreeClasses/SparseMatrixTreeFunctions_Implementation.h"

namespace TreeFunctions {
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

	template void Represent<cd>(SparseMatrixTrees<cd>& Mats, const SOP<cd>& sop,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void Represent(SOPMatrixTrees<cd>& mats, const SOP<cd>& sop,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void RepresentLayer(SparseMatrixTree<cd>& mats, const Tensor<cd>& Bra,
		const Tensor<cd>& Ket, const MLO<cd>& M, const Node& node);

	template void Contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void Contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Psi,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void Contraction(vector<SparseMatrixTree<cd>>& holes, const SparseMatrixTrees<cd>& Mats,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void Contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const SparseMatrixTree<cd>& mats, const MatrixTree<cd>& rho, const Tree& tree);

	template void Contraction(SparseMatrixTrees<cd>& holes, const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const SparseMatrixTrees<cd>& Mats,
		const MatrixTree<cd>& rho, const Tree& tree);

	template SparseMatrixTree<cd> Contraction(const TensorTree<cd>& Bra,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void Contraction<cd>(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi,
		const SparseTree& stree, bool orthogonal);

	template Tensor<cd> Apply(const SparseMatrixTree<cd>& mat, const Tensor<cd>& Phi, const MLO<cd>& M, const Node& node);

	template Tensor<cd> ApplyUpper(const SparseMatrixTree<cd>& mat, Tensor<cd> Phi, const Node& node);

	template Tensor<cd> ApplyHole(const SparseMatrixTree<cd>& holes, Tensor<cd> Phi, const Node& hole_node);


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

	template void Represent<d>(SparseMatrixTrees<d>& Mats, const SOP<d>& sop,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void Represent(SOPMatrixTrees<d>& mats, const SOP<d>& sop,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void RepresentLayer(SparseMatrixTree<d>& mats, const Tensor<d>& Bra,
		const Tensor<d>& Ket, const MLO<d>& M, const Node& node);


	template void Contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template void Contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template SparseMatrixTree<d> Contraction(const TensorTree<d>& Bra,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template void Contraction(vector<SparseMatrixTree<d>>& holes, const SparseMatrixTrees<d>& Mats,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void Contraction<d>(MatrixTree<d>& Rho, const TensorTree<d>& Psi,
		const SparseTree& stree, bool orthogonal);

	template Tensor<d> Apply(const SparseMatrixTree<d>& mat, const Tensor<d>& Phi, const MLO<d>& M, const Node& node);

	template Tensor<d> ApplyUpper(const SparseMatrixTree<d>& mat, Tensor<d> Phi, const Node& node);

	template Tensor<d> ApplyHole(const SparseMatrixTree<d>& holes, Tensor<d> Phi, const Node& hole_node);

}
