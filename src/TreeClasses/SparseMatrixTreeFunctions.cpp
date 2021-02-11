//
// Created by Roman Ellerbrock on 2/12/20.
//
#include "TreeClasses/SparseMatrixTreeFunctions_Implementation.h"

namespace TreeFunctions {
	typedef complex<double> cd;

	template void represent(SparseMatrixTree<cd>& hmat,
		const MLO<cd>& M, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const Tree& tree);

	template void represent(SparseMatrixTree<cd>& hmat, const MLO<cd>& M,
		const TensorTree<cd>& Psi, const Tree& tree);

	template SparseMatrixTree<cd> represent(const MLO<cd>& M, const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const Tree& tree);

	template SparseMatrixTree<cd> represent(const MLO<cd>& M, const TensorTree<cd>& Psi,
		const Tree& tree);

	template void represent<cd>(SparseMatrixTrees<cd>& Mats, const SOP<cd>& sop,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void represent(SOPMatrixTrees<cd>& mats, const SOP<cd>& sop,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void representLayer(SparseMatrixTree<cd>& mats, const Tensor<cd>& Bra,
		const Tensor<cd>& Ket, const MLO<cd>& M, const Node& node);

	template SparseMatrixTrees<cd> represent(const SOP<cd>& sop,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		shared_ptr<SparseTree>& stree, const Tree& tree);


	template void contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Psi,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template void contraction(vector<SparseMatrixTree<cd>>& holes, const SparseMatrixTrees<cd>& Mats,
		const TensorTree<cd>& Bra, const TensorTree<cd>& Ket, const Tree& tree);

	template void contraction(SparseMatrixTree<cd>& holes, const TensorTree<cd>& Bra, const TensorTree<cd>& Ket,
		const SparseMatrixTree<cd>& mats, const MatrixTree<cd>& rho, const Tree& tree);

	template void contraction(SparseMatrixTrees<cd>& holes, const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const SparseMatrixTrees<cd>& Mats,
		const MatrixTree<cd>& rho, const Tree& tree);

	template SparseMatrixTree<cd> contraction(const TensorTree<cd>& Bra,
		const SparseMatrixTree<cd>& mats, const Tree& tree);

	template vector<SparseMatrixTree<cd>> contraction(const TensorTree<cd>& Bra,
		const TensorTree<cd>& Ket, const vector<SparseMatrixTree<cd>>& mats,
		const MatrixTree<cd>& rho, shared_ptr<SparseTree>& stree, const Tree& tree);

	template void contraction<cd>(MatrixTree<cd>& Rho, const TensorTree<cd>& Psi,
		const SparseTree& stree, bool orthogonal);

	template Tensor<cd> apply(const SparseMatrixTree<cd>& mat, const Tensor<cd>& Phi, const MLO<cd>& M, const Node& node);

	template Tensor<cd> applyUpper(const SparseMatrixTree<cd>& mat, Tensor<cd> Phi, const Node& node);

	template Tensor<cd> applyHole(const SparseMatrixTree<cd>& holes, Tensor<cd> Phi, const Node& hole_node);


	typedef double d;
	template void represent(SparseMatrixTree<d>& hmat,
		const MLO<d>& M, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const Tree& tree);

	template void represent(SparseMatrixTree<d>& hmat, const MLO<d>& M,
		const TensorTree<d>& Psi, const Tree& tree);

	template SparseMatrixTree<d> represent(const MLO<d>& M, const TensorTree<d>& Bra,
		const TensorTree<d>& Ket, const Tree& tree);

	template SparseMatrixTree<d> represent(const MLO<d>& M, const TensorTree<d>& Psi,
		const Tree& tree);

	template void represent<d>(SparseMatrixTrees<d>& Mats, const SOP<d>& sop,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void represent(SOPMatrixTrees<d>& mats, const SOP<d>& sop,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template void representLayer(SparseMatrixTree<d>& mats, const Tensor<d>& Bra,
		const Tensor<d>& Ket, const MLO<d>& M, const Node& node);

	template SparseMatrixTrees<d> represent(const SOP<d>& sop,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		shared_ptr<SparseTree>& stree, const Tree& tree);


	template void contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra, const TensorTree<d>& Ket,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template void contraction(SparseMatrixTree<d>& holes, const TensorTree<d>& Bra,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template SparseMatrixTree<d> contraction(const TensorTree<d>& Bra,
		const SparseMatrixTree<d>& mats, const Tree& tree);

	template void contraction(vector<SparseMatrixTree<d>>& holes, const SparseMatrixTrees<d>& Mats,
		const TensorTree<d>& Bra, const TensorTree<d>& Ket, const Tree& tree);

	template vector<SparseMatrixTree<d>> contraction(const TensorTree<d>& Bra,
		const TensorTree<d>& Ket, const vector<SparseMatrixTree<d>>& mats,
		const MatrixTree<d>& rho, shared_ptr<SparseTree>& stree, const Tree& tree);

	template void contraction<d>(MatrixTree<d>& Rho, const TensorTree<d>& Psi,
		const SparseTree& stree, bool orthogonal);

	template Tensor<d> apply(const SparseMatrixTree<d>& mat, const Tensor<d>& Phi, const MLO<d>& M, const Node& node);

	template Tensor<d> applyUpper(const SparseMatrixTree<d>& mat, Tensor<d> Phi, const Node& node);

	template Tensor<d> applyHole(const SparseMatrixTree<d>& holes, Tensor<d> Phi, const Node& hole_node);

}
