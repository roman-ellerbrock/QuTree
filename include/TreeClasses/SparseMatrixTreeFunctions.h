//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef SPARSEMATRIXTREEFUNCTIONS_H
#define SPARSEMATRIXTREEFUNCTIONS_H
#include "TreeClasses/SparseMatrixTree.h"
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/MatrixTree.h"

namespace TreeFunctions {
/**
 * \namespace SparseMatrixTreeFunctions
 *
 * \ingroup Tree
 *
 * \brief This class provides functions to build SparseMatrixTrees.
 *
 * SparseMatrixTrees typically appear when representing operators
 * in a tensortree basis. Similar to MatrixTree, a SparseMatrixTree
 * can be build bottom-up or top-down in a tree. A SparsematrixTree is
 * build bottom-up if an operator is represented and top-down if
 * the representation of the operator is represented.
 * */

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Bottom-parent (Forward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Matrix<T> representUpper(const SparseMatrixTree<T>& hmat,
		const Tensor<T>& Bra, const Tensor<T>& Ket, const Node& node, Tensor<T>* work = nullptr);

	template<typename T>
	void representLayer(SparseMatrixTree<T>& mats, const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node, Tensor<T>* work = nullptr);

	template<typename T>
	void represent(SparseMatrixTree<T>& hmat, const MLO<T>& M,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const Tree& tree);

	template<typename T>
	SparseMatrixTree<T> represent(const MLO<T>& M,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const Tree& tree);

	template<typename T>
	void represent(SparseMatrixTree<T>& hmat, const MLO<T>& M,
		const TensorTree<T>& Psi, const Tree& tree);

	template<typename T>
	SparseMatrixTree<T> represent(const MLO<T>& M,
		const TensorTree<T>& Psi, const Tree& tree);

	template <typename T>
	void represent(vector<SparseMatrixTree<T>>& Mats, const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree);

	template<typename T>
	SparseMatrixTrees<T> represent(const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		shared_ptr<SparseTree>& stree, const Tree& tree);

	template <typename T>
	void represent(SOPMatrixTrees<T>& mats, const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree);

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Top-child (Backward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	SparseMatrixTree<T> contraction(const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const Tree& tree);

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const SparseMatrixTree<T>& mats,
		const Tree& tree);

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const SparseMatrixTree<T>& mats,
		const MatrixTree<T>& rho, const SparseTree& marker, const Tree& tree);

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const MatrixTree<T>& rho, const Tree& tree);

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const Tree& tree);

	template <typename T>
	void contraction(SparseMatrixTrees<T>& holes, const SparseMatrixTrees<T>& mat,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree);

	template <typename T>
	void contraction(SparseMatrixTrees<T>& holes, const SparseMatrixTrees<T>& mat,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree);

	template <typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const MatrixTree<T>& rho, const Tree& tree);

	template <typename T>
	void contraction(SparseMatrixTrees<T>& holes, const SparseMatrixTrees<T>& mat,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree);

	template <typename T>
	void contraction(SparseMatrixTrees<T>& holes, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const SparseMatrixTrees<T>& mats,
		const MatrixTree<T>& rho, const Tree& tree);

	template <typename T>
	vector<SparseMatrixTree<T>> contraction(const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const vector<SparseMatrixTree<T>>& mats,
		const MatrixTree<T>& rho, shared_ptr<SparseTree>& stree, const Tree& tree);

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi,
		const SparseTree& stree, bool orthogonal = true);

////////////////////////////////////////////////////////////////////////
/// apply MatrixTree
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Tensor<T> apply(const SparseMatrixTree<T>& mat, const Tensor<T>& Phi,
		const MLO<T>& M, const Node& node);

	template<typename T>
	void apply(Tensor<T>& hPhi, const SparseMatrixTree<T>& mat,
		const SparseMatrixTree<T> *holes, const MatrixTree<T> *rho,
		Tensor<T> Phi, const SparseTree& stree, const Node& node, int skip, Tensor<T>* work = nullptr);

	template<typename T>
	Tensor<T> apply(const SparseMatrixTree<T>& mat,
		const SparseMatrixTree<T>& holes, const MatrixTree<T> *rho,
		Tensor<T> Phi, const SparseTree& stree, const Node& node, int skip);

	template<typename T>
	Tensor<T> applyUpper(const SparseMatrixTree<T>& mat, Tensor<T> Phi, const Node& node);

	template<typename T>
	Tensor<T> applyHole(const SparseMatrixTree<T>& holes, Tensor<T> Phi, const Node& hole_node);

	template<typename T>
	void apply(Tensor<T>& hPhi, const SparseMatrixTree<T>& mat, const Tensor<T>& Phi,
		const MLO<T>& M, const Node& node);
}

#endif //SPARSEMATRIXTREEFUNCTIONS_H
