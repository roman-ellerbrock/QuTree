//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#define SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#include "Tree/SparseMatrixTreeFunctions.h"

namespace SparseMatrixTreeFunctions {
////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Bottom-Up (Forward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Matrix<T> RepresentUpper(const SparseMatrixTree<T>& hmat,
		const Tensor<T>& Bra, const Tensor<T>& Ket, const Node& node) {
		// @TODO: Optimize with switchbool trick
		Tensor<T> hKet(Ket);
		for (size_t l = 0; l < node.nChildren(); l++) {
			const Node& child = node.Down(l);
			if (!hmat.Active(child)) { continue; }
			hKet = multAB(hmat[child], hKet, node.ChildIdx()); // @TODO: CHECK THIS; product and getter
		}

		return Bra.DotProduct(hKet);
	}

	template<typename T>
	Matrix<T> RepresentBottom(const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node, const Leaf& phys) {
		Tensor<T> MKet = M.ApplyBottomLayer(Ket, phys);
		return Bra.DotProduct(MKet);
	}

	template<typename T>
	void RepresentLayer(SparseMatrixTree<T>& mats, const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node) {
		if (!mats.Active(node)) { return; }

		if (node.IsBottomlayer()) {
			mats[node] = RepresentBottom(Bra, Ket, M, node, node.PhysCoord());
		} else {
			mats[node] = RepresentUpper(mats, Bra, Ket, node);
		}
	}

	template<typename T>
	void Represent(SparseMatrixTree<T>& hmat,
		const MLO<T>& M, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const TTBasis& basis) {
		assert(Bra.size() == Ket.size());
		const TreeMarker& active = hmat.Active();
		for (size_t n = 0; n < active.size(); ++n) {
			const Node& node = active.MCTDHNode(n);
			RepresentLayer(hmat, Bra[node], Ket[node], M, node);
		}
	}

	template<typename T>
	void Represent(SparseMatrixTree<T>& hmat, const MLO<T>& M,
		const TensorTree<T>& Psi, const TTBasis& basis) {
		Represent(hmat, M, Psi, Psi, basis);
	}

	template<typename T>
	SparseMatrixTree<T> Represent(const MLO<T>& M, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const TTBasis& basis) {
		SparseMatrixTree<T> hmat(M, basis);
		Represent(hmat, M, Bra, Ket, basis);
		return hmat;
	}

	template<typename T>
	SparseMatrixTree<T> Represent(const MLO<T>& M, const TensorTree<T>& Psi,
		const TTBasis& basis) {
		SparseMatrixTree<T> hmat(M, basis);
		Represent(hmat, M, Psi, basis);
		return hmat;
	}

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Top-Down (Backward)
////////////////////////////////////////////////////////////////////////

}

#endif //SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
