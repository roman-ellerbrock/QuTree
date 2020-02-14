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

		for (auto it = basis.end(); it >= basis.begin(); it--) {
			const Node& node = *it;
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
/// Apply SparseMatrixTree to tensor tree
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Tensor<T> ApplyUpper(const SparseMatrixTree<T>& mats, Tensor<T> Phi, const Node& node) {
		Tensor<T> hPhi(Phi.Dim());
		bool switchbool = true;
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.Down(k);
			if (!mats.Active(child)) { continue; }
			if (switchbool) {
				multAB(hPhi, mats[child], Phi, true);
			} else {
				multAB(Phi, mats[child], hPhi, true);
			}
			switchbool = !switchbool;
		}
		if (switchbool) {
			return Phi;
		} else {
			return hPhi;
		}
	}

	/// Apply factor matrices locally
	template<typename T>
	Tensor<T> Apply(const SparseMatrixTree<T>& mats, const Tensor<T>& Phi,
		const MLO<T>& M, const Node& node) {
		if (!mats.Active(node)) { return Phi; }
		if (node.IsBottomlayer()) {
			const Leaf& phys = node.PhysCoord();
			return M.ApplyBottomLayer(Phi, phys);
		} else {
			return ApplyUpper(mats, Phi, node);
		}
	}

	template<typename T>
	Tensor<T> ApplyHole(const SparseMatrixTree<T>& mats, Tensor<T> Phi, const Node& hole_node) {
		assert(!hole_node.IsToplayer());
		const Node& parent = hole_node.Up();
		size_t drop = hole_node.ChildIdx();

		for (size_t k = 0; k < parent.nChildren(); ++k) {
			const Node& child = parent.Down(k);
			size_t childidx = child.ChildIdx();
			if ((childidx == drop) || (!mats.Active(child))) { continue; }
			Phi = multAB(mats[child], Phi, childidx);
		}
		return Phi;
	}

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Top-Down (Backward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats,const TreeMarker& marker, const TTBasis& basis) {

		// Swipe top-down_ but exclude topnode
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.MCTDHNode(n);
			if (!node.IsToplayer()) {
				assert(mats.Active(node));
				assert(holes.Active(node));

				const Node& parent = node.Up();
				Tensor<T> hKet = ApplyHole(mats, Ket[parent], node);
				if (!parent.IsToplayer()) {
					hKet = multStateAB(holes[parent], hKet);
				}
				holes[node] = HoleProduct(Bra[parent], hKet, node.ChildIdx());
			}
		}
	}

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const TTBasis& basis) {
		Contraction(holes, Bra, Ket, mats, holes.Active(), basis);
	}

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const TTBasis& basis) {
		Contraction(holes, Psi, Psi, mats, basis);
	}
}

#endif //SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
