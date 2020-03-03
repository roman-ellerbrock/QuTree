//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#define SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#include "TreeClasses/SparseMatrixTreeFunctions.h"

namespace SparseMatrixTreeFunctions {
////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Bottom-parent (Forward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Matrix<T> RepresentUpper(const SparseMatrixTree<T>& hmat,
		const Tensor<T>& Bra, const Tensor<T>& Ket, const Node& node) {
		Tensor<T> hKet(Ket);
		for (size_t l = 0; l < node.nChildren(); l++) {
			const Node& child = node.child(l);
			if (!hmat.Active(child)) { continue; }
			hKet = MatrixTensor(hmat[child], hKet, child.childIdx());
		}

		return Bra.DotProduct(hKet);
	}

	template<typename T>
	Matrix<T> RepresentBottom(const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node, const Leaf& leaf) {
		Tensor<T> MKet = M.ApplyBottomLayer(Ket, leaf);
		return Bra.DotProduct(MKet);
	}

	template<typename T>
	void RepresentLayer(SparseMatrixTree<T>& mats, const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node) {
		if (!mats.Active(node)) { return; }

		if (node.isBottomlayer()) {
			mats[node] = RepresentBottom(Bra, Ket, M, node, node.getLeaf());
		} else {
			mats[node] = RepresentUpper(mats, Bra, Ket, node);
		}
	}

	template<typename T>
	void Represent(SparseMatrixTree<T>& hmat,
		const MLO<T>& M, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const Tree& tree) {
		assert(Bra.size() == Ket.size());
		const SparseTree& active = hmat.Active();
		for (size_t n = 0; n < active.size(); ++n) {
			const Node& node = active.MCTDHNode(n);
			RepresentLayer(hmat, Bra[node], Ket[node], M, node);
		}
	}

	template<typename T>
	void Represent(SparseMatrixTree<T>& hmat, const MLO<T>& M,
		const TensorTree<T>& Psi, const Tree& tree) {
		Represent(hmat, M, Psi, Psi, tree);
	}

	template<typename T>
	SparseMatrixTree<T> Represent(const MLO<T>& M, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const Tree& tree) {
		SparseMatrixTree<T> hmat(M, tree);
		Represent(hmat, M, Bra, Ket, tree);
		return hmat;
	}

	template<typename T>
	SparseMatrixTree<T> Represent(const MLO<T>& M, const TensorTree<T>& Psi,
		const Tree& tree) {
		SparseMatrixTree<T> hmat(M, tree);
		Represent(hmat, M, Psi, tree);
		return hmat;
	}

	template<typename T>
	void Represent(SparseMatrixTrees<T>& Mats, const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		assert(Mats.size() == sop.size());
		for (size_t l = 0; l < sop.size(); ++l) {
			Represent(Mats[l], sop[l], Bra, Ket, tree);
		}
	}

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Top-child (Backward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const SparseTree& marker, const Tree& tree) {

		// Swipe top-down_ but exclude topnode
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.MCTDHNode(n);
			if (!node.isToplayer()) {
				assert(mats.Active(node));
				assert(holes.Active(node));

				const Node& parent = node.parent();
				Tensor<T> hKet = ApplyHole(mats, Ket[parent], node);
				if (!parent.isToplayer()) {
					hKet = multStateAB(holes[parent], hKet);
				}
				holes[node] = mHoleProduct(Bra[parent], hKet, node.childIdx());
			}
		}
	}

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const Tree& tree) {
		Contraction(holes, Bra, Ket, mats, holes.Active(), tree);
	}

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const Tree& tree) {
		Contraction(holes, Psi, Psi, mats, tree);
	}

	template<typename T>
	void Contraction(vector<SparseMatrixTree<T>>& holes, const SparseMatrixTrees<T>& Mats,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		assert(holes.size() == Mats.size());
		for (size_t l = 0; l < holes.size(); ++l) {
			Contraction(holes[l], Bra, Ket, Mats[l], tree);
		}
	}

////////////////////////////////////////////////////////////////////////
/// Apply SparseMatrixTree to tensor tree
////////////////////////////////////////////////////////////////////////

	/// Apply factor matrices locally
	template<typename T>
	Tensor<T> Apply(const SparseMatrixTree<T>& mats, const Tensor<T>& Phi,
		const MLO<T>& M, const Node& node) {
		if (!mats.Active(node)) { return Phi; }
		if (node.isBottomlayer()) {
			const Leaf& phys = node.getLeaf();
			return M.ApplyBottomLayer(Phi, phys);
		} else {
			return ApplyUpper(mats, Phi, node);
		}
	}

	template<typename T>
	Tensor<T> ApplyUpper(const SparseMatrixTree<T>& mat, Tensor<T> Phi, const Node& node) {
		Tensor<T> hPhi(Phi.shape());
		bool switchbool = true;
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (!mat.Active(child)) { continue; }
			if (switchbool) {
				MatrixTensor(hPhi, mat[child], Phi, child.childIdx(), true);
			} else {
				MatrixTensor(Phi, mat[child], hPhi, child.childIdx(), true);
			}
			switchbool = !switchbool;
		}
		if (switchbool) {
			return Phi;
		} else {
			return hPhi;
		}
	}

	template<typename T>
	Tensor<T> ApplyHole(const SparseMatrixTree<T>& mats, Tensor<T> Phi, const Node& hole_node) {
		assert(!hole_node.isToplayer());
		const Node& parent = hole_node.parent();
		size_t drop = hole_node.childIdx();

		for (size_t k = 0; k < parent.nChildren(); ++k) {
			const Node& child = parent.child(k);
			size_t childidx = child.childIdx();
			if ((childidx == drop) || (!mats.Active(child))) { continue; }
			Phi = MatrixTensor(mats[child], Phi, childidx);
		}
		return Phi;
	}

}

#endif //SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
