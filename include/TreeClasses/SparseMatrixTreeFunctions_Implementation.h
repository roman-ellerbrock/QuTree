//
// Created by Roman Ellerbrock on 2/12/20.
//

#ifndef SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#define SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeClasses/MatrixTreeFunctions.h"

namespace TreeFunctions {
////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Bottom-parent (Forward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	Matrix<T> representUpper(const SparseMatrixTree<T>& hmat,
		const Tensor<T>& Bra, const Tensor<T>& Ket, const Node& node) {
		Tensor<T> hKet(Ket);
		for (size_t l = 0; l < node.nChildren(); l++) {
			const Node& child = node.child(l);
			if (!hmat.isActive(child)) { continue; }
			hKet = matrixTensor(hmat[child], hKet, child.childIdx());
		}

		return Bra.dotProduct(hKet);
	}

	template<typename T>
	Matrix<T> RepresentBottom(const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node, const Leaf& leaf) {
		Tensor<T> MKet = M.apply(Ket, leaf);
		return Bra.dotProduct(MKet);
	}

	template<typename T>
	void representLayer(SparseMatrixTree<T>& mats, const Tensor<T>& Bra,
		const Tensor<T>& Ket, const MLO<T>& M, const Node& node) {
		if (!mats.isActive(node)) { return; }

		if (node.isBottomlayer()) {
			mats[node] = RepresentBottom(Bra, Ket, M, node, node.getLeaf());
		} else {
			mats[node] = representUpper(mats, Bra, Ket, node);
		}
	}

	template<typename T>
	void represent(SparseMatrixTree<T>& hmat,
		const MLO<T>& M, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const Tree& tree) {
		assert(Bra.size() == Ket.size());
		const SparseTree& active = hmat.sparseTree();
		for (size_t n = 0; n < active.size(); ++n) {
			const Node& node = active.node(n);
			if (!node.isToplayer()) {
				representLayer(hmat, Bra[node], Ket[node], M, node);
			}
		}
	}

	template<typename T>
	void represent(SparseMatrixTree<T>& hmat, const MLO<T>& M,
		const TensorTree<T>& Psi, const Tree& tree) {
		represent(hmat, M, Psi, Psi, tree);
	}

	template<typename T>
	SparseMatrixTree<T> represent(const MLO<T>& M, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const Tree& tree) {
		SparseMatrixTree<T> hmat(M, tree);
		represent(hmat, M, Bra, Ket, tree);
		return hmat;
	}

	template<typename T>
	SparseMatrixTree<T> represent(const MLO<T>& M, const TensorTree<T>& Psi,
		const Tree& tree) {
		SparseMatrixTree<T> hmat(M, tree);
		represent(hmat, M, Psi, tree);
		return hmat;
	}

	template<typename T>
	void represent(SparseMatrixTrees<T>& Mats, const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		assert(Mats.size() == sop.size());
		for (size_t l = 0; l < sop.size(); ++l) {
			represent(Mats[l], sop[l], Bra, Ket, tree);
		}
	}

	template<typename T>
	SparseMatrixTrees<T> represent(const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		shared_ptr<SparseTree>& stree, const Tree& tree) {

		vector<SparseMatrixTree<T>> Mats;
		for (size_t l = 0; l < sop.size(); ++l) {
			SparseMatrixTree<T> M(stree, tree);
			represent(M, sop[l], Bra, Ket, tree);
			Mats.push_back(M);
		}
		return Mats;
	}


	template<typename T>
	void represent(SOPMatrixTrees<T>& mats, const SOP<T>& sop,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		represent(mats.matrices_, sop, Bra, Ket, tree);
		contraction(mats.contractions_, mats.matrices_, Bra, Ket, tree);
	}

////////////////////////////////////////////////////////////////////////
/// Build SparseMatrixTree Top-down (Backward)
////////////////////////////////////////////////////////////////////////

	template<typename T>
	SparseMatrixTree<T> contraction(const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const Tree& tree) {
		SparseMatrixTree<T> Con(mats);
		contraction(Con, Psi, mats, tree);
		return Con;
	}

	template<typename T>
	void Contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const SparseTree& marker, const Tree& tree) {

		// Swipe top-down_ but exclude topnode
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.node(n);
			if (!node.isToplayer()) {
				assert(holes.isActive(node));

				const Node& parent = node.parent();
				Tensor<T> hKet = applyHole(mats, Ket[parent], node);
				if (!parent.isToplayer()) {
					if (!marker.isActive(parent)) {
						cerr << "Error in Contraction of operator representation:\n";
						cerr << "Missing active node at parent.\n";
						exit(1);
					}
				}
				if (marker.isActive(parent)) {
					hKet = multStateAB(holes[parent], hKet);
				}
				holes[node] = contraction(Bra[parent], hKet, node.childIdx());
			} else {
				holes[node] = identityMatrix<T>(node.shape().lastDimension());
			}
		}
	}

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const MatrixTree<T>& rho,
		const SparseTree& marker, const Tree& tree) {

		// Swipe top-down_ but exclude topnode
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.node(n);
			if (!node.isToplayer()) {
				assert(holes.isActive(node));

				const Node& parent = node.parent();
				Tensor<T> hKet = applyHole(mats, Ket[parent], node);
				if (!parent.isToplayer()) {
					if (marker.isActive(parent)) {
						hKet = multStateAB(holes[parent], hKet);
					} else {
						hKet = multStateAB(rho[parent], hKet);
					}
				}
				holes[node] = contraction(Bra[parent], hKet, node.childIdx());
			} else {
				holes[node] = identityMatrix<T>(node.shape().lastDimension());
			}
		}
	}

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const Tree& tree) {
		Contraction(holes, Bra, Ket, mats, holes.sparseTree(), tree);
	}

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Bra, const TensorTree<T>& Ket,
		const SparseMatrixTree<T>& mats, const MatrixTree<T>& rho, const Tree& tree) {
		contraction(holes, Bra, Ket, mats, rho, holes.sparseTree(), tree);
	}

	template<typename T>
	void contraction(SparseMatrixTree<T>& holes, const TensorTree<T>& Psi,
		const SparseMatrixTree<T>& mats, const Tree& tree) {
		contraction(holes, Psi, Psi, mats, tree);
	}

	template<typename T>
	void contraction(vector<SparseMatrixTree<T>>& holes, const SparseMatrixTrees<T>& Mats,
		const TensorTree<T>& Bra, const TensorTree<T>& Ket, const Tree& tree) {
		assert(holes.size() == Mats.size());
		for (size_t l = 0; l < holes.size(); ++l) {
			contraction(holes[l], Bra, Ket, Mats[l], tree);
		}
	}

	template<typename T>
	void contraction(SparseMatrixTrees<T>& holes, const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const SparseMatrixTrees<T>& mats,
		const MatrixTree<T>& rho, const Tree& tree) {
		assert(holes.size() == mats.size());
		for (size_t l = 0; l < holes.size(); ++l) {
			contraction(holes[l], Bra, Ket, mats[l], rho, tree);
		}
	}

	template <typename T>
	vector<SparseMatrixTree<T>> contraction(const TensorTree<T>& Bra,
		const TensorTree<T>& Ket, const vector<SparseMatrixTree<T>>& mats,
		const MatrixTree<T>& rho, shared_ptr<SparseTree>& stree, const Tree& tree) {
		vector<SparseMatrixTree<T>> holes;
		for (const auto& mat : mats) {
			holes.emplace_back(SparseMatrixTree<T>(stree, tree));
		}
		contraction(holes, Bra, Ket, mats, rho, tree);
		return holes;
	}

	template<typename T>
	void contraction(MatrixTree<T>& Rho, const TensorTree<T>& Psi,
		const SparseTree& stree, bool orthogonal) {
		if (!orthogonal) {
			cerr << "SparseTree contraction not implemented for non-orthogonal basis sets.\n";
			exit(1);
		}
		for (int i = stree.size() - 1; i > 0; --i) {
			const Node& node = stree.node(i);
			if (!node.isToplayer()) {
				const MatrixTree<T> *null = nullptr;
				const Node& parent = node.parent();
				contractionLocal(Rho, Psi[parent], Psi[parent], node, null);
			}
		}
	}

////////////////////////////////////////////////////////////////////////
/// apply SparseMatrixTree to tensor tree
////////////////////////////////////////////////////////////////////////

	/// apply factor matrices locally
	template<typename T>
	Tensor<T> apply(const SparseMatrixTree<T>& mats, const Tensor<T>& Phi,
		const MLO<T>& M, const Node& node) {
		if (!mats.isActive(node)) { return Phi; }
		if (node.isBottomlayer()) {
			const Leaf& phys = node.getLeaf();
			return M.apply(Phi, phys);
		} else {
			return applyUpper(mats, Phi, node);
		}
	}

	template<typename T>
	Tensor<T> applyUpper(const SparseMatrixTree<T>& mat, Tensor<T> Phi, const Node& node) {
		Tensor<T> hPhi(Phi.shape());
		bool switchbool = true;
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (!mat.isActive(child)) { continue; }
			if (switchbool) {
				matrixTensor(hPhi, mat[child], Phi, child.childIdx(), true);
			} else {
				matrixTensor(Phi, mat[child], hPhi, child.childIdx(), true);
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
	Tensor<T> applyHole(const SparseMatrixTree<T>& mats, Tensor<T> Phi, const Node& hole_node) {
		assert(!hole_node.isToplayer());
		const Node& parent = hole_node.parent();
		size_t drop = hole_node.childIdx();

		for (size_t k = 0; k < parent.nChildren(); ++k) {
			const Node& child = parent.child(k);
			size_t childidx = child.childIdx();
			if ((childidx == drop) || (!mats.isActive(child))) { continue; }
			Phi = matrixTensor(mats[child], Phi, childidx);
		}
		return Phi;
	}
}

#endif //SPARSEMATRIXTREEFUNCTIONS_IMPLEMENTATION_H
