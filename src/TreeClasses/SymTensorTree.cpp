//
// Created by Roman Ellerbrock on 1/27/21.
//

#include "TreeClasses/SymTensorTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"

void SymTensorTree::initializeFromTT(const TensorTreecd& Psi, const Tree& tree) {
	/// First set up_ to conventional representation (SPFs orthonormal, Toplayer density-weighted)
	up_ = Psi;
//	QROrthogonal(up_, tree); //

	/// Now create orthonormal down-tensors and simultaneously create weighted tensors
	weighted_ = up_;
	createWeighted(weighted_, down_, up_, tree);

	down_ = up_;
	normalizeCanonical(tree);
}

void SymTensorTree::canonicalRepresentation(const Tree& tree) {
	for (const Node& node : tree) {
		bool up = !node.isToplayer();
		bool down = !node.isBottomlayer();
		weighted_[node] = canonicalTensor(weighted_[node], up, down);
	}
}

void SymTensorTree::normalizeCanonical(const Node& node) {
	const Tensorcd& W = weighted_[node];
	/// Calculate new down-representation
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			Matrixcd w = toMatrix(W, k);
			SVDcd x = svd(w);
			auto U = get<0>(x);
			down_[child] = toTensor(U, W.shape(), k);
		}
	}
	/// Calculate new up-representation
	if (!node.isToplayer()) {
		Matrixcd w = toMatrix(W);
		SVDcd x = svd(w);
		Matrixcd U = get<0>(x);
		up_[node] = toTensor(U, W.shape(), node.parentIdx());
	}

	/// Normalize
	auto S = W.dotProduct(W);
	auto norm = sqrt(abs(S.trace()));
	weighted_[node] /= norm;

}

void SymTensorTree::normalizeCanonical(const Tree& tree) {
	canonicalRepresentation(tree);
	for (const Node& node : tree) {
		normalizeCanonical(node);
	}
}

void SymTensorTree::print(const Tree& tree) {
	for (const Node& node : tree) {
		node.info();
		weighted_[node].print();
	}
}

Tensorcd canonicalTensor(Tensorcd w, bool up, bool down) {
	if (down) {
		for (size_t k = 0; k < (w.shape().order() - 1); ++k) {
			auto rho = contraction(w, w, k);
			auto U = get<0>(svd(rho));
//		auto U = Diagonalize(rho).first;
			w = tensorMatrix(w, U, k);
		}
	}
	if (up) {
		auto rho = contraction(w, w, w.shape().lastIdx());
		auto U = get<0>(svd(rho));
//		auto U = Diagonalize(rho).first;
		w = tensorMatrix(w, U, w.shape().lastIdx());
	}
	return w;
}

Matrixcd calculateB(const Tensorcd& weighted, size_t k) {
	Matrixcd rho = contraction(weighted, weighted, k);
	return toMatrix(sqrt(diagonalize(rho)));
}

template<typename T>
void createWeighted(TensorTree<T>& weighted, TensorTree<T>& down,
	const TensorTreecd& up, const Tree& tree) {
	/** Rationale:
	 * a.) Performs an orthogonalization for Top-Down TensorTree.
	 * b.) Has to be in topdown format, i.e. every Tensor is shifted one layer down.
	 * c.) Meant to be used with MatrixTensorTree.
	 * d.) Not tested and probably not working.
	 */

	for (auto it = tree.rbegin(); it != tree.rend(); ++it) {
		const Node& node = *it;
		if (!node.isToplayer()) {
			const Node& parent = node.parent();

			Matrixcd B = calculateB(weighted[parent], node.childIdx());
			weighted[node] = tensorMatrix(weighted[node], B, node.parentIdx());
//			weighted[node] = MatrixTensor(B, weighted[node], node.nChildren());

		}
	}
}

TensorTreecd SymTensorTree::bottomUpNormalized(const Tree& tree) const {
	TensorTreecd Psi(tree);
	for (const Node& node : tree) {
		Psi[node] = up_[node];
	}
	Psi[tree.TopNode()] = weighted_[tree.TopNode()];
	return Psi;
}

SymTensorTree::SymTensorTree(mt19937& gen, const Tree& tree, bool delta_lowest) {
	TensorTreecd tmp(gen, tree, delta_lowest);
	initializeFromTT(tmp, tree);
}

SymTensorTree::SymTensorTree(mt19937& gen, const TensorTreecd& Psi,
	const SparseTree& stree, const Tree& tree, bool delta_lowest) {
	TensorTreecd tmp(gen, tree, delta_lowest);
	for (const Node& node: tree) {
		if (!stree.isActive(node)) {
			tmp[node] = Psi[node];
		}
	}
	initializeFromTT(tmp, tree);
}

namespace TreeFunctions {

	void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Bra,
		const Tensorcd& Ket, const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hKet = TreeFunctions::applyHole(mat, Ket, hchild);

		if (!parent.isToplayer() && hole.isActive(parent)) {
//			hKet = TensorMatrix(hKet, hole[parent], parent.childIdx());
			hKet = multStateAB(hole[parent], hKet);
		}
		hole[hchild] = contraction(Bra, hKet, hchild.childIdx());
	}

	void symContraction(SparseMatrixTreecd& hole, const TensorTreecd& Bra,
		const TensorTreecd& Ket, const SparseMatrixTreecd& mat, const Tree& tree) {

		int sub_topnode = mat.sparseTree().size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = mat.sparseTree().node(n);
			if (!node.isToplayer()) {
				/// Use Bra/Ket[node] instead of [parent] since Trensors are moved down
				symContractionLocal(hole, Bra[node], Ket[node], mat, node);
			}
		}
	}

	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const MLOcd& M, const Tree& tree) {
		TreeFunctions::represent(mat.first, M, Bra.up_, Ket.up_, tree);
		TreeFunctions::symContraction(mat.second, Bra.down_, Ket.down_, mat.first, tree);
	}

	void symRepresent(SymMatrixTrees& mats, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const SOPcd& S, const Tree& tree) {
		for (size_t l = 0; l < S.size(); ++l) {
			symRepresent(mats[l], Bra, Ket, S[l], tree);
		}
	}

	////////////////////////////////////////////////////////////////////////
	/// apply Operators
	////////////////////////////////////////////////////////////////////////

	Tensorcd symApplyDownNew(const Tensorcd& Phi, const SparseMatrixTreecd& hHole,
		const Node& node) {
		if (node.isToplayer() || !hHole.isActive(node)) { return Phi; }
//		return TensorMatrix(Phi, hHole[node], node.parentIdx()); @TODO: Check if this is correct
		return multStateAB(hHole[node], Phi);
	}

	Tensorcd symApply(const Tensorcd& Ket,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node) {
		Tensorcd hKet = TreeFunctions::apply(mats.first, Ket, M, node);
		return symApplyDownNew(hKet, mats.second, node);
	}

	void symApply(Tensorcd& HPhi, const Tensorcd& Phi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Node& node) {
		HPhi.zero();
		for (size_t l = 0; l < H.size(); ++l) {
			HPhi += H.coeff(l) * symApply(Phi, hmats[l], H[l], node);
		}
	}

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree) {
		const SparseTree& stree = hmats[0].first.sparseTree();

		for (const Node* node_ptr : stree) {
			const Node& node = *node_ptr;
			symApply(HPsi.weighted_[node], Psi.weighted_[node], hmats, H, node);
		}
	}

	void symContraction(MatrixTreecd& rho, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const MatrixTreecd& S, const Tree& tree) {
		for (auto it = tree.rbegin(); it != tree.rend(); ++it) {
			const Node& node = *it;
			if (node.isToplayer()) { continue; }
			const auto& bra = Bra.down_[node];
			auto ket = Ket.down_[node];
			const Node& parent = node.parent();
			size_t child_idx = node.childIdx();
			/// apply S with hole in childidx
			for (size_t k = 0; k < parent.nChildren(); ++k) {
				if (k != child_idx) {
					const Node& child = parent.child(k);
					ket = matrixTensor(S[child], ket, k);
				}
			}
			if (!parent.isToplayer()) {
				ket = multStateAB(rho[parent], ket);
			}
			rho[node] = contraction(bra, ket, node.childIdx());
		}
	}

	Matrixcd DotProduct(const Tensorcd& Bra, Tensorcd Ket,
		const MatrixTreecd& s_up, const MatrixTreecd& s_down, const Node& node) {
		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				Ket = matrixTensor(s_up[child], Ket, k);
			}
		}
		if (!node.isToplayer()) {
			Ket = multStateAB(s_down[node], Ket);
//			Ket = MatrixTensor(s_down[node], Ket, node.parentIdx());
		}
		return Bra.dotProduct(Ket);
	}

	MatrixTreecd DotProduct(const TensorTreecd& wBra, const TensorTreecd& wKet,
		const MatrixTreecd& s_up, const MatrixTreecd& s_down, const Tree& tree) {
		MatrixTreecd S(tree);
		for (const Node& node : tree) {
			S[node] = DotProduct(wBra[node], wKet[node], s_up, s_down, node);
		}
		return S;
	}

	MatrixTreecd symDotProduct(const SymTensorTree& Bra, const SymTensorTree& Ket, const Tree& tree) {
		/**
		 * Rationale:
		 * - Calculates dot-products in symmetric representation. First calculate up-, then down-overlaps
		 *   (s_up, s_down). Use these to calculate overall overlaps S. If everything works correctly, all
		 *   nodes should provide equal outputs abs(S[node_a]) == abs(S[node_b])
		 */
		MatrixTreecd s_up = TreeFunctions::dotProduct(Bra.up_, Ket.up_, tree);
		MatrixTreecd s_down(tree);
		symContraction(s_down, Bra, Ket, s_up, tree);
		MatrixTreecd S(tree);
		for (const Node& node : tree) {
			S[node] = DotProduct(Bra.weighted_[node], Ket.weighted_[node], s_up, s_down, node);
		}
		return S;
	}
}

double epsUpNormalization(const SymTensorTree& Psi, const Tree& tree) {
	/// Check Bottom-up
	double r = 0.;
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.up_[node];
		auto S = Phi.dotProduct(Phi);
		r += abs(residual(S, identityMatrixcd(S.dim1())));
	}
	return r;
}

double epsDownNormalization(const SymTensorTree& Psi, const Tree& tree) {
	/// Check top-down
	double r = 0.;
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.down_[node];
		auto S = contraction(Phi, Phi, node.childIdx());
		r += abs(residual(S, identityMatrixcd(S.dim1())));
	}
	return r;
}

bool isWorking(const SymTensorTree& Psi, const Tree& tree) {
	bool works = true;
	double eps = 1e-7;
	double r_up = epsUpNormalization(Psi, tree);
	double r_down = epsDownNormalization(Psi, tree);
	if (r_up > eps) {
		cerr << "parent-normalization failed.\n";
		works = false;
	}
	if (r_down > eps) {
		cerr << "Down-normalization failed.\n";
		works = false;
	}
	return works;
}

