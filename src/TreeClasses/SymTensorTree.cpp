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
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			Matrixcd w = toMatrix(W, k);
			SVDcd x = svd(w);
			auto U = get<0>(x);
			down_[child] = toTensor(U, W.shape(), k);
		}
	}
	if (!node.isToplayer()) {
		Matrixcd w = toMatrix(W);
		SVDcd x = svd(w);
		Matrixcd U = get<0>(x);
		up_[node] = toTensor(U, W.shape(), node.parentIdx());
	}
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
			auto rho = Contraction(w, w, k);
			auto U = get<0>(svd(rho));
//		auto U = Diagonalize(rho).first;
			w = TensorMatrix(w, U, k);
		}
	}
	if (up) {
		auto rho = Contraction(w, w, w.shape().lastIdx());
		auto U = get<0>(svd(rho));
//		auto U = Diagonalize(rho).first;
		w = TensorMatrix(w, U, w.shape().lastIdx());
	}
	return w;
}

Matrixcd calculateB(const Tensorcd& weighted, size_t k) {
	Matrixcd rho = Contraction(weighted, weighted, k);
	return toMatrix(sqrt(Diagonalize(rho)));
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
			weighted[node] = TensorMatrix(weighted[node], B, node.parentIdx());
//			weighted[node] = MatrixTensor(B, weighted[node], node.nChildren());

		}
	}
}

TensorTreecd SymTensorTree::bottomUpNormalized(const Tree& tree) {
	TensorTreecd Psi(tree);
	for (const Node& node : tree) {
		Psi[node] = up_[node];
	}
	Psi[tree.TopNode()] = weighted_[tree.TopNode()];
	return Psi;
}

namespace TreeFunctions {

	void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Bra,
		const Tensorcd& Ket, const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hKet = TreeFunctions::ApplyHole(mat, Ket, hchild);

		if (!parent.isToplayer() && hole.Active(parent)) {
			hKet = TensorMatrix(hKet, hole[parent], parent.childIdx());
		}
		hole[hchild] = Contraction(Bra, hKet, hchild.childIdx());
	}

	void symContraction(SparseMatrixTreecd& hole, const TensorTreecd& Bra,
		const TensorTreecd& Ket, const SparseMatrixTreecd& mat, const Tree& tree) {

		int sub_topnode = mat.Active().size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = mat.Active().MCTDHNode(n);
			if (!node.isToplayer()) {
				symContractionLocal(hole, Bra[node], Ket[node], mat, node);
			}
		}
	}

	void symRepresent(SymMatrixTree& mat, const SymTensorTree& Bra,
		const SymTensorTree& Ket,
		const MLOcd& M, const Tree& tree) {
		TreeFunctions::Represent(mat.first, M, Bra.up_, Ket.up_, tree);
		TreeFunctions::symContraction(mat.second, Bra.down_, Ket.down_, mat.first, tree);
	}

	void symRepresent(SymMatrixTrees& mats, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const SOPcd& S, const Tree& tree) {
		for (size_t l = 0; l < S.size(); ++l) {
			symRepresent(mats[l], Bra, Ket, S[l], tree);
		}
	}

	////////////////////////////////////////////////////////////////////////
	/// Apply Operators
	////////////////////////////////////////////////////////////////////////

	Tensorcd symApply(const Tensorcd& Ket,
		const SymMatrixTree& mats, const MLOcd& M, const Node& node) {
		Tensorcd hKet = TreeFunctions::Apply(mats.first, Ket, M, node);
		return symApplyDown(hKet, mats.second, node);
	}

	void symApply(Tensorcd& HPhi, const Tensorcd& Phi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Node& node) {
		HPhi.Zero();
		for (size_t l = 0; l < H.size(); ++l) {
			HPhi += H.Coeff(l) * symApply(Phi, hmats[l], H[l], node);
		}
	}

	void symApply(SymTensorTree& HPsi, const SymTensorTree& Psi,
		const SymMatrixTrees& hmats,
		const SOPcd& H, const Tree& tree) {
		for (const Node& node : tree) {
			symApply(HPsi.weighted_[node], Psi.weighted_[node], hmats, H, node);
		}
	}
}

double epsUpNormalization(const SymTensorTree& Psi, const Tree& tree) {
	/// Check Bottom-up
	double r = 0.;
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.up_[node];
		auto S = Phi.DotProduct(Phi);
		r += abs(Residual(S, IdentityMatrixcd(S.Dim1())));
	}
	return r;
}

double epsDownNormalization(const SymTensorTree& Psi, const Tree& tree) {
	/// Check top-down
	double r = 0.;
	for (const Node& node : tree) {
		if (node.isToplayer()) { continue; }
		const auto& Phi = Psi.down_[node];
		auto S = Contraction(Phi, Phi, node.childIdx());
		r += abs(Residual(S, IdentityMatrixcd(S.Dim1())));
	}
	return r;
}

bool isWorking(const SymTensorTree& Psi, const Tree& tree) {
	bool works = true;
	double eps = 1e-7;
	double r_up = epsUpNormalization(Psi, tree);
	double r_down = epsDownNormalization(Psi, tree);
	if (r_up > eps) {
		cerr << "Up-normalization failed.\n";
		works = false;
	}
	if (r_down > eps) {
		cerr << "Down-normalization failed.\n";
		works = false;
	}
	return works;
}

