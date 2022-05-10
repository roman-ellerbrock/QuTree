//
// Created by Roman Ellerbrock on 1/27/21.
//

#include "TreeClasses/SymTensorTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Core/Tensor_Extension.h"

void SymTensorTree::initialize(const Tree& tree) {
	weighted_.initialize(tree);
	up_.initialize(tree);
	down_.initialize(tree);
}

SymTensorTree::SymTensorTree(mt19937& gen, const Tree& tree, bool delta_lowest)
	: SymTensorTree() {
	initialize(tree);

	for (const Node& node : tree) {
		Tensor_Extension::generate(weighted_[node], gen);

		Tensorcd& A = weighted_[node];
		if (delta_lowest) {
			for (size_t i = 0; i < node.shape().lastBefore(); ++i) {
				A(i) = 0;
			}
			A(0) = 1.;
		}
		gramSchmidt(A);
	}
	normalizeUp(tree);
	normalizeDown(tree);
}

void SymTensorTree::orthogonalUp(const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const Tensorcd& A = weighted_[node];
			size_t k = node.parentIdx();

			Tensorcd& Q = up_[node];
			Q = qr(A, k);

			/// rotate parent
			auto R = contraction(A, Q, k);
			const Node& parent = node.parent();
			Tensorcd& Apar = weighted_[parent];
			Apar = matrixTensor(R, Apar, node.childIdx());
		}
	}
}

void SymTensorTree::orthogonalDown(const Tree& tree) {
	for (auto it = tree.rbegin(); it != tree.rend(); ++it) {
		const Node& node = *it;
		if (node.isToplayer()) { continue; }
		const Node& parent = node.parent();
		Tensorcd& A = weighted_[parent];
		size_t k = node.childIdx();

		Tensorcd& Q = down_[node];
		Q = qr(A, k);

		auto R = contraction(A, Q, k);
		Tensorcd& Apar = weighted_[node];
		Apar = matrixTensor(R, Apar, node.parentIdx());
	}
}

void SymTensorTree::orthogonal(const Tree& tree) {
	orthogonalUp(tree);
	orthogonalDown(tree);
}

SymTensorTree::SymTensorTree(TensorTreecd Psi, const Tree& tree) {
	initialize(tree);
	for (const Node& node : tree) {
		weighted_[node] = Psi[node];
		if (!node.isToplayer()) {
			up_[node] = Psi[node];
			down_[node] = Psi[node.parent()];
		}
	}

	orthogonal(tree);

}

void SymTensorTree::normalizeUp(const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			up_[node] = Tensor_Extension::normalize(weighted_[node], node.parentIdx(), eps_);
		}
	}

}

void SymTensorTree::normalizeDown(const Tree& tree) {
	for (const Node& node : tree) {
		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			down_[node] = Tensor_Extension::normalize(weighted_[parent], node.childIdx(), eps_);
		}
	}
}

void SymTensorTree::normalize(const Tree& tree) {
	normalizeUp(tree);
	normalizeDown(tree);
}

namespace TreeFunctions {

	///////////////////////////////////////////////////////////////////////
	/// Tree Contractions
	///////////////////////////////////////////////////////////////////////

	void contractionUpLocal(MatrixTreecd& S, const Tensorcd& Bra, Tensorcd Ket,
		const Node& node) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			if (!node.isBottomlayer()) {
				const Node& child = node.child(k);
				Ket = matrixTensor(S[child], Ket, k);
			}
		}
		S[node] = Bra.dotProduct(Ket);
	}

	void contractionUp(MatrixTreecd& S, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const Tree& tree) {
		for (const Node& node : tree) {
			contractionUpLocal(S, Bra.up_[node], Ket.up_[node], node);
		}
	}

	void contractionDownLocal(MatrixTreecd& Sdown, const Tensorcd& Bra, Tensorcd Ket,
		const MatrixTreecd& S, const Node& child) {
		const Node& node = child.parent();
		assert(!node.isToplayer());
		size_t hole = child.childIdx();
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& ochild = node.child(k);
			if (k != hole) {
				Ket = matrixTensor(S[ochild], Ket, k);
			}
		}
		if (!node.isToplayer()) {
			const Node& parent = node.parent();
			Ket = matrixTensor(Sdown[parent], Ket, node.parentIdx());
		}
		Sdown[child] = contraction(Bra, Ket, hole);
	}

	void contractionDown(MatrixTreecd& Sdown, const SymTensorTree& Bra,
		const SymTensorTree& Ket, const MatrixTreecd& S, const Tree& tree) {
		for (int i = tree.nNodes() - 1; i >= 0; --i) {
			const Node& node = tree.getNode(i);
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				contractionDownLocal(Sdown, Bra.down_[parent], Ket.down_[parent], S, node);
			}
		}
	}

	Tensorcd symApply(Tensorcd Ket, const MatrixTreecd& Sup, const MatrixTreecd& Sdown,
		const Node& node) {
		if (!node.isBottomlayer()) {
			for (size_t k = 0; k < node.nChildren(); ++k) {
				const Node& child = node.child(k);
				Ket = matrixTensor(Sup[child], Ket, k);
			}
		}
		if (!node.isToplayer()) {
			Ket = matrixTensor(Sdown[node.parent()], Ket, node.parentIdx());
		}
		return Ket;
	}

	 vector<double> dotProduct(const SymTensorTree& Bra, SymTensorTree Ket,
		const Tree& tree) {
		MatrixTreecd Sup(tree); // wazzzz suuuuuuuup???!!!
		MatrixTreecd Sdown(tree);
		contractionUp(Sup, Bra, Ket, tree);
		contractionDown(Sdown, Bra, Ket, Sup, tree);

		vector<double> eps;
		for (const Node& node : tree) {
			Tensorcd SKet = symApply(Ket.weighted_[node], Sup, Sdown, node);
			Matrixcd s = Bra.weighted_[node].dotProduct(SKet);
			eps.push_back(abs(s.trace()));
		}
		return eps;
	}

	double error(const vector<double>& vec) {
		double val = 0.;
		for (auto x : vec) {
			val += pow(x, 2);
		}
		return sqrt(val);
	}

	///////////////////////////////////////////////////////////////////////
	/// Tree Functions
	///////////////////////////////////////////////////////////////////////

	void symContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Bra,
		const Tensorcd& Ket, const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hKet = TreeFunctions::applyHole(mat, Ket, hchild);

		if (!parent.isToplayer() && hole.isActive(parent)) {
			//			hKet = TensorMatrix(hKet, hole[parent], parent.childIdx());
			hKet = multStateAB(hole[parent], hKet);
			//			hKet = matrixTensor(hole[parent], hKet, hchild.parentIdx());
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
	/**
	 * Idea:
	 * - always perform normalize up right before building up-/down-matrices
	 * - only perform once per operator to avoid multiple evaluations
	 * - weighted Tensors store the important information
	 * - allows to perform different 'normalization' for CDVR
	 * 		- this can incorporate finding optimized x-basis sets
	 */

	/**
	 *
	 * @param HPsi
	 * @param Psi
	 * @param hmats
	 * @param H
	 * @param node
	 */
	void iterate(SymTensorTree& HPsi, SymMatrixTrees& mats, const SymTensorTree& Psi,
		const SOPcd& H, const Tree& tree) {

		symApply(HPsi, Psi, mats, H, tree);
		HPsi.normalize(tree);
		symRepresent(mats, HPsi, Psi, H, tree);

	}

	void ApplySCF(SymTensorTree& HPsi, SymMatrixTrees& mats, const SymTensorTree& Psi,
		const SOPcd& H, const Tree& tree, double eps, size_t max_iter, ostream* os) {

		auto HPsi_last = HPsi;
		for (size_t i = 0; i < max_iter; ++i) {
			iterate(HPsi, mats, Psi, H, tree);
			double change = error(dotProduct(HPsi_last, HPsi, tree));
			if (os) {
				*os << "i = " << i << "\t" << change << endl;
			}
			if (change < eps) { break; }
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

		for (const Node *node_ptr : stree) {
			const Node& node = *node_ptr;
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

