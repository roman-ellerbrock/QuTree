//
// Created by Roman Ellerbrock on 1/3/21.
//

#include "TreeClasses/MatrixTensorTreeFunctions.h"

namespace TreeFunctions {

	void Represent(SparseMatrixTreecd& mat, const MatrixTensorTree& Psi, const MLOcd& M,
		const Tree& tree) {

		const TensorTreecd& Psi_up = Psi.BottomUpNormalized(tree);
		TreeFunctions::Represent(mat, M, Psi_up, tree);
	}

	void ContractionLocal(SparseMatrixTreecd& hole, const Tensorcd& Phi,
		const SparseMatrixTreecd& mat, const Node& hchild) {
		assert(!hchild.isToplayer());
		const Node& parent = hchild.parent();

		Tensorcd hPhi = TreeFunctions::ApplyHole(hole, Phi, hchild);

		if (!parent.isToplayer() && hole.Active(parent)) {
			hPhi = tensorMatrix(hPhi, hole[parent], parent.childIdx());
		}
		hole[hchild] = contraction(Phi, hPhi, hchild.childIdx());
	}

	void Contraction(SparseMatrixTreecd& hole, const MatrixTensorTree& Psi,
		const SparseMatrixTreecd& mat, const SparseTree& marker, const Tree& tree) {

		const TensorTreecd& Psi_down = Psi.TopDownNormalized(tree);
		int sub_topnode = marker.size() - 1;
		for (int n = sub_topnode; n >= 0; --n) {
			const Node& node = marker.MCTDHNode(n);
			if (!node.isToplayer()) {
				ContractionLocal(hole, Psi_down[node], mat, node);
			}
		}

	}

	void Represent(SparseMatrixTreePaircd& mats,
		const MatrixTensorTree& Psi, const MLOcd& M,
		const Tree& tree) {

		Represent(mats.first, Psi, M, tree);
		Contraction(mats.second, Psi, mats.first, mats.second.Active(), tree);

	}

	void Represent(SparseMatrixTreePairscd& matset, const MatrixTensorTree& Psi,
		const SOPcd& H, const Tree& tree) {

		for (size_t l = 0; l < H.size(); ++l) {
			Represent(matset[l], Psi, H[l], tree);
		}

	}

	Tensorcd symApplyDown(const Tensorcd& Phi, const SparseMatrixTreecd& hHole,
		const Node& node) {
		if (node.isToplayer() || !hHole.Active(node)) { return Phi; }
		return tensorMatrix(Phi, hHole[node], node.parentIdx());
		/// @TODO: WAS THIS USED?? For Applying operators, this is wrong!
	}

/*	Tensorcd symApply(const Tensorcd& Phi,
		const SparseMatrixTreePaircd& mats, const MLOcd& M, const Node& node) {
		Tensorcd hPhi = TreeFunctions::Apply(mats.first, Phi, M, node);
		return symApplyDown(hPhi, mats.second, node);
	}

	void symApply(MatrixTensorTree& HPsi, const MatrixTensorTree& Psi,
		const SparseMatrixTreePaircd& hpair, const MLOcd& M, complex<double> c, const Tree& tree) {
		const SparseMatrixTreecd& hup = hpair.first;
		const SparseTree& stree = hup.Active();
		for (const Node* node_ptr : stree) {
			const Node& node = *node_ptr;
			HPsi.nodes()[node] += c * symApply(Psi.nodes()[node], hpair, M, node);
		}
	}

	void symApply(MatrixTensorTree& HPsi, const MatrixTensorTree& Psi,
		const SparseMatrixTreePairscd& hMats,
		const SOPcd& H, const Tree& tree) {

		for (size_t i = 0; i < H.size(); ++i) {
			/// @TODO: Coefficients?
			symApply(HPsi, Psi, hMats[i], H[i], H.Coeff(i), tree);
		}
	}
 */

/*
	Tensorcd symApply(Tensorcd Phi,
		const SparseMatrixTreePairscd& hMatSet,
		const SOPcd& H, const Node& node) {

		Tensorcd hPhi(Phi.shape());
		for (size_t l = 0; l < H.size(); ++l) {
			Tensorcd Psi = H.Coeff(l) * Phi;
			hPhi += symApply(Psi, hMatSet[l], H[l], node);
		}
		return hPhi;
	}

	void symApply(MatrixTensorTree& Chi, const MatrixTensorTree& Psi,
		const SparseMatrixTreePairscd& hMats,
		const SOPcd& H, const Tree& tree) {

		if (hMats.empty()) {
			cerr << "Cannot Apply h-Matrices: Empty Operator Representation\n";
			exit(1);
		}
		const SparseTree& stree = hMats.front().first.Active();
		for (const Node* node_ptr : stree) {
			const Node& node = *node_ptr;
			Chi.first[node] = symApply(Psi.first[node], hMats, H, node);
		}
	}
*/
}

