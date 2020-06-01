//
// Created by Roman Ellerbrock on 5/21/20.
//

#include "TreeOperators/TreeStructured/MatrixListTree.h"


void MatrixListTree::Initialize(const TreeSOP& op, const Tree& tree) {
	attributes_.resize(tree.nNodes());
}

void MatrixListTree::Calculate(const TensorTreecd& Psi,
	const TreeSOP& op, const Tree& tree) {
	/// Calculate the MatrixListTree at all layers
	for (const Node& node : tree) {
		operator[](node) = Calculate(Psi[node], op, node);
	}
}

MatrixList MatrixListTree::Calculate(const Tensorcd& Phi,
	const TreeSOP& H, const Node& node) {
	/// Calculate the MatrixListTree for a NodeSOPlist at an upper layer
	MatrixList spolist;
	for (const NodeSOP& S : H[node]) {
		spolist.emplace_back(Calculate(Phi, S, H, node));
	}
	return spolist;
}

Matrixcd MatrixListTree::Calculate(const Tensorcd& Phi, const NodeSOP& S,
	const TreeSOP& H, const Node& node) {
	/// Calculate SPO of the the MatrixListTree for a NodeSOP at an upper layer
	Tensorcd SPhi = Apply(Phi, S, H, node);
	Matrixcd s = Phi.DotProduct(SPhi);
	return s;
//	return SPOcd(s, node.childIdx());
}

Tensorcd MatrixListTree::Apply(const Tensorcd& Phi, const NodeSOP& S,
	const TreeSOP& H, const Node& node) const {
	/// Calculate SPO of the the MatrixListTree for a NodeSOP at an upper layer
	Tensorcd SPhi(Phi.shape());
	for (size_t l = 0; l < S.size(); ++l) {
		auto coeff = S.Coeff(l);
		// M = x1*x1
		const NodeProductOperator& M = S[l];
		Tensorcd MPhi = ApplyNodeProductOperator(Phi, M, H, node);
		MPhi *= coeff;
		SPhi += MPhi;
	}
	return SPhi;
}

Tensorcd MatrixListTree::ApplyNodeProductOperator(const Tensorcd& Phi,
	const NodeProductOperator& M, const TreeSOP& H,
	const Node& node) const {
	if (node.isBottomlayer()) {
		return ApplyNodeProductOperator_bottom(Phi, M, H.leafOperatorLib(node), node);
	} else {
		return ApplyNodeProductOperator_upper(Phi, M, node);
	}
}

Tensorcd MatrixListTree::ApplyNodeProductOperator_upper(Tensorcd hPhi,
	const NodeProductOperator& M, const Node& node) const {
	for (const NodeOperator& hl : M) {
		const Node& child = node.child(hl.Mode());
		const MatrixList& hs = operator[](child);
		const Matrixcd& h = hs[hl.Part()];
		hPhi = MatrixTensor(h, hPhi, child.childIdx());
//		hPhi = h * hPhi;
	}
	return hPhi;
}

Tensorcd MatrixListTree::ApplyNodeProductOperator_bottom(Tensorcd Phi,
	const NodeProductOperator& M, const LeafOperatorLib& lib,
	const Node& node) const {
	const Leaf& leaf = node.getLeaf();
	const LeafInterface& grid = leaf.PrimitiveGrid();
	Tensorcd hPhi(Phi.shape());
	for (const NodeOperator& h : M) {
		const shared_ptr<LeafOperatorcd>& b_h = lib[h.Part()];
		b_h->Apply(grid, hPhi, Phi);
		Phi = hPhi;
	}
	return hPhi;
}

Tensorcd MatrixListTree::ApplyHole(const Tensorcd& Phi, const NodeSOP& S,
	const Node& node, const NodeOperator& hole_S) const {

	Tensorcd SPhi(Phi.shape());
	for (size_t L = 0; L < S.size(); ++L) {
		const NodeProductOperator& M = S[L];
		if (IsActive(hole_S, M)) {
			Tensorcd MPhi = ApplyHole(Phi, M, node, hole_S);
			MPhi *= S.Coeff(L);
			SPhi += MPhi;
		}
	}
	return SPhi;
}

Tensorcd MatrixListTree::ApplyHole(Tensorcd MPhi, const NodeProductOperator& M,
	const Node& node, const NodeOperator& hole_S) const {
	for (const NodeOperator& hl : M) {
		if (!(hl == hole_S)) {
			const Node& child = node.child(hl.Mode());
			const MatrixList& hs = operator[](child);
			const Matrixcd& h = hs[hl.Part()];
//			MPhi = h * MPhi;
			MPhi = MatrixTensor(h, MPhi, child.childIdx());
		}
	}
	return MPhi;
}

void MatrixListTree::print(const Tree& tree, ostream& os) const {

	for (const Node& node : tree) {
		node.info();
		const MatrixList& hs = operator[](node);
		for (const Matrixcd& h : hs) {
			h.print();
		}
	}
}

