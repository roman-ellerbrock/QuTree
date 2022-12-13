//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "DVR/MatrixTensorTree.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "Core/TensorBLAS.h"

MatrixTensorTree::MatrixTensorTree(const Wavefunction& Psi,
	const Tree& tree, bool orthogonal) {
	initialize(Psi, tree, orthogonal);
}

void MatrixTensorTree::initialize(const Wavefunction& Psi,
	const Tree& tree, bool orthogonal) {

	/// Note: Requires orthogonal wavefunction representation (typically given)
	assert(orthogonal);

	buildNodes(Psi, tree);
	buildEdges(tree);
}

void MatrixTensorTree::buildNodes(TensorTreecd Psi, const Tree& tree) {
	orthogonal(Psi, tree);
	nodes() = Psi;
}

void MatrixTensorTree::buildEdges(const Tree& tree) {
	Wavefunction& nodes_ = first;
	MatrixTreecd& edges_ = second;

	/// Get edge matrices (B's)
	MatrixTreecd rho = TreeFunctions::contraction(nodes_, tree, true);
	auto B = sqrt(rho, tree);

	/// Build node representation (A^\tilde's)
	for (const Edge& e : tree.edges()) {
		const Node& node = e.down();
		Tensorcd& A = nodes_[node];
		/// Basically like multiplying with sqrt(rho)'s
		A = matrixTensorBLAS(B[e], A, node.nChildren());
	}

	edges_ = inverse(B, tree);
}

void MatrixTensorTree::buildFromWeighted(const Tree& tree) {
	/// Re-orthonormalize bottom-up
	buildNodes(bottomUpNormalized(tree), tree);
	/// TopDown Orthonormal
	buildEdges(tree);
}

TensorTreecd MatrixTensorTree::topDownNormalized(const Tree& tree) const {
	/// Contraction-normalized representation
	/// Build A^{(p\circ k) p}
	/// Note: Tensors get moved one layer down!
	TensorTreecd Psi(nodes());
	for (const Edge& e : tree.edges()) {
		const Node& node = e.down();
		const Node& parent = node.parent();
		Psi[node] = matrixTensorBLAS(edges()[e].transpose(), Psi[parent], node.childIdx());
	}

	return Psi;
}

TensorTreecd MatrixTensorTree::bottomUpNormalized(const Tree& tree) const {
	/// Dot-Product normalized reperesentation
	/// Build A^{p\circ k (p)}
	/// This is the conventional wavefunction representation.
	TensorTreecd Psi(nodes());

	for (const Edge& e : tree.edges()) {
		const Node& node = e.down();
		Psi[node] = matrixTensorBLAS(edges()[e], Psi[node], node.nChildren());
	}

	return Psi;
}

bool IsWorking_bottomup(const MatrixTensorTree& Psi, const Tree& tree, double eps) {
	auto bottomup = Psi.bottomUpNormalized(tree);
	for (const Node& node : tree) {
		auto x = contraction(bottomup[node], bottomup[node], node.nChildren());
		auto r = residual(x, identityMatrixcd(x.dim1()));
		if (r > eps) {
			cerr << "bottom-up normalization failed.\n";
			return false;
		}
	}
	return true;
}

bool IsWorking_topdown(const MatrixTensorTree& Psi, const Tree& tree, double eps) {
	auto topdown = Psi.topDownNormalized(tree);
	for (const Edge& e : tree.edges()) {
		const Node& node = e.down();
		auto x = contraction(topdown[node], topdown[node], node.childIdx());
		auto r = residual(x, identityMatrixcd(x.dim1()));
		if (r > eps) {
			cerr << "top-down normalization failed.\n";
			return false;
		}
	}
	return true;
}

bool IsWorking(const MatrixTensorTree& Psi, const Tree& tree, double eps) {
	return (IsWorking_bottomup(Psi, tree, eps) && IsWorking_topdown(Psi, tree, eps));
}
