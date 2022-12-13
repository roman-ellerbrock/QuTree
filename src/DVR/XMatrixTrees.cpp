//
// Created by Roman Ellerbrock on 3/11/20.
//
#include "DVR/XMatrixTrees.h"
#include <random>

Wavefunction Regularize(Wavefunction Psi, const Tree& tree, double eps) {

	normal_distribution<double> dist;
	mt19937 gen(23949);
	for (const Node& node : tree) {
		Tensorcd& A = Psi[node];
		for (size_t i = 0; i < A.shape().totalDimension(); ++i) {
			A(i) += eps * dist(gen);
		}
		gramSchmidt(A);
	}
	return Psi;
}

//a.) introduce scale factors
//b.) Use explicit edge representation

SOPcd Xsop(const Tree& tree) {
	LeafFuncd x = &LeafInterface::applyX;
	LeafFuncd I = &LeafInterface::identity;
	SOPcd xops;
	for (size_t l = 0; l < tree.nLeaves(); ++l) {
		const Leaf& leaf = tree.getLeaf(l);
		size_t mode = leaf.mode();
		MLOcd M(x, mode);
		for (size_t i = 0; i < tree.nLeaves(); ++i) {
			M.push_back(I, i);
		}
		xops.push_back(M, 1.);
	}
	return xops;
}

Tensorcd XMatrixTrees::Optimize(const Tensorcd& Phi, const Matrixcd& rho,
	const Node& node, const Node& node_small) const {

	size_t n_occ = node_small.shape().lastDimension();

	/// x rho x
	auto X = BuildX(Phi, rho, mats_, node, eps_);
	/// (1-P) x rho x ( 1-P)
	X = UnProject(n_occ, X, Phi);
	auto xspec = diagonalize(X);
	auto oPhi = Occupy(Phi, xspec.first, n_occ, node);

	return oPhi;
}

Wavefunction XMatrixTrees::Optimize(Wavefunction Psi,
	const MatrixTreecd& rho, const Tree& tree, const Tree& tree_small) {

	Wavefunction Chi(Psi);
	for (const Node& node : tree) {
		const Node& node2 = tree_small.getNode(node.address());
		if (!node.isToplayer() && (node.shape() != node2.shape())) {
			Chi[node] = Optimize(Psi[node], rho[node],
				node, node2);
			Update(Chi, tree);
		}
	}
	return Chi;
}

Matrixcd BuildX(const Tensorcd& Phi, const Matrixcd& rho,
	const SparseMatrixTreescd& xmats, const Node& node, double eps) {

	if (node.isBottomlayer()) {
		const Leaf& leaf = node.getLeaf();
		const LeafInterface& grid = leaf.interface();
		auto w = 1. / abs(rho.trace()) * rho;
		w = regularize(w, eps);
		Tensorcd xPhi(Phi.shape());
		grid.applyX(xPhi, Phi);
		auto wxPhi = matrixTensor(w, xPhi, node.parentIdx());
		return Tensor_Extension::outerProduct(wxPhi, xPhi);
	} else {
		size_t bef = node.shape().lastBefore();
		Matrixcd X(bef, bef);
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			for (const SparseMatrixTreecd& x : xmats) {
				if (x.isActive(child)) {
					Tensorcd xPhi = matrixTensor(x[child], Phi, k);
					// Build needed Tensors
					Tensor_Extension::weightedOuterProductAdd(X, xPhi, xPhi, rho);
				}
			}
		}
		return X;
	}
}

Matrixcd UnProject(size_t n_occupied, const Matrixcd& X,
	const Tensorcd& Phi) {
	/**
	 * Calculate (1 - P) A (1 - P) from A
	 */
	size_t dimpart = Phi.shape().lastBefore();
	size_t nlast = Phi.shape().lastDimension();
	assert (Phi.shape().lastDimension() >= n_occupied);

	// Build unit
	Matrixcd Pu(dimpart, dimpart);
	for (size_t i = 0; i < dimpart; ++i)
		Pu(i, i) = 1;

	// Build P (it doesnt matter if you use ntesor or small_ntensor
	size_t start = nlast - n_occupied;
	for (size_t i = 0; i < dimpart; ++i) {
		for (size_t j = 0; j < dimpart; ++j) {
			for (size_t n = start; n < nlast; ++n) {
				Pu(j, i) -= Phi(j, n) * conj(Phi(i, n));
			}
		}
	}
	/*
	// Build P (it doesnt matter if you use ntesor or small_ntensor
	for (size_t i = 0; i < dimpart; ++i) {
		for (size_t j = 0; j < dimpart; ++j) {
			for (size_t n = 0; n < n_occupied; ++n) {
				Pu(j, i) -= Phi(j, n) * conj(Phi(i, n));
			}
		}
	}
	*/
	// (1 - P) A (1 - P)
	Matrixcd B = X * Pu;
	return Pu * B;
}

Tensorcd Occupy(const Tensorcd& Phi, const Matrixcd& trafo,
	size_t n_occupied, const Node& node) {

	const TensorShape& shape = node.shape();
	size_t dimpart = shape.lastBefore();
	size_t ntensor = shape.lastDimension();
	assert(ntensor >= n_occupied);

	// Copy the highest value eigenvectors to a new Phi
	Tensorcd oPhi(Phi);

	// (The following ints have to be integers, not size_t)
	int last = (int) dimpart - 1;
	for (int n = 0; n < (ntensor - n_occupied); ++n) {
		int mat_idx = last - n;
		int t_idx = (int) n;
		for (size_t i = 0; i < dimpart; i++) {
			oPhi(i, (size_t) t_idx) = trafo(i, (size_t) mat_idx);
		}
	}

/*	int last = (int) dimpart - 1;
	for (int n = 0; n < (ntensor - n_occupied); ++n) {
		int mat_idx = last - n;
		int t_idx = (int) n_occupied + n;
		for (size_t i = 0; i < dimpart; i++) {
			oPhi(i, (size_t) t_idx) = trafo(i, (size_t) mat_idx);
		}
	}*/

//	GramSchmidt(oPhi);
	return oPhi;
}

////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////

