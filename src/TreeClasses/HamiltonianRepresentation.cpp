//
// Created by Roman Ellerbrock on 3/3/20.
//
#include "TreeClasses/HamiltonianRepresentation.h"

Tensorcd Apply(const Hamiltonian& H, const Tensorcd& Phi,
	const HamiltonianRepresentation& hRep,
	const Node& node) {

	Tensorcd dPhi(Phi.shape());
	Tensorcd Psi(Phi.shape());
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.hMats_[l];
		if (!hmat.isActive(node)) { continue; }
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.coeff(l) * Phi[i];
		}

		Psi = TreeFunctions::apply(hmat, Psi, H[l], node);

		// Multiply with hole-matrix
		const SparseMatrixTreecd& hhole = hRep.hContractions_[l];
		if (!node.isToplayer() && hhole.isActive(node)) {
			multStateAB(dPhi, hhole[node], Psi, false);
		} else {
			dPhi += Psi;
		}
	}

	if (H.hasV) {
		dPhi += hRep.cdvr_.apply(Phi, hRep.rho_decomposition_[node], node);
//		node.info();
//		cout << "x:\n";
//		auto x = hRep.cdvr_.apply(Phi, hRep.rho_decomposition_[node], node);
//		x.print();
//		cout << "y:\n";
//		auto y = hRep.cdvr_.applySym(Phi, hRep.rho_decomposition_[node], node);
//		y.print();
//		getchar();
	}

	return dPhi;
}

Matrixcd Expectation(const HamiltonianRepresentation& hRep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {

	const Node& top = tree.topNode();
	auto HPhi = Apply(H, Psi[top], hRep, top);
	return Psi[top].dotProduct(HPhi);
}

void LayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

	// dPhi = -i*rho^-1 * (1-P) * <Phi|H|Phi>
	// Apply Hamiltonian
	dPhi = Apply(H, Phi, hRep, node);

	// Inverse Densitymatrix and (1-P) projector for SPF-type EOM
	if (!node.isToplayer()) {
		// (1-P) projector
		dPhi = projectOut(dPhi, Phi);

		// Multiply with inverse single-particle density matrix
		const auto& rhoinv = hRep.rho_inverse_[node];
		dPhi = multStateAB(rhoinv, dPhi);
	}

	dPhi *= -propagation_phase * QM::im;
}

void Derivative(Wavefunction& dPsi, HamiltonianRepresentation& hRep,
	double time, const Wavefunction& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {
	hRep.build(H, Psi, tree, time);
	for (const Node& node: tree) {
		LayerDerivative(dPsi[node], time, Psi[node], H, hRep, node);
	}
}

////////////////////////////////////////////////////////////////////////
///
////////////////////////////////////////////////////////////////////////

void symLayerDerivative(Tensorcd& dPhi, double time, const Tensorcd& Phi,
	const Hamiltonian& H, const HamiltonianRepresentation& hRep,
	const Node& node, complex<double> propagation_phase) {

//	dPhi = TreeFunctions::symApply(Phi, hRep.hMatSets_, H, node);
	dPhi *= -propagation_phase * QM::im;
}

void symDerivative(MatrixTensorTree& dPsi, HamiltonianRepresentation& hRep,
	double time, const MatrixTensorTree& Psi, const Hamiltonian& H,
	const Tree& tree, complex<double> propagation_phase) {

	for (const Node& node: tree) {
		symLayerDerivative(dPsi.first[node], time, Psi.first[node],
			H, hRep, node, propagation_phase);
	}
}

void output(const HamiltonianRepresentation& hrep,
	const Wavefunction& Psi, const Hamiltonian& H, const Tree& tree) {

	cout << "Energies:\n";
	auto h_matrix = Expectation(hrep, Psi, H, tree);
	auto S = TreeFunctions::dotProduct(Psi, Psi, tree);
	auto s_top = S[tree.topNode()];
	for (size_t i = 0; i < h_matrix.dim1(); ++i) {
		double e = real(h_matrix(i, i) / s_top(i, i));
		double e0 = real(h_matrix(0, 0) / s_top(0, 0));
		cout << i << ":\t" << e * QM::cm << " 1/cm\t";
		if (i > 0) { cout << (e - e0) * QM::cm << " \t"; } else { cout << "- g.s. - \t"; }
		cout << e << " a.u." << endl;
	}
}

size_t nActives(const SparseMatrixTreecd& hmat, const MLOcd& M,
	const SparseMatrixTreecd& hcon, const Node& node) {

	/// Count how many neighbors of node are active in hmat & hcon
	size_t n{0};
	if (node.isBottomlayer()) {
		if (M.isActive(node.getLeaf().mode())) { n++; }
	} else {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			if (hmat.isActive(child)) { n++; }
		}
	}

	if (node.isToplayer()) { return n; }
	if (hcon.isActive(node)) { n++; }
	return n;
}

void HamiltonianRepresentation::buildUpCorrection(
	const Hamiltonian& H, const Wavefunction& Psi, const Node& node) {
	/// assert that tail = false for hMats & hCons

	/// reset hCorr matrices
	hCorr_[node].zero();

	/// Add underlying corrections
	if (!node.isBottomlayer()) {
		Tensorcd hA(node.shape());
		for (size_t k = 0; k < node.nChildren(); ++k) {
			const Node& child = node.child(k);
			hA += matrixTensor(hCorr_[child], Psi[node], k);
		}
		hCorr_[node] = contraction(Psi[node], hA, node.parentIdx());
	}

	/// sum h-mats into a single correction matrix
	if (!node.isToplayer()) {
		const Node& parent = node.parent();
		for (size_t l = 0; l < H.size(); ++l) {
//			size_t n_active = nActives(hMats_[l], H[l], hContractions_[l], parent);
			size_t n_active = nActives_[parent][l];
			if (n_active != 1) { continue; } /// only calculate corrections if only 1 neighbor is active
			if (hMats_[l].isActive(node)) {
				hCorr_[node] += H.coeff(l) * hMats_[l][node];
			}
		}
	}
}

void HamiltonianRepresentation::buildDownCorrection(
	const Hamiltonian& H, const Wavefunction& Psi, const Node& hole) {

	/// assert that tail = false for hMats & hCons
	if (hole.isToplayer()) { return; }
	const Node& node = hole.parent();

	/// build upwards correction
	Tensorcd hA(node.shape());
	for (size_t l = 0; l < H.size(); ++l) {
//		size_t n_active = nActives(hMats_[l], H[l], hContractions_[l], node);
		size_t n_active = nActives_[node][l];
		if (n_active < 2) { continue; }
		if (hMats_[l].isActive(hole)) { continue; }
		auto tmp = TreeFunctions::applyHole(hMats_[l], Psi[node], hole);
		if (!node.isToplayer() && hContractions_[l].isActive(node)) {
			tmp = matrixTensor(hContractions_[l][node], tmp, node.parentIdx());
		}
		hA += H.coeff(l) * tmp;
	}

	/// Add underlying corrections
	if (!node.isBottomlayer()) {
		for (size_t k = 0; k < node.nChildren(); ++k) {
			if (k == hole.childIdx()) { continue; }
			const Node& child = node.child(k);
			hA += matrixTensor(hCorr_[child], Psi[node], k);
		}
	}

	if (!node.isToplayer()) {
		hA += matrixTensor(hConCorr_[node], Psi[node], node.parentIdx());
	}

	hConCorr_[hole] = contraction(Psi[node], hA, hole.childIdx());
}

void HamiltonianRepresentation::buildCorrection(const Hamiltonian& H,
	const Wavefunction& Psi, const Tree& tree) {
	for (const Node& node: tree) {
		buildUpCorrection(H, Psi, node);
	}

	for (int n = tree.nNodes() - 1; n >= 0; --n) {
		const Node& hole = tree.getNode(n);
		buildDownCorrection(H, Psi, hole);
	}
}
