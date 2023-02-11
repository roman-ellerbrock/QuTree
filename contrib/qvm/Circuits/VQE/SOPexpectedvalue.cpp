//
// Created by Roman Ellerbrock on 8/20/21.
//
#include "TreeClasses/SOPMatrixTrees.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"

Tensorcd Apply(const SOPcd& H, const Tensorcd& Phi,
	const MatrixTreescd& hRep,
	const Node& node) {

	Tensorcd HPhi(Phi.shape());
	Tensorcd Psi(Phi.shape());
	/// For every term in the Hamiltonian
	for (size_t l = 0; l < H.size(); l++) {
		const auto& hmat = hRep.matrices_[l];
		if (!hmat.isActive(node)) { continue; }

		/// Initialize Psi via |Psi> = c_l * |Phi>
		for (size_t i = 0; i < Psi.shape().totalDimension(); ++i) {
			Psi[i] = H.coeff(l) * Phi[i];
		}

		/// Apply matrix representation of H for every term l
		Psi = TreeFunctions::apply(hmat, Psi, H[l], node);

		/// Multiply with Mean-field matrix
		const auto& hhole = hRep.contractions_[l];
		if (!node.isToplayer() && hhole.isActive(node)) {
			multStateAB(HPhi, hhole[node], Psi, false);
		} else {
			HPhi += Psi;
		}
	}

	return HPhi;
}

Matrixcd expectation(const SOPcd& H,
	const TensorTreecd& Psi, const Tree& tree) {

	MatrixTreescd hRep(H, tree);
	TreeFunctions::represent(hRep, H, Psi, Psi, tree);

	const Node& top = tree.topNode();
	auto HPhi = Apply(H, Psi[top], hRep, top);
	return Psi[top].dotProduct(HPhi);
}

double expectationVal(const SOPcd& H,
	const TensorTreecd& Psi, const Tree& tree) {
	return real(expectation(H, Psi, tree)(0,0));
}