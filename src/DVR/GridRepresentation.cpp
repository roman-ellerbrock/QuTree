//
// Created by Roman Ellerbrock on 1/2/21.
//
#include "GridRepresentation.h"
#include "DVR/MatrixTensorTree.h"

namespace OptimizeGrid{
	typedef pair<Matrixcd, size_t> MatrixIdx;
	Matrixcd BuildX(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node) {

		if (node.isBottomlayer()) {

			const Leaf& leaf = node.getLeaf();
			const LeafInterface& grid = leaf.PrimitiveGrid();
			auto w = 1. / abs(rho.Trace()) * rho;
			w = Regularize(w, 1e-5);
			Tensorcd xPhi(Phi.shape());
			grid.applyX(xPhi, Phi);
			auto wxPhi = MatrixTensor(w, xPhi, node.parentIdx());
			return Tensor_Extension::OuterProduct(wxPhi, xPhi);
		} else {

		}

	}

	Tensorcd OptimizeUp(const Tensorcd& Phi, const Matrixcd& rho,
		const Node& node, const Node& node_small) {

		size_t n_occ = node_small.shape().lastDimension();

		auto X = BuildX(Phi, rho, node);
		X = UnProject(n_occ, X, Phi);
		auto xspec = Diagonalize(X);
		auto oPhi = Occupy(Phi, xspec.first, n_occ, node);

		return oPhi;
	}

	ExplicitEdgeWavefunction OptimizeUp(ExplicitEdgeWavefunction Psi,
		const Tree& tree, const Tree& tree_small) {

		ExplicitEdgeWavefunction Chi(Psi);
		for (const Node& node : tree) {
			const Node& node2 = tree_small.GetNode(node.Address());
			if (!node.isToplayer() && (node.shape() != node2.shape())) {
				Chi[node] = Optimize(Psi[node], rho[node],
					node, node2);
				Update(Chi, tree);
			}
		}
		return Chi;
	}


}