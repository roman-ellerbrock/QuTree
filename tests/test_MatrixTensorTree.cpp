//
// Created by Roman Ellerbrock on 4/3/20.
//

#include "UnitTest++/UnitTest++.h"
#include "TreeClasses//MatrixTensorTree.h"
#include "TreeShape/TreeFactory.h"

SUITE (ExplicitEdgeWavefunction) {
	double eps = 1e-7;

	TEST (Init) {
		Tree tree = TreeFactory::BalancedTree(10, 2, 3);
		mt19937 gen(34676949);
		TensorTreecd Psi(gen, tree, false);

		/// Transform to symmetric representation
		MatrixTensorTree Chi(Psi, tree, true);

		/// Check re-obtaining wavefunction
		TensorTreecd Psi2 = Chi.BottomUpNormalized(tree);
		auto S = TreeFunctions::DotProduct(Psi, Psi2, tree);
		for (const Node& node : tree) {
			const auto& s = S[node];
				CHECK_CLOSE(0., residual(s * s, identityMatrixcd(s.dim1())), eps);
		}

		/// Check top-down
		const MatrixTreecd rho = TreeFunctions::Contraction(Psi, tree, true);
		const MatrixTreecd T = TreeFunctions::Contraction(Psi2, Psi, S, tree);
		const TensorTreecd& Atilde = Chi.nodes();
		for (const Node& node : tree) {
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				auto x = contraction(Atilde[parent], Atilde[parent], node.childIdx());
				auto v1 = diagonalize(x).second;
				auto v2 = diagonalize(rho[node]).second;
				auto r = residual(v1, v2);
					CHECK_CLOSE(0., r, eps);
			}
		}
		/// @TODO: Add strong test for top-down.
			CHECK_EQUAL(true, IsWorking(Chi, tree, eps));
	}
}

