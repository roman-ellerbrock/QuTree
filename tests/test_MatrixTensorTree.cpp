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

		/// Check bottom-up
		{
			const TensorTreecd& Atilde = Chi.nodes();
			auto rho = TreeFunctions::Contraction(Psi, tree, true);
			for (const Node& node : tree) {
				if (!node.isToplayer()) {
					auto x = Contraction(Atilde[node], Atilde[node], node.nChildren());
					auto r = Residual(x, rho[node]);
						CHECK_CLOSE(0., r, eps);
				}
			}

			/// Check top-down
			for (const Node& node : tree) {
				if (!node.isBottomlayer()) {
					for (size_t k = 0; k < node.nChildren(); ++k) {
						const Node& child = node.child(k);
						auto x = Contraction(Atilde[node], Atilde[node], k);
						auto r = Residual(x, rho[child]);
							CHECK_CLOSE(0., r, eps);
					}
				}
			}
		}

			CHECK_EQUAL(true, IsWorking(Chi, tree, eps));
	}


}

