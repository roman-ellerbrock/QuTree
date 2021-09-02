//
// Created by Roman Ellerbrock on 5/23/20.
//

#include <UnitTest++/UnitTest++.h>
#include "TreeOperators/TensorOperators/TensorTreeOperator.h"
#include "TreeShape/TreeFactory.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeShape/LeafTypes/SpinGroup.h"
#include "TreeClasses/SpectralDecompositionTree.h"
#include "TreeOperators/TensorOperators/TTOrepresentation.h"
#include "TreeOperators/TensorOperators/TTOcontraction.h"
#include "Util/GateOperators.h"
#include "TreeOperators/TensorOperators/contractSOP.h"
#include "TreeOperators/TensorOperators/contractCircuit.h"

SUITE (TensorTreeOperator) {

	double eps = 1e-7;

	class HelperFactory {
	public:
		HelperFactory() {
			Initialize();
		}

		~HelperFactory() = default;
		mt19937 rng_;
		Tree tree_;
		Tree optree_;
		TensorTreed Psi_;
		TensorTreeOperatord H_;
		LeafMatrixd leafI_;
		LeafMatrixd leafX_;

		void Initialize() {
			rng_ = mt19937(1990);
			tree_ = TreeFactory::balancedTree(12, 2, 2);
			optree_ = TreeFactory::operatorTree(tree_);
			Psi_ = TensorTreed(rng_, tree_);

			Matrixd I = identityMatrixd(2);
			Matrixd X(2, 2);
			X(0, 1) = 1.;
			X(1, 0) = 1.;
			leafI_ = LeafMatrixd(I);
			leafX_ = LeafMatrixd(X);

			H_ = TensorTreeOperatord(optree_);
			H_.occupy(optree_);
			for (const Node& node : optree_) {
				if (node.isBottomlayer()) {
					H_.setLeafOperator(leafI_, 0, node);
					H_.setLeafOperator(leafX_, 1, node);
				}
			}
		}
	};

/*	TEST (OperatorTree) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		Tree optree = TreeFactory::operatorTree(tree);
		for (const Node& node : optree) {
			if (node.isBottomlayer()) {
					CHECK_EQUAL(2, node.shape().order());
					CHECK_EQUAL(4, node.shape().lastDimension());
					CHECK_EQUAL(4, node.shape().lastBefore());
			}
		}
	}*/

	TEST (TTNO) {
		Tree tree = TreeFactory::balancedTree(12, 2, 2);
		Tree optree = TreeFactory::operatorTree(tree);

		mt19937 gen(2348);
		TensorTreeOperatord A(optree);
		A.occupy(optree, gen);
		for (const Node& node : optree) {
			const Tensord& B = A[node];
			for (size_t i = 0; i < node.shape().totalDimension(); ++i) {
				double r = abs(B[i]);
					CHECK_EQUAL((r > 1e-15), true);
			}
		}
	}

	TEST (TTNOrep) {
		SOPd S;
		Tree tree = TreeFactory::balancedTree(32, 2, 3);
		Tree optree = TreeFactory::operatorTree(tree);

		for (size_t l = 0; l < tree.nLeaves(); l++) {
			Matrixd sigma = JordanWigner::sigmaX();
			MLOd M(sigma, l);
			S.push_back(M, 1.);
		}

		mt19937 gen(time(NULL));
		TensorTreeOperatord A(optree, gen);

		orthogonal(A, optree);
		orthonormal(A, optree);

		double err0 = error(A, S, optree);
			CHECK_EQUAL(1, (err0 > 1e-1));
		contractSOP(A, S, 2, optree, nullptr);
		double err = error(A, S, optree);
			CHECK_EQUAL(1, (err < 1e-12));
	}

	TEST (TTNOrep_nonhermitian) {
		SOPd S;
		Tree tree = TreeFactory::balancedTree(6, 2, 3);
		Tree optree = TreeFactory::operatorTree(tree, 3);

		for (size_t l = 0; l < tree.nLeaves(); l++) {
			Matrixd sigma = JordanWigner::sigmaPlus();
			if (l % 2) { sigma = JordanWigner::sigmaMinus(); }
			MLOd M(sigma, l);
			for (size_t i = 0; i < l; ++i) {
				if (l != i) { M.push_back(JordanWigner::sigmaZ(), i); }
			}
			S.push_back(M, 1. / (double) (l + 1.));
		}

		mt19937 gen(time(NULL));
		TensorTreeOperatord A(optree, gen);

		orthogonal(A, optree);
		orthonormal(A, optree);

		double err0 = error(A, S, optree);
			CHECK_EQUAL(1, (err0 > 1e-1));
		contractSOP(A, S, 3, optree, nullptr);
		double err = error(A, S, optree);
			CHECK_EQUAL(1, (err < 1e-12));

		TTOrepresentationd rep(tree, optree);
		TensorTreed Psi(gen, tree);
		rep.calculate(Psi, A, Psi, optree);
			CHECK_CLOSE(0., rep[optree.topNode()][0][0], eps);

		TTOcontractiond con(tree, optree);
		con.calculate(Psi, A, rep, Psi, optree);
	}

	TEST (Circuit_contraction) {
		mt19937 gen(time(NULL));
		Tree tree = TreeFactory::balancedTree(3, 2, 4);
		Tree optree = TreeFactory::operatorTree(tree, 4);

		/// Build using regular scheme
		TensorTreeOperatorcd A(optree, gen);
		{
			SOPcd circ(CNot(0, 1));
			for (size_t i = 1; i < tree.nLeaves() - 1; ++i) {
				circ = CNot(i, i + 1) * circ;
			}
			contractSOP(A, circ, 10, optree, nullptr);
		}

		/// Build from new scheme
		TensorTreeOperatorcd U(optree, gen);
		{
			SOPcd circ2 = (CNot(0, 1));
			contractSOP(U, circ2, 10, optree, nullptr);
		}

		SOPcd cnot = CNot(1, 2);
		U = contractCircuit(cnot, U, 2, optree, nullptr);
		canonicalTransformation(U, optree, true);

		Matrixcd E(1, 1);
		{
			auto S_aa = TreeFunctions::dotProduct(A, A, optree);
			E += S_aa[optree.topNode()];
		}
		{
			auto S_ab = TreeFunctions::dotProduct(A, U, optree);
			E -= 2. * S_ab[optree.topNode()];
		}
		{
			auto S_bb = TreeFunctions::dotProduct(U, U, optree);
			E += S_bb[optree.topNode()];
		}
			CHECK_CLOSE(0., real(E(0, 0)), 1e-10);
	}

	TEST(contractCircuit_CNotChain) {
		mt19937 gen(time(NULL));
		Tree tree = TreeFactory::balancedTree(6, 2, 4);
		Tree optree = TreeFactory::operatorTree(tree, 4);

		/// Build using regular scheme
		TensorTreeOperatorcd A(optree, gen);
		SOPVectorcd circuit;
		SOPcd circ(CNot(0, 1));
		circuit.push_back(CNot(0, 1));
		for (size_t i = 1; i < tree.nLeaves() - 1; ++i) {
			circ = CNot(i, i + 1) * circ;
			circuit.push_back(CNot(i, i + 1));
		}

		contractSOP(A, circ, 10, optree, nullptr);

		/// Perform same operators iteratively
		auto U = contractCircuit(circuit, 5, optree, &cout);
		Matrixcd E(1, 1);
		{
			auto S_aa = TreeFunctions::dotProduct(A, A, optree);
			E += S_aa[optree.topNode()];
		}
		{
			auto S_ab = TreeFunctions::dotProduct(A, U, optree);
			E -= 2. * S_ab[optree.topNode()];
		}
		{
			auto S_bb = TreeFunctions::dotProduct(U, U, optree);
			E += S_bb[optree.topNode()];
		}
			CHECK_CLOSE(0., real(E(0, 0)), 1e-10);
	}

	TEST(QFT) {
		/// @TODO: switch qubit order to even / uneven
		size_t nleaves = 4;
		Tree tree = TreeFactory::balancedTree(nleaves, 2, 4);
		Tree optree = TreeFactory::operatorTree(tree, 4);
		auto qft = QFT(0, nleaves);

		auto A = contractCircuit(qft, 10, optree, &cout);
		canonicalTransformation(A, optree);
		cout << "QFT: =============>\n";
		A.print(optree);
		cout << "<============= :QFT\n";
		auto rho = TreeFunctions::contraction(A, optree, true);
		rho.print(optree);
	}
}

