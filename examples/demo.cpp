//
// Created by Roman Ellerbrock on 6/15/20.
//
#include "TreeClasses/TensorTree.h"
#include "TreeClasses/MatrixTreeFunctions.h"
#include "TreeClasses/SparseMatrixTreeFunctions.h"
#include "TreeShape/TreeFactory.h"
#include "TreeShape/Tree.h"

namespace Demo {

	void TensorAlgebra() {

		/// Dimension-object for a Tensor
		TensorShape shape({3, 4, 5});
		/// Create tensor with complex<double> elements (cd)
		Tensorcd A(shape);
		Tensorcd B = A;
		Tensorcd C = A - B;
		Matrixcd Overlap = A.dotProduct(B);
		/// Tensor-contraction, hole in index 0
		Matrixcd rho = contraction(A, B, 0);


		B(0) = 1.;
	}

	void TreeAlgebra() {
		size_t number_leaves = 7;
		size_t dim_leaves = 2;
		size_t dim_nodes = 2;
		mt19937 RNG(012);


		{
			Tree tree = TreeFactory::BalancedTree(251, 2, 2);
			TensorTreecd Psi(RNG, tree);
			MatrixTreecd Overlap = TreeFunctions::dotProduct(Psi, Psi, tree);
			MatrixTreecd RDMs = TreeFunctions::contraction(Psi, tree, true);
		}


		{

			/// Order 5
			TensorShape shape({2, 2, 2, 2, 2});

		}
		{

			/// dim 200
			TensorShape shape({200, 200, 200});

									 /// Number Leaves 1024
			Tree tree = TreeFactory::BalancedTree(1024, 2, 2);

		}
	}
}