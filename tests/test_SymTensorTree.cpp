//
// Created by Roman Ellerbrock on 1/31/21.
//

#include "UnitTest++/UnitTest++.h"
#include "TreeClasses/SymTensorTree.h"
#include "TreeShape/TreeFactory.h"
#include "Core/Tensor_Extension.h"

SUITE(SymTensorTree) {
	double eps = 1e-7;

	TEST(TensorFlatten) {
		TensorShape shape({2, 3, 4});
		mt19937 gen(19239123);
		Tensorcd A(shape);

		size_t k = 1;
		auto Aflat = toMatrix(A, k);
		auto A2 = toTensor(Aflat, A.shape(), k);
		CHECK_CLOSE(0., Residual(A, A2), eps);
	}

	TEST(solveSLE) {
		TensorShape shape({2, 3, 4});
		mt19937 gen(19239123);
		/// Generate Random complex tensor
		Tensorcd A(shape);
		uniform_real_distribution<double>dist (-1., 1.);
		for (size_t I = 0; I < shape.totalDimension(); ++I) {
			double x = dist(gen);
			double y = dist(gen);
			A(I) = complex<double>(x, y);
//			A(I) = x;
		}
		/// Test for one contraction
		size_t k = 1;
		Matrixcd rho = Contraction(A, A, k);
		Matrixcd B = toMatrix(sqrt(Diagonalize(rho)));
		/// At the moment this solves Anorm * B = A for Anorm.
		/// Should I solve: B * Anorm = A?
		auto Anorm = solveSLE(B, A, k);
		auto s = Contraction(Anorm, Anorm, k);
			CHECK_CLOSE(0., Residual(s, IdentityMatrixcd(s.Dim1())), eps);
//		auto A2 = MatrixTensor(B, Anorm, k);
		auto A2 = TensorMatrix(Anorm, B, k); /// Do I want this ot the version above?
		auto r = Residual(A, A2);
			CHECK_CLOSE(0., r, eps);
	}

	TEST(CreateSym) {
		Tree tree = TreeFactory::BalancedTree(10, 2, 3);
		mt19937 gen(34676949);
		TensorTreecd Psi(gen, tree, false);

		SymTensorTree sPsi(Psi, tree);
		/// Check top-down against QR scheme
		for (const Node& node : tree) {
			/// Top-down normalized
			if (!node.isToplayer()) {
				const Node& parent = node.parent();
				parent.info();
				const auto& w = sPsi.weighted_[parent];
				auto Anorm = normalizedTensor(w, node.childIdx());
				auto s = Anorm.DotProduct(sPsi.down_[node]);
				s.print();
			}
		}
		for (const Node& node : tree) {
			node.info();
			const auto& w = sPsi.weighted_[node];
			auto Anorm = normalizedTensor(w, node.parentIdx());
			auto s = Anorm.DotProduct(sPsi.up_[node]);
			s.print();
		}
		getchar();
	}


}
