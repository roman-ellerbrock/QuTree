//
// Created by Roman Ellerbrock on 2020-01-24.
//
#include "UnitTest++/UnitTest++.h"
#include "LeafOperator.h"
#include "LeafMatrix.h"
#include "MultiLeafOperator.h"
#include "HO_Basis.h"
#include "TreeShape/Tree.h"

SUITE (Operators) {
	class HelperFactory {
	public:
		HelperFactory() {
			Initialize();
		}
		~HelperFactory() = default;
		Matrixcd X;
		LeafMatrixcd x;
		HO_Basis ho;

		void Initialize() {
			// Generate an bit-flip operator and Fmatrix
			X = Matrixcd(2, 2);
			X(0, 1) = 1.;
			X(1, 0) = 1.;
			x = LeafMatrixcd(X);
			ho = HO_Basis(10);
		}
	};

	TEST_FIXTURE (HelperFactory, SPO_HO) {
		HO_Basis ho(10);
		ho.Initialize(1., 0., 0., 1.);
		TensorDim tdim(vector<size_t>({10, 1}));
		Tensorcd A(tdim);
		ho.InitSPF(A);
		auto xA = ho.applyX(A);
		string file("HO_Applyx.dat");
		xA.Write(file);
		Tensorcd B(file);
			CHECK_EQUAL(xA, B);
	}

	TEST_FIXTURE (HelperFactory, SPO_SPOM) {
		// Check that applying a Matrix to a tensor
		// or the corresponding SPOM does the same

		TensorDim tdim(vector<size_t>({2, 1}));
		Tensorcd A(tdim);
		A(0) = 1.;
		Tensorcd XA = multAB(X, A, 0);
		Tensorcd xA(tdim);
		x.Apply(ho, xA, A);
			CHECK_EQUAL(XA, xA);
	}

	TEST_FIXTURE (HelperFactory, MPO_1) {
		MLOcd M(x, 1);
		mt19937 gen(time(nullptr));
		Tree tree(4, 2, 2);
		TensorTreecd Chi(tree);
		Chi.FillRandom(gen, tree, false);
		auto Psi = M.Apply(Chi, tree);

		/// Checking
		string filename("MLO.Apply.dat");
		Psi.Write(filename);
		TensorTreecd Xi(filename);
			CHECK_EQUAL(Xi.size(), Psi.size());
		for (const Node& node : tree) {
				CHECK_EQUAL(Xi[node], Psi[node]);
		}
	}
}