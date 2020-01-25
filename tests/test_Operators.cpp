//
// Created by Roman Ellerbrock on 2020-01-24.
//
#include "UnitTest++/UnitTest++.h"
#include "SingleParticleOperator.h"
#include "SingleParticleOperatorFunction.h"
#include "SingleParticleOperatorMatrix.h"
#include "MultiParticleOperator.h"
#include "SumOfProductsOperator.h"
#include "HO_Basis.h"
#include "TensorTree/TensorTreeBasis.h"

SUITE(Operators) {

	TEST(SPO_HO) {
		HO_Basis ho(10);
		ho.Initialize(1., 0., 0., 1.);
		TensorDim tdim({10}, 1);
		Tensorcd A(tdim);
		ho.InitSPF(A);
		auto xA = ho.applyX(A);
		string file("HO_Applyx.dat");
//		xA.Write(file);
		Tensorcd B(file);
		CHECK_EQUAL(xA, B);
	}

	TEST(SPO_SPOM) {
		// Check that applying a Matrix to a tensor
		// or the corresponding SPOM does the same
		// We need a dummy-basis
		HO_Basis ho(10);
		// Generate an bit-flip operator
		FactorMatrixcd M(2, 1);
		M(0, 1) = 1.;
		M(1, 0) = 1.;
		SPOMcd x(M);

		TensorDim tdim({3, 2, 4}, 1);
		Tensorcd A(tdim);
		A(0) = 1.;
		auto MA = M * A;
		Tensorcd xA(tdim);
		x.Apply(ho, xA, A);
		CHECK_EQUAL(MA, xA);
	}

	TEST(MPO_1) {
		TTBasis basis(2, 2, 2);
		TensorTreecd Psi(basis);
		mt19937 gen(time(nullptr));
		Psi.Generate(basis, gen);
	}

	TEST(SOP_1) {

	}

}