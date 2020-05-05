//
// Created by Roman Ellerbrock on 5/5/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Util/RandomMatrices.h"


SUITE (RMT) {

	TEST (testQ) {
		mt19937 gen(1239);
		size_t dim = 15;
		Matrixcd U(dim, dim);
		Vectord ew(dim);
		U = RandomMatrices::GUE(dim, gen);

		size_t rank = 5;
		for (size_t r = 0; r < rank; ++r) {
			ew(r) = 1.;
		}

		SpectralDecompositioncd x(U, ew);
		Matrixcd A = BuildMatrix(x);

		auto x2 = Diagonalize(A);
		cout << "ew:\n";
		x2.second.print();

		auto Aprime = RandomMatrices::RandomProjection(A, rank + 3, gen);
		auto x3 = Diagonalize(Aprime);
		cout << "approx ew:\n";
		x3.second.print();
	}
}