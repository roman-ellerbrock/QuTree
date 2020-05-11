//
// Created by Roman Ellerbrock on 5/5/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Util/RandomMatrices.h"


SUITE (RMT) {

	TEST (testQ) {
		mt19937 gen(1239);
		size_t dim = 20;
		size_t rank = 5;
		size_t p = 3;

		/// Build a test matrix
		Matrixcd U(dim, dim);
		Vectord ew(dim);
		U = RandomMatrices::GUE(dim, gen);
		for (size_t r = 0; r < rank; ++r) {
			ew(r) = 1.;
		}
		SpectralDecompositioncd x(U, ew);
		Matrixcd A = BuildMatrix(x);

		/// Diagonalize
		auto x2 = Diagonalize(A);
		auto x3 = RandomMatrices::DiagonalizeRandom(A, rank, p, gen);

		size_t last3 = rank + p - 1;
		size_t last2 = A.Dim1() - 1;

		Vectord ew_acc = reverse(x2.second);
		Vectord ew_approx = reverse(x3.second);
		auto r = Residual(ew_approx, ew_acc);
		CHECK_CLOSE(0., r, 1e-12);
	}
}