//
// Created by Roman Ellerbrock on 5/5/20.
//
#include "UnitTest++/UnitTest++.h"
#include "Util/RandomMatrices.h"
#include "Core/Matrix_Extension.h"


SUITE (RMT) {
	using namespace RandomMatrices;

	Matrixcd BuildRankReduced(const SpectralDecompositioncd& x, size_t rank) {
		auto ew = x.second;
		for (size_t i = rank; i < ew.Dim(); ++i) {
			ew(i) = 0.;
		}
		return toMatrix(SpectralDecompositioncd({x.first, ew}));
	}

	TEST (LowRankDiagonalization) {
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
		Matrixcd A = toMatrix(x);
		double diag = abs(A.Trace());
		A /= diag;

		/// Diagonalize
		auto x2 = Diagonalize(A);
		auto x3 = RandomMatrices::DiagonalizeRandom(A, rank, p, gen);

		Vectord ew_acc = reverse(x2.second);
		Vectord ew_approx = reverse(x3.second);
		auto r = Residual(ew_approx, ew_acc);
			CHECK_CLOSE(0., r, 1e-12);
	}

	TEST (EVStatisticalProperties) {
		mt19937 gen(1239);
		size_t dim = 100;
		size_t rank = 10;
		size_t p = 5;

		/// Build a test matrix
		Matrixcd U(dim, dim);
		Vectord ew(dim);
		U = RandomMatrices::GUE(dim, gen);
		for (size_t r = 0; r < dim; ++r) {
			ew(r) = 1. / (2. * pow(r + 1., 1));
		}
		SpectralDecompositioncd x(U, ew);
		Matrixcd A = toMatrix(x);
		double diag = abs(A.Trace());
		A /= diag;

		/// Diagonalize
		auto x2 = Diagonalize(A);
		auto x3 = RandomMatrices::DiagonalizeRandom(A, rank, p, gen);

		Vectord ew_acc = reverse(x2.second);
		Vectord ew_app = reverse(x3.second);
		Vectord ew_svd(rank + p);
		for (size_t k = 0; k < rank + p; ++k) {
			ew_svd(k) = ew_acc(k);
		}

		/// Check moments
		auto mom_acc = sum(ew_acc);
		auto mom_app = sum(ew_app);
		auto mom_svd = sum(ew_svd);
		ew_app /= mom_app;
		ew_svd /= mom_svd;

		double S_acc = entropy(ew_acc);
		double S_app = entropy(ew_app);
		double S_svd = entropy(ew_svd);

		auto r = Residual(ew_app, ew_acc);
//			CHECK_CLOSE(0., r, 1e-3);
	}

	TEST (StatisticalProperties) {
		mt19937 gen(1239);
		size_t dim = 200;
		size_t rank = 15;
		size_t p = 5;

		/// Build a test matrix
		Matrixcd U(dim, dim);
		Vectord ew(dim);
		U = RandomMatrices::GUE(dim, gen);
		for (size_t r = 0; r < dim; ++r) {
//			ew(r) = 1./(2.*pow(r+1., 0.5));
			ew(r) = 1. / (2. * log(r + 2.));
		}
		SpectralDecompositioncd y({U, ew});
		Matrixcd A = toMatrix(y);

		auto x = RandomMatrices::DiagonalizeRandom(A, rank, p, gen);
		auto Aprime = toMatrix(x);

		auto Ared = BuildRankReduced(y, rank + p);

		auto p_acc = probabilitiyDist(A);
		auto p_app = probabilitiyDist(Aprime);
		auto p_red = probabilitiyDist(Ared);

		double S_acc = entropy(p_acc);
		double S_app = entropy(p_app);
		double S_red = entropy(p_red);
		double H_app = crossEntropy(p_app, p_acc);
		double H_red = crossEntropy(p_red, p_acc);
		double alpha_app = (H_app - S_acc) / S_acc;
		double alpha_red = (H_red - S_acc) / S_acc;

		if (false) {
			cout << "S(p_accurate) = " << S_acc << endl;
			cout << "S(p_random) = " << S_red << endl;
			cout << "S(p_rankr) = " << S_red << endl;
			cout << "H(p_app, p_acc) = " << H_app << endl;
			cout << "H(p_red, p_acc) = " << H_red << endl;
			cout << "DeltaH(p_app, p_acc) = " << H_app - S_acc << endl;
			cout << "DeltaH(p_red, p_acc) = " << H_red - S_acc << endl;
			cout << "alpha_app = " << alpha_app << endl;
			cout << "alpha_red = " << alpha_red << endl;
			cout << "probability distributions:\n";
			p_acc.print();
			p_app.print();
			p_red.print();
		}

			CHECK_CLOSE(0., alpha_app, 1e-2);
			CHECK_CLOSE(0., alpha_red, 1e-2);
	}

	TEST (randomSVD) {
		mt19937 gen(1239);
		size_t dim = 10;
		size_t rank = 3;
		size_t p = 3;

		/// Build a test matrix
		Vectord ew(dim);
		auto U1 = RandomMatrices::GUE(dim, gen);
		auto U2 = RandomMatrices::GUE(dim, gen);
		for (size_t r = 0; r < rank; ++r) {
			ew(r) = 1.;
		}
		SVDcd x(U1, U2, ew);
		auto A = toMatrix(x);

		auto svd_acc = svd(A);
		auto C = toMatrix(svd_acc);
			CHECK_CLOSE(0., Residual(A, C), 1e-12);

		auto svd_rand = svdRandom(A, rank + p, gen);
		auto B = toMatrix(svd_rand);
			CHECK_CLOSE(0., Residual(A, B), 1e-12);
	}

	TEST (ProjectRandomGUE) {
		mt19937 gen(1239);
		size_t dim = 200;
		size_t rank = 100;
		Vectord ew(dim);
		auto U1 = RandomMatrices::GUE(dim, gen);
		auto U2 = RandomMatrices::GUE(dim, gen);
		auto A = U1 * U2;
		A /= 1. * A.Dim1();
		/// B = A * P
		auto P = RandomProjector(A.Dim2(), rank, gen);
		auto B = A * P;

		double avg = 0.;
		for (size_t r = 0.; r < A.Dim1(); ++r) {
			auto ui = A.row(r);
			auto vi = B.row(r);
			double contr = (ui.Norm() - vi.Norm()) / ui.Norm();
			avg += contr;
		}
		avg /= 1. * A.Dim1();
			CHECK_CLOSE(0., avg, 1e-2);
	}

	TEST (Krylov) {
		size_t dim = 128;
		size_t rank = 20;
//		uniform_real_distribution<double> dist(0., 1.);
		normal_distribution<double> dist;
		mt19937 gen(2343854);

		/// Build a test matrix
		Vectord ew(dim);
		auto U1 = RandomMatrices::GUE(dim, gen);
		auto U2 = RandomMatrices::GUE(dim, gen);
		for (size_t r = 0; r < dim; ++r) {
//			ew(r) = 1. / (r + 1.);
			ew(r) = dist(gen);
		}
		SpectralDecompositioncd x(U1, ew);
		auto A = toMatrix(x);
		x = Diagonalize(A);
		ew = x.second;
		double shift = ew(0);
		for (size_t i = 0; i < dim; ++i) {
//			ew(i) -= shift;
		}
		for (size_t i = 0; i < dim; ++i) {
//			ew(i) /= ew(dim - 1);
		}
		x.second = ew;
		A = toMatrix(x);
		A /= A.FrobeniusNorm();

		/// Build a start vector
		auto r = GaussVector(dim, gen);
		auto r2 = GaussVector(dim, gen);

		auto space = RandomMatrices::BuildKrylovSpace(r, A, rank);
		for (size_t i = 0; i < space.size(); ++i) {
			for (size_t j = 0; j < space.size(); ++j) {
				double dot = abs(space[i] * space[j]);
				if (i == j) {
						CHECK_CLOSE(1., dot, 1e-9);
				} else {
						CHECK_CLOSE(0., dot, 1e-9);
				}
			}
		}

		auto P = toMatrix(space);
//		auto P = RandomProjector(dim, rank, gen);
//		auto P = RandomGauss(dim, rank, gen);
		auto PP = P * P.Adjoint();

		/// CHECK moment
		auto y = r;
		for (size_t i = 0; i < space.size(); ++i) {
			normalize(y);
			auto PPy = y;
			PPy = PP * PPy;
			normalize(PPy);
			double res = Residual(y, PPy);
				CHECK_CLOSE(0., res, 1e-8);
			cout << "Residual: " << res << endl;
			y = A * y;
		}

		if (true) {
			auto PPA = PP * A;
			PPA /= PPA.FrobeniusNorm();

			auto Arr_x = reduceRank(x, rank);
			auto Arr = toMatrix(Arr_x);
			Arr /= PPA.FrobeniusNorm();

			cout << "entropy(A) = " << entropy(A) << endl;
			cout << "entropy(PP * A) = " << entropy(PPA) << endl;
			cout << "entropy(Arr) = " << entropy(Arr) << endl;
			cout << "cross-entropy(PP* A,  A) = " << crossEntropy(PPA, A) << endl;
			cout << "cross-entropy(Arr,  A) = " << crossEntropy(Arr, A) << endl;
			cout << "cross-entropy-diff(PP* A,  A) = " << crossEntropyDifference(PPA, A) << endl;
		}
	}
}

