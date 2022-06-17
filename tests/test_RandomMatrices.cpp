//
// Created by Roman Ellerbrock on 5/5/20.
//
#include <gtest/gtest.h>
#include "Util/RandomMatrices.h"
#include "Core/Matrix_Extension.h"
#include "Util/RandomProjector.h"


using namespace RandomMatrices;

Matrixcd BuildRankReduced(const SpectralDecompositioncd& x, size_t rank) {
    auto ew = x.second;
    for (size_t i = rank; i < ew.dim(); ++i) {
        ew(i) = 0.;
    }
    return toMatrix(SpectralDecompositioncd({x.first, ew}));
}

TEST (RMT, LowRankDiagonalization) {
    mt19937 gen(1239);
    size_t dim = 20;
    size_t rank = 5;
    size_t p = 3;

    /// Build a test matrix
    Matrixcd U(dim, dim);
    Vectord ew(dim);
    U = RandomMatrices::gue(dim, gen);
    for (size_t r = 0; r < rank; ++r) {
        ew(r) = 1.;
    }
    SpectralDecompositioncd x(U, ew);
    Matrixcd A = toMatrix(x);
    double diag = abs(A.trace());
    A /= diag;

    /// Diagonalize
    auto x2 = diagonalize(A);
    auto x3 = RandomMatrices::diagonalizeRandom(A, rank, p, gen);

    Vectord ew_acc = reverse(x2.second);
    Vectord ew_approx = reverse(x3.second);
    auto r = residual(ew_approx, ew_acc);
    ASSERT_NEAR(0., r, 1e-12);
}

TEST (RMT, EVStatisticalProperties) {
    mt19937 gen(1239);
    size_t dim = 50;
    size_t rank = 10;
    size_t p = 5;

    /// Build a test matrix
    Matrixcd U(dim, dim);
    Vectord ew(dim);
    U = RandomMatrices::gue(dim, gen);
    for (size_t r = 0; r < dim; ++r) {
        ew(r) = 1. / (2. * pow(r + 1., 1));
    }
    SpectralDecompositioncd x(U, ew);
    Matrixcd A = toMatrix(x);
    double diag = abs(A.trace());
    A /= diag;

    /// Diagonalize
    auto x2 = diagonalize(A);
    auto x3 = RandomMatrices::diagonalizeRandom(A, rank, p, gen);

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

    auto r = residual(ew_app, ew_acc);
//			ASSERT_NEAR(0., r, 1e-3);
}

TEST (RMT, StatisticalProperties) {
    mt19937 gen(1239);
    size_t dim = 50;
    size_t rank = 15;
    size_t p = 5;

    /// Build a test matrix
    Matrixcd U(dim, dim);
    Vectord ew(dim);
    U = RandomMatrices::gue(dim, gen);
    for (size_t r = 0; r < dim; ++r) {
//			ew(r) = 1./(2.*pow(r+1., 0.5));
        ew(r) = 1. / (2. * log(r + 2.));
    }
    SpectralDecompositioncd y({U, ew});
    Matrixcd A = toMatrix(y);

    auto x = RandomMatrices::diagonalizeRandom(A, rank, p, gen);
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

    ASSERT_NEAR(0., alpha_app, 1e-2);
    ASSERT_NEAR(0., alpha_red, 1e-2);
}

TEST (RMT, randomSVD) {
    mt19937 gen(1239);
    size_t dim = 10;
    size_t rank = 3;
    size_t p = 3;

    /// Build a test matrix
    Vectord ew(dim);
    auto U1 = RandomMatrices::gue(dim, gen);
    auto U2 = RandomMatrices::gue(dim, gen);
    for (size_t r = 0; r < rank; ++r) {
        ew(r) = 1.;
    }
    SVDcd x(U1, U2, ew);
    auto A = toMatrix(x);

    auto svd_acc = svd(A);
    auto C = toMatrix(svd_acc);
    ASSERT_NEAR(0., residual(A, C), 1e-12);

    auto svd_rand = svdRandom(A, rank + p, gen);
    auto B = toMatrix(svd_rand);
    ASSERT_NEAR(0., residual(A, B), 1e-12);
}

TEST (RMT, ProjectRandomGUE) {
    mt19937 gen(1239);
    size_t dim = 50;
    size_t rank = 10;
    Vectord ew(dim);
    auto U1 = RandomMatrices::gue(dim, gen);
    auto U2 = RandomMatrices::gue(dim, gen);
    auto A = U1 * U2;
    A /= 1. * A.dim1();
    /// B = A * P
    auto P = randomProjector(A.dim2(), rank, gen);
    auto B = A * P;

    double avg = 0.;
    for (size_t r = 0.; r < A.dim1(); ++r) {
        auto ui = A.row(r);
        auto vi = B.row(r);
        double contr = (ui.norm() - vi.norm()) / ui.norm();
        avg += contr;
    }
    avg /= 1. * A.dim1();
    ASSERT_NEAR(0., avg, 1e-1);
}

TEST (RMT, Krylov) {
    size_t dim = 50;
    size_t rank = 20;
//		uniform_real_distribution<double> dist(0., 1.);
    normal_distribution<double> dist;
    mt19937 gen(2343854);

    /// Build a test matrix
    Vectord ew(dim);
    auto U1 = RandomMatrices::gue(dim, gen);
    auto U2 = RandomMatrices::gue(dim, gen);
    for (size_t r = 0; r < dim; ++r) {
//			ew(r) = 1. / (r + 1.);
        ew(r) = dist(gen);
    }
    SpectralDecompositioncd x(U1, ew);
    auto A = toMatrix(x);
    x = diagonalize(A);
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
    A /= A.frobeniusNorm();

    /// Build a start vector
    auto r = gaussVector(dim, gen);
    auto r2 = gaussVector(dim, gen);

    auto space = RandomMatrices::buildKrylovSpace(r, A, rank);
    for (size_t i = 0; i < space.size(); ++i) {
        for (size_t j = 0; j < space.size(); ++j) {
            double dot = abs(space[i] * space[j]);
            if (i == j) {
                ASSERT_NEAR(1., dot, 1e-9);
            } else {
                ASSERT_NEAR(0., dot, 1e-9);
            }
        }
    }

    auto P = toMatrix(space);
//		auto P = RandomProjector(dim, rank, gen);
//		auto P = RandomGauss(dim, rank, gen);
    auto PP = P * P.adjoint();

    /// CHECK moment
    auto y = r;
    for (size_t i = 0; i < space.size(); ++i) {
        normalize(y);
        auto PPy = y;
        PPy = PP * PPy;
        normalize(PPy);
        double res = residual(y, PPy);
        ASSERT_NEAR(0., res, 1e-8);
        y = A * y;
    }

    if (false) {
        auto PPA = PP * A;
        PPA /= PPA.frobeniusNorm();

        auto Arr_x = reduceRank(x, rank);
        auto Arr = toMatrix(Arr_x);
        Arr /= PPA.frobeniusNorm();

        cout << "entropy(A) = " << entropy(A) << endl;
        cout << "entropy(PP * A) = " << entropy(PPA) << endl;
        cout << "entropy(Arr) = " << entropy(Arr) << endl;
        cout << "cross-entropy(PP* A,  A) = " << crossEntropy(PPA, A) << endl;
        cout << "cross-entropy(Arr,  A) = " << crossEntropy(Arr, A) << endl;
        cout << "cross-entropy-diff(PP* A,  A) = " << crossEntropyDifference(PPA, A) << endl;
    }
}

TEST (RMT, DiagRandom) {

    /*
    size_t dim = 200;
    size_t rank = 25;
//		uniform_real_distribution<double> dist(0., 1.);
    normal_distribution<double> dist;
    mt19937 gen(2343854);

    /// Build a test matrix
    Vectord ew(dim);
    auto U1 = RandomMatrices::GUE(dim, gen);
    for (size_t r = 0; r < dim; ++r) {
        ew(r) = 1. / (r + 1.);
//			ew(r) = dist(gen);
    }
    SpectralDecompositioncd preA(U1, ew);
    auto A = toMatrix(preA);
    auto x = Diagonalize(A);
    cout << "eigenvalues:\n";
    for (int i = x.second.Dim() - 1; i >= 0; --i) {
        cout << x.second(i) << " ";
    }
    cout << endl;

    cout << "randomDiag eigenvalues:\n";
    vector<Vectord> ews;
    for (size_t qq = 10; qq >= 1; --qq) {
        auto y = Random::DiagonalizeRandom<complex<double>, Matrixcd>
            (A, rank, qq, gen);
        cout << "power: " << qq << endl;
        ews.push_back(y.second);
        for (int i = rank - 1; i >= 0; --i) {
            size_t shift = x.second.Dim() - y.second.Dim();
            cout << (x.second(i + shift) - y.second(i)) / x.second(i + shift) << " ";
        }
        cout << endl;
    }

    ofstream os("randomDiag."+to_string(rank)+".dat");
    os << "# " << dim << "x" << dim << " Matrix with rank = " << rank << "\n";
    size_t q = ews.size();
    for (const auto& e : ews) {
        os << q << " ";
        q--;
        for (size_t i = 0; i < rank; ++i) {
            size_t last = rank - 1;
            size_t lastfull = dim - 1;
            double z = x.second(lastfull - i);
            double zz = abs(e(last - i) - z) / z;
            os << " " << zz;
        }
        os << endl;
    }
    */
}

