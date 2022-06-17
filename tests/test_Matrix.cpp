//
// Created by Roman Ellerbrock on 2020-01-17.
//
#include <gtest/gtest.h>
//#include "Core/Matrix.h"
#include "Core/Matrix_Implementation.h"
#include <random>
#include "Util/RandomMatrices.h"


class MatrixFactory : public ::testing::Test {
protected:
    MatrixFactory() {
        CreateMatrices();
    }

    Matrixcd A;
    Matrixcd B;
    Matrix<int> C_;
    Matrixcd realH2x2_;
    Matrixcd realG2x2_;
    Matrixcd imagH2x2_;
    Matrixcd imagG2x2_;
    Vectorcd v_;

    void CreateMatrixA() {
        A = Matrixcd(3, 3);
        for (size_t i = 0; i < A.dim1(); ++i) {
            for (size_t j = 0; j < A.dim2(); ++j) {
                A(j, i) = i + j;
            }
        }
    }

    void CreateMatrixB() {
        B = Matrixcd(3, 3);
        for (size_t i = 0; i < B.dim1(); ++i) {
            for (size_t j = 0; j < B.dim2(); ++j) {
                B(j, i) = i * j;
            }
        }
    }

    void CreateMatrixC() {
        C_ = Matrix<int>(3, 3);
        for (size_t i = 0; i < C_.dim1(); ++i) {
            for (size_t j = 0; j < C_.dim2(); ++j) {
                C_(j, i) = j * 3 + i; /// ordered according to output
            }
        }
    }

    void CreateMatrices() {
        CreateMatrixA();
        CreateMatrixB();
        CreateMatrixC();

        realH2x2_ = Matrixcd(2, 2);
        realH2x2_(0, 0) = 1;
        realH2x2_(1, 0) = 2;
        realH2x2_(0, 1) = 2;
        realH2x2_(1, 1) = 1;

        complex<double> im(0., 1.);
        imagH2x2_ = im * realH2x2_;

        realG2x2_ = Matrixcd(2, 2);
        realG2x2_(0, 0) = 1;
        realG2x2_(1, 0) = 2;
        realG2x2_(0, 1) = 3;
        realG2x2_(1, 1) = 4;

        imagG2x2_ = im * realG2x2_;

        v_ = Vectorcd(2);
        v_(0) = 1;
        v_(1) = 2;
    }
};

Matrixcd Create() {
    Matrixcd tmp(3, 3);
    for (size_t i = 0; i < tmp.dim1(); ++i) {
        for (size_t j = 0; j < tmp.dim2(); ++j) {
            tmp(j, i) = i + j;
        }
    }
    return tmp;
}

constexpr double eps = 1e-7;

TEST_F (MatrixFactory, Matrix_FileIO) {
    /// Test Matrix I/O
    A.write("matrix1.tmp.dat");
    Matrixcd N("matrix1.tmp.dat");
    bool success = A == N;
    ASSERT_EQ(success, true);
}

TEST_F (MatrixFactory, Matrix_Add) {
    auto S = A + B;
    S.write("matrix_add.dat");
    Matrixcd S_read("matrix_add.dat");
    double r = residual(S, S);
    ASSERT_NEAR(r, 0., eps);
}

TEST_F (MatrixFactory, Matrix_Subst) {
    auto D = A - B;
    D.write("matrix_subst.dat");
    Matrixcd D_read("matrix_subst.dat");
    ASSERT_NEAR(residual(D, D_read), 0., eps);
}

TEST_F (MatrixFactory, Matrix_Prod) {
    auto D = A * B;
    D.write("matrix_prod.dat");
    Matrixcd D_read("matrix_prod.dat");
    ASSERT_NEAR(residual(D, D_read), 0., eps);
}


TEST_F (MatrixFactory, Matrix_Diagonalization) {
    auto x = A.cDiag();
    const Matrixcd& Ua = x.first;
    const Vectord& la = x.second;
    Ua.write("matrix_cdiag_trafo.dat");
    la.write("matrix_cdiag_ev.dat");

    Matrixcd U("matrix_cdiag_trafo.dat");
    Vectord lambda("matrix_cdiag_ev.dat");
    auto residual_U = residual(U, x.first);
    auto residual_L = residual(lambda, x.second);
    ASSERT_NEAR(residual_U, 0., eps);
    ASSERT_NEAR(residual_L, 0., eps);
}

TEST_F (MatrixFactory, Matrix_RoF) {
    {
        // Copy asignment operator
        auto Aca = A;
        double r = residual(A, Aca);
        ASSERT_NEAR(r, 0., eps);
    }

    {
        // Copy constructor
        auto Acc(A);
        double r = residual(A, Acc);
        ASSERT_NEAR(r, 0., eps);
    }

    {
        // Move asignment operator
        auto Ama = move(Create());
        double r = residual(A, Ama);
        ASSERT_NEAR(r, 0., eps);
    }

    {
        // Move constructor
        auto Amc(move(Create()));
        ASSERT_NEAR(residual(A, Amc), 0., eps);
    }
}

TEST (Matrix, Matrix_Rebuild) {
    mt19937 gen(1990);
    size_t dim = 10;
    Matrixcd A = RandomMatrices::gue(dim, gen);
    auto x = diagonalize(A);
    Matrixcd B = toMatrix(x);
    auto res = residual(A, B);
    ASSERT_NEAR(0., res, eps);
}

TEST (Matrix, Matrix_BuildInverse) {
    mt19937 gen(1990);
    size_t dim = 10;
    Matrixcd A = RandomMatrices::gue(dim, gen);
    Matrixcd Adag = A.adjoint();
    A = A * Adag;
    auto x = diagonalize(A);

    Matrixcd B = toMatrix(inverse(x, 1e-10));
    Matrixcd I_test = A * B;
    Matrixcd I = identityMatrix<complex<double>>(A.dim2());
    auto res = residual(I, I_test);
    ASSERT_NEAR(0., res, eps);
}

TEST_F (MatrixFactory, add) {
    Matrixcd res = realH2x2_ + realH2x2_;
    ASSERT_NEAR(2., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(2., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, addEqual) {
    Matrixcd res = realH2x2_;
    res += realH2x2_;
    ASSERT_NEAR(2., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(2., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, substract) {
    Matrixcd res = realH2x2_ - realH2x2_;
    ASSERT_NEAR(0., abs(res(0, 0)), eps);
    ASSERT_NEAR(0., abs(res(1, 0)), eps);
    ASSERT_NEAR(0., abs(res(0, 1)), eps);
    ASSERT_NEAR(0., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, substractEqual) {
    Matrixcd res = realH2x2_;
    res -= realH2x2_;
    ASSERT_NEAR(0., abs(res(0, 0)), eps);
    ASSERT_NEAR(0., abs(res(1, 0)), eps);
    ASSERT_NEAR(0., abs(res(0, 1)), eps);
    ASSERT_NEAR(0., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, multiplyScalarleft) {
    Matrixcd res = 2. * realH2x2_;
    ASSERT_NEAR(2., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(2., abs(res(1, 1)), eps);
}


TEST_F (MatrixFactory, multiplyScalarRight) {
    Matrixcd res = realH2x2_ * 2.;
    ASSERT_NEAR(2., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(2., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, multiplyScalarEqual) {
    Matrixcd res = realH2x2_;
    res *= 2;
    ASSERT_NEAR(2., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(2., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, divideScalarEqual) {
    Matrixcd res = realH2x2_;
    res /= 2.;
    ASSERT_NEAR(0.5, abs(res(0, 0)), eps);
    ASSERT_NEAR(1., abs(res(1, 0)), eps);
    ASSERT_NEAR(1., abs(res(0, 1)), eps);
    ASSERT_NEAR(0.5, abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, matrixVector) {
    auto res = realH2x2_ * v_;
    ASSERT_NEAR(5., abs(res(0)), eps);
    ASSERT_NEAR(4., abs(res(1)), eps);
}

TEST_F (MatrixFactory, matrixMatrix) {
    auto res = realH2x2_ * realH2x2_;
    ASSERT_NEAR(5., abs(res(0, 0)), eps);
    ASSERT_NEAR(4., abs(res(1, 0)), eps);
    ASSERT_NEAR(4., abs(res(0, 1)), eps);
    ASSERT_NEAR(5., abs(res(1, 1)), eps);
}

TEST_F (MatrixFactory, equal) {
    auto Acopy = A;
    ASSERT_EQ(true, Acopy == A);
    ASSERT_EQ(false, A == B);
}

TEST_F (MatrixFactory, notEqual) {
    auto Acopy = A;
    ASSERT_EQ(false, Acopy != A);
    ASSERT_EQ(true, A != B);
}

TEST_F (MatrixFactory, frobeniusNorm) {
    double norm = realH2x2_.frobeniusNorm();
    ASSERT_NEAR(10., norm * norm, eps);
}

TEST_F (MatrixFactory, trace) {
    complex<double> trace = realH2x2_.trace();
    ASSERT_NEAR(2., real(trace), eps);
    ASSERT_NEAR(0., imag(trace), eps);
}

TEST_F (MatrixFactory, adjoint) {
    Matrixcd adj = imagH2x2_.adjoint();
    auto res = residual(imagH2x2_, -1. * adj);
    ASSERT_NEAR(0., res, eps);
}

TEST_F (MatrixFactory, transpose) {
    Matrixcd tra = realG2x2_.transpose();
    ASSERT_NEAR(0., abs(tra(0, 0) - 1.), eps);
    ASSERT_NEAR(0., abs(tra(1, 0) - 3.), eps);
    ASSERT_NEAR(0., abs(tra(0, 1) - 2.), eps);
    ASSERT_NEAR(0., abs(tra(1, 1) - 4.), eps);
}

TEST_F (MatrixFactory, zero) {
    Matrixcd zero = realG2x2_;
    zero.zero();
    ASSERT_EQ(0., zero(0, 0));
    ASSERT_EQ(0., zero(1, 0));
    ASSERT_EQ(0., zero(0, 1));
    ASSERT_EQ(0., zero(1, 1));
}

TEST_F (MatrixFactory, identityMatrix) {
    Matrixcd I = identityMatrixcd(2);
    ASSERT_EQ(1., I(0, 0));
    ASSERT_EQ(0., I(1, 0));
    ASSERT_EQ(0., I(0, 1));
    ASSERT_EQ(1., I(1, 1));
}

TEST_F (MatrixFactory, conjugate) {
    Matrixcd I = imagG2x2_ + realG2x2_;
    I = I.conjugate();
    ASSERT_EQ(complex<double>(1., -1.), I(0, 0));
    ASSERT_EQ(complex<double>(2., -2.), I(1, 0));
    ASSERT_EQ(complex<double>(3., -3.), I(0, 1));
    ASSERT_EQ(complex<double>(4., -4.), I(1, 1));
}

TEST_F (MatrixFactory, submatrix) {
    Matrixcd sub = subMatrix(A, 2, 2);
    ASSERT_EQ(0., sub(0, 0));
    ASSERT_EQ(1., sub(1, 0));
    ASSERT_EQ(1., sub(0, 1));
    ASSERT_EQ(2., sub(1, 1));
}

TEST_F (MatrixFactory, residual) {
    ASSERT_NEAR(0., residual(realG2x2_, realG2x2_), eps);
    ASSERT_NEAR(0., residual(imagG2x2_, imagG2x2_), eps);
}

TEST_F (MatrixFactory, row) {
    auto row = C_.row(1);
    ASSERT_EQ(3, row(0));
    ASSERT_EQ(4, row(1));
    ASSERT_EQ(5, row(2));
}

TEST_F (MatrixFactory, col) {
    auto col = C_.col(1);
    ASSERT_EQ(1, col(0));
    ASSERT_EQ(4, col(1));
    ASSERT_EQ(7, col(2));
}

TEST (Matrix, svd) {
    mt19937 gen(1293123);
    size_t dim = 20;
    Matrixcd A = RandomMatrices::gue(dim, gen);
    auto Asvd = svd(A);
    auto B = toMatrix(Asvd);
    ASSERT_NEAR(0., residual(A, B), 1e-8);
}
