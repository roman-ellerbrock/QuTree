//#include "Core/Tensor.h"
#include "Core/Tensor_Implementation.h"
#include <iostream>
#include <gtest/gtest.h>
#include "Util/QMConstants.h"
#include "Core/Tensor_Extension.h"

using namespace std;

class TensorFactory : public ::testing::Test {
protected:
    TensorFactory() {
        CreateTensors();
    }

    Tensorcd A;
    Tensorcd B;
    Tensorcd C_;
    Tensorcd C2_;
    TensorShape shape_c_;

    void CreateTensorA() {
        TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
        A = Tensorcd(tdim);
        for (size_t i = 0; i < tdim.totalDimension(); ++i) {
            A(i) = i;
        }
    }

    void CreateTensorB() {
        TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
        B = Tensorcd(tdim);
        for (size_t i = 0; i < tdim.totalDimension(); ++i) {
            B(i) = i % 3;
        }
    }

    void CreateTensorC() {
        shape_c_ = TensorShape({2, 2, 2});
        // C
        C_ = Tensorcd(shape_c_);
        for (size_t bef = 0; bef < shape_c_.before(1); ++bef) {
            for (size_t act = 0; act < shape_c_[1]; ++act) {
                for (size_t aft = 0; aft < shape_c_.after(1); ++aft) {
                    C_(bef, act, aft, 1) = (double) act * QM::im;
                }
            }
        }
        // C2
        C2_ = Tensorcd(shape_c_);
        for (size_t bef = 0; bef < shape_c_.before(1); ++bef) {
            for (size_t act = 0; act < shape_c_[1]; ++act) {
                for (size_t aft = 0; aft < shape_c_.after(1); ++aft) {
                    C2_(bef, act, aft, 1) = (double) aft * QM::im;
                }
            }
        }
    }

    void CreateTensors() {
        CreateTensorA();
        CreateTensorB();
        CreateTensorC();
    }
};

Tensorcd NewTensor() {
    TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
    Tensorcd tmp(tdim);
    for (size_t i = 0; i < tdim.totalDimension(); ++i) {
        tmp(i) = i;
    }
    return tmp;
}

constexpr double eps = 1e-7;

TEST (Tensor, TensorDim_FileIO) {
    /// Create a TensorDim, write to file, read in again
    TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
    tdim.write("tdim.dat");
    TensorShape odim("tdim.dat");
    bool same = odim == tdim;
    ASSERT_EQ(same, true);
}

TEST (Tensor, TensorDim_Getters) {
    /// Check Getters and initialization
    bool success = true;
    TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
    ASSERT_EQ(3 * 4 * 5 * 2, tdim.totalDimension());
    ASSERT_EQ(3 * 4 * 5, tdim.lastBefore());
    ASSERT_EQ(2, tdim.lastDimension());
}

TEST (Tensor, Tensor_Constructor) {
    TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
    Tensorcd A(tdim);
    Tensorcd B(tdim);
    auto same = A - B;
    auto s = same.dotProduct(same);
    auto delta = s.frobeniusNorm();
    ASSERT_NEAR(delta, 0., eps);
}

TEST_F (TensorFactory, Tensor_FileIO) {
    /// Test Tensor I/O
    A.write("tensor1.dat");
    Tensorcd B("tensor1.dat");
    Tensorcd C = A - B;
    Matrixcd s = C.dotProduct(C);
    double residual = abs(s.trace());
    ASSERT_NEAR(residual, 0., eps);
}

TEST_F (TensorFactory, Tensor_Product) {
    Matrixcd x = contraction(A, B, 0);
    x.write("Tensor_Product.dat");
    Matrixcd s("Tensor_Product.dat");
    auto r = residual(s, x);
    ASSERT_NEAR(0., r, eps);
}

TEST_F (TensorFactory, Tensor_Matrix_Product) {
    Matrixcd x = contraction(A, B, 1);
    x.write("Tensor_Product_0.dat");
    Matrixcd s("Tensor_Product_0.dat");
    auto r = residual(x, s);
    ASSERT_NEAR(0., r, eps);
}

TEST_F(TensorFactory, Tensor_RoF) {
    {
        // Copy asignment operator
        auto Aca = A;
        ASSERT_NEAR(0., residual(A, Aca), eps);
    }

    {
        // Copy constructor
        auto Acc(A);
        ASSERT_NEAR(0., residual(A, Acc), eps);
    }

    {
        // Move asignment operator
        auto Ama = move(NewTensor());
        ASSERT_NEAR(0., residual(A, Ama), eps);
    }

    {
        // Move constructor
        auto Amc(NewTensor());
        ASSERT_NEAR(0., residual(A, Amc), eps);
    }
}

TEST_F (TensorFactory, AdjustDimension_inc) {
    gramSchmidt(A);
    size_t leaf = 1;
    size_t dim = A.shape()[leaf];
    size_t inc_dim = dim + 1;

    auto C = A.adjustActiveDim(inc_dim, leaf);
    A = A.adjustActiveDim(inc_dim, leaf);
    auto s = contraction(C, A, C.shape().lastIdx());
    auto res = residual(s, identityMatrix<complex<double>>(2));
    ASSERT_NEAR(0., res, eps);
}

TEST_F (TensorFactory, AdjustDimension_inc_dec) {
    gramSchmidt(A);
    size_t leaf = 1;
    size_t dim = A.shape()[leaf];
    size_t inc_dim = dim + 1;

    auto C = A.adjustActiveDim(inc_dim, leaf);
    C = C.adjustActiveDim(dim, leaf);
    auto s = contraction(C, A, C.shape().lastIdx());
    auto res = residual(s, identityMatrix<complex<double>>(2));
    ASSERT_NEAR(0., res, eps);
}

TEST_F (TensorFactory, HoleProduct) {
    Matrixcd s = contraction(C_, C_, 1);
    ASSERT_EQ(shape_c_[1], s.dim1());
    ASSERT_EQ(shape_c_[1], s.dim2());
    double dim = shape_c_.before(1) * shape_c_.after(1);
    for (size_t i = 0; i < shape_c_[1]; ++i) {
        for (size_t j = 0; j < shape_c_[1]; ++j) {
            ASSERT_NEAR(dim * (double) i * (double) j, abs(s(i, j)), eps);
        }
    }
}

TEST_F (TensorFactory, DotProduct) {
    Matrixcd s = C2_.dotProduct(C2_);
    ASSERT_EQ(shape_c_[2], s.dim1());
    ASSERT_EQ(shape_c_[2], s.dim2());
    for (size_t i = 0; i < shape_c_[2]; ++i) {
        for (size_t j = 0; j < shape_c_[2]; ++j) {
            ASSERT_NEAR(shape_c_.before(2) * (double) i * (double) j, abs(s(i, j)), eps);
        }
    }
}

TEST (Tensor, DirectSum) {
    TensorShape Ashape({2, 2});
    TensorShape Bshape({3, 3});
    Tensord A(Ashape);
    for (size_t I = 0; I < Ashape.totalDimension(); ++I) {
        A(I) = I * I + 1;
    }
    Tensord B(Bshape);
    for (size_t I = 0; I < Bshape.totalDimension(); ++I) {
        B(I) = I + 1;
    }
    Tensord C = Tensor_Extension::directSum(A, B, true, true);
    for (size_t I = 0; I < Ashape.totalDimension(); ++I) {
        auto Ibreak = indexMapping(I, Ashape);
        size_t L = indexMapping(Ibreak, C.shape());
        ASSERT_NEAR(A(I), C(L), eps);
    }

    for (size_t I = 0; I < Bshape.totalDimension(); ++I) {
        auto Ibreak = indexMapping(I, Bshape);
        for (size_t k = 0; k < Ashape.order(); ++k) {
            Ibreak[k] += Ashape[k];
        }
        size_t L = indexMapping(Ibreak, C.shape());
        ASSERT_NEAR(B(I), C(L), eps);
    }
}

TEST(Tensor, QRTensor) {
    TensorShape shape({2, 3, 4, 5});
    mt19937 gen(1239);
    Tensorcd A(shape);
    Tensor_Extension::generate(A, gen);

    /// Test standard QR
    Tensorcd Q = qr(A);
    auto S = Q.dotProduct(Q);
    ASSERT_NEAR(0., residual(S, identityMatrixcd(S.dim1())), eps);

    /// Test QR for other than last mode
    for (size_t i = 0; i < shape.order(); ++i) {
        Tensorcd Q2 = qr(A, i);
        auto S1 = contraction(Q2, Q2, i);
        ASSERT_NEAR(0., residual(S1, identityMatrixcd(S1.dim1())), eps);
    }
}

TEST(Tensor, ones) {
	TensorShape shape({3});
	Tensord a = ones<double>(shape);
	ASSERT_NEAR(1., a(0), 1e-12);
	ASSERT_NEAR(1., a(1), 1e-12);
	ASSERT_NEAR(1., a(2), 1e-12);
}


