//
// Created by Roman Ellerbrock on 11/18/21.
//

#include <gtest/gtest.h>
#include <iostream>
#include "Tensor/TensorBLAS2.h"
#include "Tensor/TensorSlow.h"
#include "Util/QMConstants.h"


class TensorBLAS1Factory : public ::testing::Test {
public:
	TensorBLAS1Factory() {
		A = arangecd({2, 3, 4, 2});

		B = Tensorcd(A.shape_);
		for (size_t i = 0; i < B.shape_.totalDimension(); ++i) {
			if (i % 3 == 0) {
				B(i) = 1.;
			}
		}

		C_ = aranged({5, 10});
	}

	Tensorcd A;
	Tensorcd B;
	Tensord C_;
};

constexpr double eps = 1e-7;

TEST_F (TensorBLAS1Factory, nrm2) { /// ||A||_2
	double norm = pow(nrm2(B), 2);
		EXPECT_NEAR(16., norm, eps);

}


TEST_F (TensorBLAS1Factory, axpy) { /// vec add
	Tensorcd D(A.shape_);
	complex<double> alpha = -0.5;
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		D(i) = alpha * A(i) + B(i);
	}

	axpy(A, B, alpha);

	double res = 0.;
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		res += pow(abs(D(i) - B(i)), 2);
	}
	res = sqrt(res);

		EXPECT_NEAR(0., res, eps);
}

TEST_F (TensorBLAS1Factory, add) {
	Tensorcd D(A.shape_);
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		D(i) = A(i) + B(i);
	}

	A += B;

	double res = 0.;
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		res += pow(abs(D(i) - A(i)), 2);
	}
	res = sqrt(res);

		EXPECT_NEAR(0., res, eps);
}

TEST_F (TensorBLAS1Factory, subst) {
	Tensorcd D(A.shape_);
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		D(i) = A(i) - B(i);
	}

	A -= B;

	double res = 0.;
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		res += pow(abs(D(i) - A(i)), 2);
	}
	res = sqrt(res);

		EXPECT_NEAR(0., res, eps);
}

TEST_F (TensorBLAS1Factory, residual) {
	Tensorcd D(A);
	D -= B;

	double res = 0.;
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		res += pow(abs(D(i) - A(i)), 2);
	}
	res = sqrt(res);

	double res2 = abs(residual(D, A));
		EXPECT_NEAR(0., abs(res - res2), eps);
}

TEST_F (TensorBLAS1Factory, plus_op) {
	Tensorcd AaddB = A + B;
	Tensorcd ApB = A;
	ApB += B;
		EXPECT_NEAR(0., abs(residual(AaddB, ApB)), eps);
}


TEST_F (TensorBLAS1Factory, minus_op) {
	Tensorcd AB = A - B;
	Tensorcd AmB = A;
	AmB -= B;
		EXPECT_NEAR(0., abs(residual(AB, AmB)), eps);
}

TEST_F (TensorBLAS1Factory, scaleeq_op) {
	Tensorcd aA = A;
	aA *= 2.;
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) *= 2;
	}
		EXPECT_NEAR(0., abs(residual(aA, A)), eps);
}


TEST_F (TensorBLAS1Factory, scaleleft_op) {
	Tensorcd aA = 2. * A;
	A *= 2.;
		EXPECT_NEAR(0., abs(residual(aA, A)), eps);
}

TEST_F (TensorBLAS1Factory, scaleright_op) {
	Tensorcd Aa = A * 2.;
	A *= 2.;
		EXPECT_NEAR(0., abs(residual(Aa, A)), eps);
}

TEST_F (TensorBLAS1Factory, diveq_op) {
	Tensorcd Aa = A;
	Aa /= 2.;
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) /= 2;
	}
		EXPECT_NEAR(0., abs(residual(Aa, A)), eps);
}

TEST_F (TensorBLAS1Factory, divight_op) {
	Tensorcd Aa = A / 2.;
	A /= 2.;
		EXPECT_NEAR(0., abs(residual(Aa, A)), eps);
}

TEST_F (TensorBLAS1Factory, adjust_inc_dec) {
	size_t leaf = 1;
	size_t dim = A.shape_[leaf];
	size_t inc_dim = dim + 1;

	auto B = A.adjustActiveDim(inc_dim, leaf);
	B = B.adjustActiveDim(dim, leaf);
		EXPECT_NEAR(0., residual(A, B), eps);
}

TEST_F (TensorBLAS1Factory, productElementwise) {
	auto C = productElementwise(A, B);
	Tensorcd D(A);
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		D(i) *= B(i);
	}
		EXPECT_NEAR(0., residual(C, D), eps);
}

TEST_F (TensorBLAS1Factory, conj) {
	complex<double> im(0., 1.);
	A *= im;
	Tensorcd D(A);
	for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
		D(i) *= -1.;
	}
	A = conj(A);
		EXPECT_NEAR(0., residual(A, D), eps);
}

TEST_F (TensorBLAS1Factory, diag) {
	Tensord mat = aranged({3, 3});
	Tensord diag = diagonal(mat);
		EXPECT_EQ(3, diag.shape_[0]);
		EXPECT_EQ(1, diag.shape_.order());
		EXPECT_NEAR(0., diag(0), eps);
		EXPECT_NEAR(4., diag(1), eps);
		EXPECT_NEAR(8., diag(2), eps);
}

TEST_F (TensorBLAS1Factory, diagC) {
	auto diag = diagonal(C_);
	EXPECT_EQ(5, diag.shape_[0]);
	EXPECT_NEAR(0., diag(0), eps);
	EXPECT_NEAR(6., diag(1), eps);
	EXPECT_NEAR(12., diag(2), eps);
	EXPECT_NEAR(18., diag(3), eps);
	EXPECT_NEAR(24., diag(4), eps);
}

TEST_F (TensorBLAS1Factory, off_diag) {
	Tensord mat = aranged({3, 3});
	Tensord diag = diagonal(mat);
	Tensord off = mat;
	off(0,0) = 0.;
	off(1,1) = 0.;
	off(2,2) = 0.;
	Tensord off2(mat.shape_);
	offDiagonal(off2, mat);
	EXPECT_NEAR(0., residual(off, off2), eps);
}

TEST_F (TensorBLAS1Factory, trace) {
		EXPECT_NEAR(0., abs(trace(B) - 1.), eps);
}

TEST_F (TensorBLAS1Factory, transpose_1) {
	Tensorcd B(A.shape_);
	size_t dim1 = A.shape_.lastBefore();
	size_t dim2 = A.shape_.lastDimension();
	transpose(B.data(), A.data(), dim1, dim2, (complex<double>) 0.);

	Tensorcd C({dim2, dim1});
	for (size_t i = 0; i < C.shape_.lastDimension(); ++i) {
		for (size_t j = 0; j < C.shape_.lastBefore(); ++j) {
			C(j, i) = A(i, j);
		}
	}
		EXPECT_NEAR(0., residual(B, C), eps);
}

TEST_F (TensorBLAS1Factory, transpose_AB) {
	Tensorcd B(A.shape_);
	size_t k = 2;
	size_t a = A.shape_.before(k);
	size_t b = A.shape_[k];
	size_t c = A.shape_.after(k);
	transposeAB(B.data(), A.data(), a, b, c);

	Tensorcd C({b, a, c});
	for (size_t aft = 0; aft < c; ++aft) {
		for (size_t act = 0; act < b; ++act) {
			for (size_t bef = 0; bef < a; ++bef) {
				C(act, bef, aft, 1) = A(bef, act, aft, k);
			}
		}
	}
		EXPECT_NEAR(0., residual(B, C), eps);
}


TEST_F (TensorBLAS1Factory, transpose_BC) {
	Tensorcd B(A.shape_);
	size_t k = 2;
	size_t a = A.shape_.before(k);
	size_t b = A.shape_[k];
	size_t c = A.shape_.after(k);
	transposeBC(B.data(), A.data(), a, b, c);

	Tensorcd C({a, c, b});
	for (size_t aft = 0; aft < c; ++aft) {
		for (size_t act = 0; act < b; ++act) {
			for (size_t bef = 0; bef < a; ++bef) {
				C(bef, aft, act, 1) = A(bef, act, aft, k);
			}
		}
	}
		EXPECT_NEAR(0., residual(B, C), eps);
}

TEST (TensorBLAS, shape_transpose) {
	size_t k = 2;
	TensorShape shape({2, 3, 4, 5, 6});

	shape = transpose(shape, k, false);
	TensorShape shape_there({2, 3, 5, 6, 4});
		EXPECT_EQ(true, (shape == shape_there));

	shape = transpose(shape, k, true);
	TensorShape shape_back({2, 3, 4, 5, 6});
		EXPECT_EQ(shape_back, shape);
}

TEST (TensorBLAS, shape_transposeToFront) {
	size_t k = 2;
	TensorShape shape({2, 3, 4, 5, 6});

	shape = transposeToFront(shape, k, false);
	TensorShape shape_there({4, 2, 3, 5, 6});
		EXPECT_EQ(true, (shape == shape_there));

	shape = transposeToFront(shape, k, true);
	TensorShape shape_back({2, 3, 4, 5, 6});
		EXPECT_EQ(shape_back, shape);
}

TEST_F (TensorBLAS1Factory, transpose_tens) {
	Tensorcd B(A.shape_);
	size_t k = 2;
	transpose(B, A, k);

	Tensorcd C(A.shape_);
	size_t a = A.shape_.before(k);
	size_t b = A.shape_[k];
	size_t c = A.shape_.after(k);
	transposeBC(C.data(), A.data(), a, b, c);
		EXPECT_NEAR(0., residual(B, C), eps);
}

TEST_F (TensorBLAS1Factory, transpose_tens_back) {
	Tensorcd B(A.shape_);
	size_t k = 2;
	transpose(B, A, k);

	Tensorcd C(A.shape_);
	transpose(C, B, k, true);
		EXPECT_NEAR(0., residual(C, A), eps);
}

TEST_F (TensorBLAS1Factory, transpose_return) {
	Tensorcd ATref(A.shape_);
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		transpose(ATref, A, k);

		Tensorcd AT = transpose(A, k);
			EXPECT_NEAR(0., residual(AT, ATref), eps);
			EXPECT_EQ(true, (ATref.shape_ == AT.shape_));
	}
}

TEST_F (TensorBLAS1Factory, transpose_mat) {
	size_t last = A.shape_.lastDimension();
	size_t befo = A.shape_.lastBefore();
	Tensorcd R(A.shape_);
	transpose(R.data(), A.data(), befo, last);

	Tensorcd C(A.shape_);
	transpose(C, A);
		EXPECT_NEAR(0., residual(C, R), eps);

	auto D = transpose(A);
		EXPECT_NEAR(0., residual(D, R), eps);

	TensorShape newshape({2, 2, 3, 4});
		EXPECT_EQ(C.shape_, newshape);
		EXPECT_EQ(D.shape_, newshape);
}

TEST_F (TensorBLAS1Factory, adjoint_mat) {
	size_t last = A.shape_.lastDimension();
	size_t befo = A.shape_.lastBefore();
	A *= QM::im;
	Tensorcd R(A.shape_);
	transpose(R.data(), A.data(), befo, last);
	R *= -1.;

	Tensorcd C(A.shape_);
	adjoint(C, A);
		EXPECT_NEAR(0., residual(C, R), eps);

	auto D = adjoint(A);
		EXPECT_NEAR(0., residual(D, R), eps);

	TensorShape newshape({2, 2, 3, 4});
		EXPECT_EQ(C.shape_, newshape);
		EXPECT_EQ(D.shape_, newshape);
}

TEST_F (TensorBLAS1Factory, reshape) {
	TensorShape xshape({4, 3, 2, 2});
	A.reshape(xshape);
		EXPECT_EQ(xshape, A.shape_);
}

TEST_F (TensorBLAS1Factory, resize) {
	TensorShape xshape({5, 4, 3, 3});
	A.resize(xshape);
		EXPECT_EQ(xshape, A.shape_);
	for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
		A(i) = i;
			EXPECT_NEAR(0., abs(A(i) - 1. * i), eps);
	}
}

TEST (TensorBLAS, nRowsCols) {
	TensorShape shape({5, 10});
		EXPECT_EQ(5, nrows(shape));
		EXPECT_EQ(10, ncols(shape));

		EXPECT_EQ(10, nrows(shape, blas::Op::Trans));
		EXPECT_EQ(5, ncols(shape, blas::Op::Trans));

		EXPECT_EQ(10, nrows(shape, blas::Op::ConjTrans));
		EXPECT_EQ(5, ncols(shape, blas::Op::ConjTrans));
}


TEST (TensorBLAS, vectorTensor) {
	Tensord B({3, 3, 3});
	for (size_t i = 0; i < 3; ++i) {
		B({i, i, i}) = i;
	}
	Tensord twoB = 2. * B;
	Tensord vec({3});
	vec(0) = vec(1) = vec(2) = 2.;
	vectorTensor(B, vec, 1);
		EXPECT_NEAR(0, residual(B, twoB), eps);
}

