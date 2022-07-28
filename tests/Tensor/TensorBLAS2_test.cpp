//
// Created by Roman Ellerbrock on 11/7/21.
//

#include <gtest/gtest.h>
#include <iostream>
#include "Tensor/TensorBLAS2.h"
#include "Tensor/TensorSlow.h"
#include "Util/QMConstants.h"


using namespace std;

/**
 * Rationale:
 * The functions matrix-tensor and contraction for tensors
 * are tested for consistency by comparing the results for
 * a straightforward implementation with the BLAS implementation
 */

static double eps = 1e-7;

class TensorFactory2 : public ::testing::Test {
public:
	TensorFactory2() {
		A = arangecd({2, 3, 4, 2});
		CreateTensorB();
	}

	Tensorcd A;
	Tensorcd B;

	Tensorcd constructMatrix(const Tensorcd& A, size_t k) {
		size_t dim = A.shape_[k];
		Tensorcd mat({dim, dim});
		return mat;
	}

	void fill(Tensorcd& mat) {
		for (size_t i = 0; i < mat.shape_[1]; ++i) {
			for (size_t j = 0; j < mat.shape_[0]; ++j) {
				mat(j, i) = i * j + i;
			}
		}
	}

	Tensorcd createMatrix(const Tensorcd& A, size_t k) {
		auto mat = constructMatrix(A, k);
		for (size_t i = 0; i < mat.shape_[1]; ++i) {
			for (size_t j = 0; j < mat.shape_[0]; ++j) {
				mat(j, i) = i * j + i;
			}
		}
		return mat;
	}

	void CreateTensorB() {
		TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
		B = Tensorcd(tdim);
		for (size_t i = 0; i < tdim.totalDimension(); ++i) {
			B(i) = i % 3;
		}
	}
};

TEST_F (TensorFactory2, gemm) {
	auto op_as = {blas::Op::NoTrans, blas::Op::ConjTrans};
	auto op_bs = {blas::Op::NoTrans, blas::Op::ConjTrans};
	complex<double> alpha = 1.5;
	complex<double> beta = 0.5;

	for (auto op_a : op_as) {
		for (auto op_b : op_bs) {
			Tensorcd a({5, 7});
			fill(a);
			if (op_a == blas::Op::ConjTrans) {
				a = adjoint(a);
			}
			Tensorcd b({7, 9});
			fill(b);
			if (op_b == blas::Op::ConjTrans) {
				b = adjoint(b);
			}

			Tensorcd c({5, 9});
			auto c2 = c;
			gemm(c, a, b, alpha, beta, op_a, op_b);
			gemmRef(c2, a, b, alpha, beta, op_a, op_b);
				EXPECT_NEAR(0., residual(c, c2), eps);
		}
	}
}

TEST_F (TensorFactory2, gemm_return) {
	auto op_as = {blas::Op::NoTrans, blas::Op::ConjTrans};
	auto op_bs = {blas::Op::NoTrans, blas::Op::ConjTrans};
	complex<double> alpha = 1.5;

	for (auto op_a : op_as) {
		for (auto op_b : op_bs) {
			Tensorcd a({5, 7});
			fill(a);
			if (op_a == blas::Op::ConjTrans) {
				a = adjoint(a);
			}
			Tensorcd b({7, 9});
			fill(b);
			if (op_b == blas::Op::ConjTrans) {
				b = adjoint(b);
			}

			auto c = gemm(a, b, alpha, op_a, op_b);
			auto c2 = gemmRef(a, b, alpha, op_a, op_b);
				EXPECT_NEAR(0., residual(c, c2), eps);
		}
	}
}

TEST_F(TensorFactory2, OperatorProduct) {
	Tensorcd a({5, 7});
	fill(a);
	Tensorcd b({7, 9});
	fill(b);

	Tensorcd c_ref = gemm(a, b);
	Tensorcd c = a * b;

	EXPECT_NEAR(0., residual(c_ref, c), eps);
	EXPECT_EQ(c_ref.shape_, c.shape_);
}

TEST (TensorBLAS, unitarySimilarityTransform) {
	auto a = arangecd({3, 3});
	auto u = arangecd({3, 3});
	auto res = unitarySimilarityTrafo(a, u);
	auto resRef = gemm(a, u);
	resRef = gemm(adjoint(u), resRef);
		EXPECT_NEAR(0., residual(res, resRef), eps);
}

TEST_F (TensorFactory2, matrixTensor_plain) {
	Tensorcd hA(A.shape_);
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = createMatrix(A, k);
		matrixTensor(hA, mat, A, k);
		Tensorcd hA2(A.shape_);
		matrixTensorRef(hA2, mat, A, k);
			EXPECT_NEAR(0., residual(hA, hA2), eps);
	}
}

TEST_F (TensorFactory2, matrixTensor_alpha) {
	Tensorcd hA(A.shape_);
	complex<double> alpha = 2.;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = createMatrix(A, k);
		matrixTensor(hA, mat, A, k, alpha);
		Tensorcd hA2(A.shape_);
		matrixTensorRef(hA2, alpha * mat, A, k);
			EXPECT_NEAR(0., residual(hA, hA2), eps);
	}
}

TEST_F (TensorFactory2, matrixTensor_beta) {
	Tensorcd hA(A);
	complex<double> alpha = 1.;
	complex<double> beta = 1.;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		hA = A;
		auto mat = createMatrix(A, k);
		matrixTensor(hA, mat, A, k, alpha, beta);
		Tensorcd hA2 = A;
		matrixTensorRef(hA2, mat, A, k, false);
			EXPECT_NEAR(0., residual(hA, hA2), eps);
	}
}

TEST_F (TensorFactory2, matrixTensor_trans) {
	Tensorcd hA(A.shape_);
	complex<double> alpha = 1.;
	complex<double> beta = 0.;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		hA = A;
		Matrixcd mat = createMatrix(A, k);
		matrixTensor(hA, mat, A, k, alpha, beta, blas::Op::Trans);
		Tensorcd hA2 = A;
		matrixTensorRef(hA2, transpose(mat), A, k);
			EXPECT_NEAR(0., residual(hA, hA2), eps);
	}
}

TEST_F (TensorFactory2, contraction_plain) {
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = constructMatrix(A, k);
		contraction(mat, A, B, k);

		auto matRef = constructMatrix(A, k);
		contractionRef(matRef, A, B, k);
			EXPECT_NEAR(0., residual(mat, matRef), eps);
	}
}

TEST_F (TensorFactory2, contraction_alpha) {
	complex<double> alpha = 0.5;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = constructMatrix(A, k);
		contraction(mat, A, B, k, alpha);

		auto matRef = constructMatrix(A, k);
		contractionRef(matRef, A, B, k, alpha);
			EXPECT_NEAR(0., residual(mat, matRef), eps);
	}
}

TEST_F (TensorFactory2, contraction_beta) {
	complex<double> alpha = 0.5;
	complex<double> beta = 1.5;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = createMatrix(A, k);
		contraction(mat, A, B, k, alpha, beta);

		auto matRef = createMatrix(A, k);
		contractionRef(matRef, A, B, k, alpha, beta);
			EXPECT_NEAR(0., residual(mat, matRef), eps);
	}
}

TEST_F (TensorFactory2, contraction_return) {
	complex<double> alpha = 0.5;
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto mat = contraction(A, B, k, alpha);
		auto matRef = constructMatrix(A, k);
		contraction(matRef, A, B, k, alpha);
			EXPECT_NEAR(0., residual(mat, matRef), eps);
	}
}

TEST_F (TensorFactory2, generalContractionHole) {
	for (size_t k = 0; k < A.shape_.order(); ++k) {
		auto matRef = contraction(A, B, k);
		auto mat = constructMatrix(A, k);
		vector<size_t> hole{k};
		contraction(mat, A, B, hole);
			EXPECT_NEAR(0., residual(mat, matRef), eps);
	}
}

TEST_F (TensorFactory2, generalContractionFull) {
	vector<size_t> hole{};
	Tensorcd T({1});
	contraction(T, B, B, hole);

	Tensorcd Tref({1});
	Tref(0) = 80.;
		EXPECT_NEAR(0., residual(T, Tref), eps);
}


