//
// Created by Roman Ellerbrock on 11/18/21.
//
#include <gtest/gtest.h>
#include <iostream>
#include "Tensor/TensorLapack.h"
#include "Tensor/TensorSlow.h"
#include "Util/QMConstants.h"

/**
 * Rationale:
 * The functions matrix-tensor and contraction for tensors
 * are tested for consistency by comparing the results for
 * a straightforward implementation with the BLAS implementation
 */

static double eps = 1e-7;

class TensorFactory : public ::testing::Test {
public:
	TensorFactory()
	{
		A = arangecd({2, 3, 4, 2});
		CreateTensorB();
	}

	Tensorcd A;
	Tensorcd B;

	void CreateTensorB() {
		TensorShape tdim(vector<size_t>({2, 3, 4, 2}));
		B = Tensorcd(tdim);
		for (size_t i = 0; i < tdim.totalDimension(); ++i)
		{
			B(i) = i % 3;
		}
	}
};

/// QR
TEST_F(TensorFactory, qr_tensor) {
	Tensorcd Q(A.shape_);
	qr(Q, A);
	auto x = contraction(Q, Q, A.shape_.lastIdx());
	EXPECT_NEAR(0., isCloseToIdentity(x), eps);
}

TEST_F(TensorFactory, qr_tensor_return)
{
	auto Q = qr(A);
	auto x = contraction(Q, Q, A.shape_.lastIdx());
	EXPECT_NEAR(0., isCloseToIdentity(x), eps);
}

TEST_F(TensorFactory, qr_tensor_k)
{
	for (size_t k = 0; k < A.shape_.order(); ++k)
	{
		auto Q = qr(A, k);
		auto x = contraction(Q, Q, k);
		EXPECT_NEAR(0., isCloseToIdentity(x), eps);
	}
}

TEST_F(TensorFactory, gramschmidt)
{
	gramSchmidt(A);
	auto s = contraction(A, A, A.shape_.lastIdx());
	EXPECT_NEAR(0., isCloseToIdentity(s), eps);
}

TEST_F(TensorFactory, gramschmidt_k)
{
	for (size_t k = 0; k < A.shape_.order(); ++k)
	{
		Tensorcd B(A);
		gramSchmidt(B, k);
		auto s = contraction(B, B, k);
		EXPECT_NEAR(0., isCloseToIdentity(s), eps);
	}
}

/// SVD
TEST(TesorLAPACK, SVD_constructor)
{
	TensorShape shape({10, 5});
	SVDcd x(shape);
	EXPECT_EQ(10, x.U().shape_[0]);
	EXPECT_EQ(5, x.U().shape_[1]);

	EXPECT_EQ(5, x.VT().shape_[0]);
	EXPECT_EQ(5, x.VT().shape_[1]);

	EXPECT_EQ(5, x.sigma().shape_[0]);
}

TEST_F(TensorFactory, svd_U_VT_unitary)
{
	Tensorcd mat = arangecd({10, 5});
	SVDcd x(mat.shape_);
	svd(x, mat);

	const Tensorcd &U = x.U();
	auto s = gemm(adjoint(U), U);
	EXPECT_NEAR(0., isCloseToIdentity(s), eps);

	const Tensorcd &VT = x.VT();
	auto y = gemm(VT, adjoint(VT));
	EXPECT_NEAR(0., isCloseToIdentity(y), eps);
}

TEST_F(TensorFactory, svd_toTensor)
{
	Tensorcd mat = arangecd({10, 5});
	SVDcd x(mat.shape_);
	svd(x, mat);

	auto mat2 = toTensor(x);
	EXPECT_NEAR(0., residual(mat, mat2), eps);
}

TEST_F(TensorFactory, svdTensor)
{
	Tensorcd A = arangecd({3, 4, 5});
	for (size_t k = 0; k < A.shape_.order(); ++k)
	{
		SVDcd x = svd(A, k);
		Tensorcd A2 = toTensor(x, k);
		EXPECT_NEAR(0., residual(A, A2), eps);
	}
}

TEST_F(TensorFactory, regularize)
{
	Tensorcd A = arangecd({3, 3, 3});
	for (size_t k = 0; k < A.shape_.order(); ++k)
	{
		SVDcd x = svd(A, k);
		x.sigma()(1) = 0.;
		x.sigma()(2) = 0.;
		Tensorcd Aref = toTensor(x, k);

		SVDcd xreg = x;
		regularize(xreg, k, eps, rng::gen);
		Tensorcd Areg = toTensor(xreg, k);

		EXPECT_EQ(true, residual(x.U(), xreg.U()) > eps);
		EXPECT_NEAR(0., residual(Aref, Areg), eps);
	}
}

TEST_F(TensorFactory, normalize)
{
	Tensorcd A = arangecd({3, 3, 3});
	for (size_t k = 0; k < A.shape_.order(); ++k)
	{
		SVDcd xA = svd(A, k);
		SVDcd xU = xA;
		normalize(xU, k, eps, rng::gen);

		for (size_t i = 0; i < A.shape_[k]; ++i)
		{
			EXPECT_NEAR(1., xU.sigma()(i), eps);
			if (xA.sigma()(i) > eps)
			{
				xU.sigma()(i) = xA.sigma()(i);
			}
			else
			{
				xU.sigma()(i) = xA.sigma()(i) = 0.;
			}
		}
		Tensorcd AU = toTensor(xU, k);
		EXPECT_NEAR(0, residual(A, AU), eps);
	}
}

TEST_F(TensorFactory, toTensor_Eigen)
{
	Tensorcd U = arangecd({5, 5});
	Tensord ev = aranged({5});
	SpectralDecompositioncd diag({5, 5});
	diag.U() = U;
	diag.ev() = ev;
	auto mat = toTensor(diag);
	auto matRef = toTensor(diag);
	EXPECT_NEAR(0., residual(mat, matRef), eps);
}

TEST_F(TensorFactory, phaseconvention)
{
	Tensorcd U = arangecd({3, 3});
	U(0, 0) = 1.;
	Tensorcd URef = U;
	for (size_t i = 0; i < 3; ++i)
	{
		U(i, 0) *= -1.;
		U(i, 2) *= -1.;
	}
	phaseConvention(U);
	EXPECT_NEAR(0., residual(URef, U), eps);
}

TEST_F(TensorFactory, heev)
{
	Tensorcd mat = arangecd({5, 5});
	mat = 0.5 * (mat + adjoint(mat));
	SpectralDecompositioncd diag(mat.shape_);
	diag.U() = mat;
	heev(diag);
	Tensorcd mat2 = toTensor(diag);
	EXPECT_NEAR(0., residual(mat, mat2), eps);
}

TEST_F(TensorFactory, heev_return)
{
	Tensorcd mat = arangecd({5, 5});
	mat = 0.5 * (mat + adjoint(mat));

	SpectralDecompositioncd diag_ref(mat.shape_);
	diag_ref.U() = mat;
	heev(diag_ref);
	Tensorcd res_ref = toTensor(diag_ref);

	auto diag = heev(mat);
	Tensorcd res = toTensor(diag);
	EXPECT_NEAR(0., residual(res_ref, res), eps);
}