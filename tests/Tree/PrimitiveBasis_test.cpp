//
// Created by Roman Ellerbrock on 11/20/21.
//

#include <gtest/gtest.h>
#include "Tensor/TensorBLAS.h"
#include "Tensor/TensorLapack.h"
#include "Tree/PrimitiveBasis/BasisAPI.h"
#include "Util/QMConstants.h"

using namespace testing;

class Prim : public ::testing::Test
{
public:
	Prim()
	{
		BasisParameters par_ho = {1., 1., 1., 1.};
		ho_.initialize(par_ho);
		BasisParameters par_lp = {1., 0.3, 1., 1.};
		lp_.initialize(par_lp);
		BasisParameters par_ft = {0., 1, 0.5, 1.};
		ft_.initialize(par_ft);

		par_ = BasisParameters({1., 0., 0., 1.,
								10, 2, 0});
	}

	~Prim() = default;

	HarmonicOscillator ho_;
	LegendrePolynomials lp_;
	FFTGrid ft_;

	BasisParameters par_;

	double peps_{1e-7};
};

Tensorcd genP(const Tensord &x, double x0, double p)
{
	size_t dim = x.shape_.front();
	Tensorcd A({dim, 1});
	for (size_t i = 0; i < dim; ++i)
	{
		A(i) = exp(QM::im * p * x(i));
	}
	gramSchmidt(A);
	return A;
}

TEST_F(Prim, BasisParameters)
{
	BasisParameters par;
	EXPECT_NEAR(1., par.omega(), peps_);
	EXPECT_NEAR(0., par.r0(), peps_);
	EXPECT_NEAR(0., par.wfr0(), peps_);
	EXPECT_NEAR(1., par.wfomega(), peps_);
}

TEST_F(Prim, BasisParameters_list)
{
	BasisParameters par = {1., 2., 3., 4.};
	EXPECT_NEAR(1., par.omega(), peps_);
	EXPECT_NEAR(2., par.r0(), peps_);
	EXPECT_NEAR(3., par.wfr0(), peps_);
	EXPECT_NEAR(4., par.wfomega(), peps_);
}

// ==================================================
// ==== Harmonic Oscillator =========================
// ==================================================

TEST_F(Prim, HarmonicOscillator_x)
{
	Tensord x({3});
	double r0 = 1.;
	x(0) = -1.22474 + r0;
	x(1) = 0. + r0;
	x(2) = 1.22474 + r0;
	EXPECT_NEAR(0., residual(x, ho_.x_), 1e-5);
}

TEST_F(Prim, HarmonicOscillator_U)
{
	Tensorcd U({3, 3});
	U(0, 0) = 0.408248;
	U(1, 0) = -0.707107;
	U(2, 0) = 0.57735;

	U(0, 1) = 0.816497;
	U(1, 1) = 0.;
	U(2, 1) = -0.57735;

	U(0, 2) = 0.408248;
	U(1, 2) = 0.707107;
	U(2, 2) = 0.57735;
	EXPECT_NEAR(0., residual(U, ho_.trafo_), 1e-5);
}

TEST_F(Prim, HarmonicOscillator_P)
{
	Tensorcd P({3, 3});
	P(1, 0) = 0.816497;
	P(2, 0) = -0.408248;

	P(0, 1) = -0.816497;
	P(2, 1) = 0.816497;

	P(0, 2) = 0.408248;
	P(1, 2) = -0.816497;
	P *= QM::im;
	EXPECT_NEAR(0., residual(P, ho_.p_), 1e-5);
}

TEST_F(Prim, HarmonicOscillator_Kin)
{
	Tensorcd T({3, 3});
	T(0, 0) = 0.66666666;
	T(1, 0) = -0.4166666;
	T(2, 0) = -0.083333333;

	T(0, 1) = -0.416666666;
	T(1, 1) = 0.916666666;
	T(2, 1) = -0.416666666;

	T(0, 2) = -0.08333333;
	T(1, 2) = -0.41666666;
	T(2, 2) = 0.66666666;
	EXPECT_NEAR(0., residual(T, ho_.kin_), 1e-5);
}

TEST_F(Prim, HarmonicOscillator_w)
{
	Tensorcd w({3});
	w(0) = 0.864262;
	w(1) = 0.816497;
	w(2) = 0.864262;
	EXPECT_NEAR(0., residual(w, ho_.w_), 1e-5);
}

TEST_F(Prim, HarmonicOscillator_occupy)
{
	ho_.initialize(BasisParameters({1., 2., 1.1, 2.1}));
	Tensorcd phi({3, 1});
	ho_.occupy(phi);
	Tensorcd phiR({3, 1});
	phiR(0) = 0.911595;
	phiR(1) = 0.410993;
	phiR(2) = 0.00889652;
	EXPECT_NEAR(0., residual(phiR, phi), 1e-5);
}

// ==================================================
// ==== Legendre Polynomials ========================
// ==================================================
TEST_F(Prim, Legendre_U)
{
	Tensorcd U({3, 3});
	U(0, 0) = 0.527046;
	U(1, 0) = -1. / sqrt(2.);
	U(2, 0) = 0.471405;

	U(0, 1) = 2. / 3.;
	U(1, 1) = 0.;
	U(2, 1) = -0.745356;

	U(0, 2) = U(0, 0);
	U(1, 2) = -U(1, 0);
	U(2, 2) = U(2, 0);
	EXPECT_NEAR(0., residual(U, lp_.trafo_), 1e-5);
}

TEST_F(Prim, Legendre_X)
{
	Tensord x({3});
	double r0 = 0.3;
	x(0) = cos(acos(-0.774597) + r0);
	x(1) = cos(acos(1e-17) + r0);
	x(2) = cos(acos(0.774597) + r0);
	EXPECT_NEAR(0., residual(x, lp_.x_), 1e-5);
}

TEST_F(Prim, Legendre_T)
{
	Tensorcd T({3, 3});
	T(0, 0) = 7. / 6.;
	T(1, 0) = -1.05409;
	T(2, 0) = 1. / 6.;

	T(0, 1) = conj(T(1, 0));
	T(1, 1) = 5. / 3.;
	T(2, 1) = T(1, 0);

	T(0, 2) = conj(T(2, 0));
	T(1, 2) = conj(T(2, 1));
	T(2, 2) = T(0, 0);
	EXPECT_NEAR(0., residual(T, lp_.kin_), 1e-5);
}

TEST_F(Prim, Legendre_W)
{
	Tensorcd w({3});
	w(0) = 0.711438;
	w(1) = 2. / 3;
	w(2) = w(0);
	EXPECT_NEAR(0., residual(w, lp_.w_), 1e-5);
}

// ==================================================
// ==== FFT Grid/Basis ==============================
// ==================================================
// @TODO: check for consistency with my/uwe's code
TEST_F(Prim, FFT_x)
{
	Tensord x = ft_.buildXvec(3);
	Tensord xR({3});
	xR(0) = 1. / 6.;
	xR(1) = 0.5;
	xR(2) = 5. / 6.;
	EXPECT_NEAR(0., residual(xR, x), 1e-5);
}

TEST_F(Prim, FFT_p)
{
	Tensorcd p = ft_.buildP(3);
	Tensorcd pR({3, 3});
	pR(0, 0) = -QM::two_pi;
	pR(1, 1) = 0.;
	pR(2, 2) = QM::two_pi;
	EXPECT_NEAR(0., residual(pR, p), 1e-5);
}

TEST_F(Prim, FFT_U)
{
	Tensorcd U = ft_.trafo_;
	Tensorcd R({3, 3});

	R(0, 0) = complex<double>(0.288675, 0.5);
	R(1, 0) = 0.57735;
	R(2, 0) = conj(R(0, 0));

	R(0, 1) = -R(1, 0);
	R(1, 1) = R(1, 0);
	R(2, 1) = -R(1, 0);

	R(0, 2) = conj(R(0, 0));
	R(1, 2) = R(1, 1);
	R(2, 2) = R(0, 0);
	EXPECT_NEAR(0., residual(R, U), 1e-5);
}

TEST_F(Prim, FFT_expect_p)
{
	BasisParameters par_ft = {-3.8, 3.8, 0., 1., 22, 2, 0};
	FFTGrid ft;
	ft.initialize(par_ft);
	double pval = 1.2401;
	Tensorcd A = genP(ft.getX(), 0., pval);
	auto p = ft.p_;

	Tensorcd pA(A.shape_);
	ft.applyP(pA, A);
	Tensorcd expect = contraction(A, pA, 1);
	EXPECT_NEAR(pval, real(expect(0)), 1e-4);
}

TEST_F(Prim, FFT_T)
{
	Tensorcd T = ft_.buildKin(3);
	Tensorcd P = ft_.buildP(3);
	Tensorcd R = 0.5 * gemm(P, P);
	EXPECT_NEAR(0., residual(R, T), peps_);
}

TEST_F(Prim, FFT_unitary)
{
	//		const Tensorcd& U = ft_.trafo_;
	//		const Tensorcd& UT = adjoint(ft_.trafo_);
	//		Tensorcd I = gemm<complex<double>>(U, UT, al);
	//		size_t dim = U.shape_[0];
	size_t dim = 10;
	TensorShape shape({dim, dim});
	Tensor<double> I;
	//		gemm<complex<double>>(I, U, UT);
	//		Tensorcd R = identitycd({3, 3});
	//			EXPECT_NEAR(0., residual(I, R), peps_);
}

// ==================================================
// ==== Spins =======================================
// ==================================================

// ---

// ==================================================
// ==== Leaf API ====================================
// ==================================================

TEST_F(Prim, LeafAPI)
{
	BasisAPI api;
	api.initialize(par_);
	EXPECT_EQ(1., api.ptr()->par_.omega());
	EXPECT_EQ(10, api.ptr()->par_.dim_);
}

TEST_F(Prim, LeafAPI_copy)
{
	BasisAPI api;
	api.initialize(par_);
	BasisAPI api2(api);
	EXPECT_EQ(1., api2.ptr()->par_.omega());
	EXPECT_EQ(10, api2.ptr()->par_.dim_);
}

TEST_F(Prim, LeafAPI_equal)
{
	BasisAPI api;
	api.initialize(par_);
	BasisAPI api2 = api;
	EXPECT_EQ(1., api2.ptr()->par_.omega());
	EXPECT_EQ(10, api2.ptr()->par_.dim_);
}

TEST_F(Prim, LeafAPI_par)
{
	BasisParameters par({1., 0., 0., 1.,
						 10, 2, 3});
	BasisAPI api;
	api.initialize(par);
	EXPECT_EQ(1., api.ptr()->par_.omega());
	EXPECT_EQ(0., api.ptr()->par_.r0());
	EXPECT_EQ(0., api.ptr()->par_.wfr0());
	EXPECT_EQ(1., api.ptr()->par_.wfomega());
	EXPECT_EQ(10, api.ptr()->par_.dim_);
	EXPECT_EQ(2, api.ptr()->par_.type_);
	EXPECT_EQ(3, api.ptr()->par_.mode_);
}
