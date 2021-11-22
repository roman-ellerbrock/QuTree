//
// Created by Roman Ellerbrock on 11/20/21.
//

#include <UnitTest++/UnitTest++.h>
#include "Tree/PrimitiveBasis/HarmonicOscillator.h"
#include "Tensor/Tensor"
#include "Util/QMConstants.h"

SUITE (PrimitiveBasis) {
	double eps = 1e-7;

	class HarmOss {
	public:
		HarmOss() {
			BasisParameters par = {1., 0., 0., 1.};
			ho_.initialize(3, par);
		}

		~HarmOss() = default;

		HarmonicOscillator ho_;
	};

	TEST (BasisParameters) {
		BasisParameters par;
			CHECK_CLOSE(0., par.omega(), eps);
			CHECK_CLOSE(0., par.r0(), eps);
			CHECK_CLOSE(0., par.wfomega(), eps);
			CHECK_CLOSE(0., par.wfr0(), eps);
	}

	TEST (BasisParameters_list) {
		BasisParameters par = {1., 2., 3., 4.};
			CHECK_CLOSE(1., par.omega(), eps);
			CHECK_CLOSE(2., par.r0(), eps);
			CHECK_CLOSE(3., par.wfomega(), eps);
			CHECK_CLOSE(4., par.wfr0(), eps);
	}

	TEST (HarmonicOscillator_x) {
		HarmonicOscillator ho;
		BasisParameters par = {1., 0., 0., 1.};
		ho.initialize(3, par);
		Tensord x({3});
		x(0) = -1.22474;
		x(1) = 0.;
		x(2) = 1.22474;
			CHECK_CLOSE(0., residual(x, ho.x_), 1e-5);
	}

	TEST_FIXTURE (HarmOss, HarmonicOscillator_U) {
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
			CHECK_CLOSE(0., residual(U, ho_.trafo_), 1e-5);
	}

	TEST_FIXTURE (HarmOss, HarmonicOscillator_P) {
		Tensorcd P({3, 3});
		P(1, 0) = 0.816497;
		P(2, 0) = -0.408248;

		P(0, 1) = -0.816497;
		P(2, 1) = 0.816497;

		P(0, 2) = 0.408248;
		P(1, 2) = -0.816497;
		P *= QM::im;
			CHECK_CLOSE(0., residual(P, ho_.p_), 1e-5);
	}

	TEST_FIXTURE (HarmOss, HarmonicOscillator_Kin) {
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
			CHECK_CLOSE(0., residual(T, ho_.kin_), 1e-5);
	}

	TEST_FIXTURE (HarmOss, HarmonicOscillator_w) {
		Tensorcd w({3});
		w(0) = 0.864262;
		w(1) = 0.816497;
		w(2) = 0.864262;
			CHECK_CLOSE(0., residual(w, ho_.w_), 1e-5);
	}

	TEST_FIXTURE (HarmOss, HarmonicOscillator_occupy) {
		ho_.initialize(3, BasisParameters({1., 2., 1.1, 2.1}));
		Tensorcd phi({3, 1});
		ho_.occupy(phi);
		phi.print();
		Tensorcd phiR({3, 1});
		phiR(0) = 0.337127;
		phiR(1) = 0.831582;
		phiR(2) = 0.441380;
			CHECK_CLOSE(0., residual(phiR, phi), 1e-5);
	}

}
