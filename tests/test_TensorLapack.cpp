//
// Created by Roman Ellerbrock on 11/18/21.
//
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tensor/TensorLapack.h"
#include "Tensor/TensorSlow.h"
#include "Util/QMConstants.h"


SUITE (TensorLapack) {
	/**
	 * Rationale:
	 * The functions matrix-tensor and contraction for tensors
	 * are tested for consistency by comparing the results for
	 * a straightforward implementation with the BLAS implementation
	 */

	static double eps = 1e-7;

	class TensorFactory {
	public:
		TensorFactory() {
			A = arangecd({2, 3, 4, 2});
			CreateTensors();
		}

		Tensorcd A;
		Tensorcd B;
		Tensorcd C_;
		Tensorcd C2_;
		Tensorcd mat_;
		TensorShape shape_c_;

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
			CreateTensorB();
			CreateTensorC();
		}
	};

	/// QR
	TEST_FIXTURE (TensorFactory, qr_tensor) {
		Tensorcd Q(A.shape_);
		qr(Q, A);
		auto x = contraction(Q, Q, A.shape_.lastIdx());
			CHECK_CLOSE(0., isCloseToIdentity(x), eps);
	}

	TEST_FIXTURE (TensorFactory, qr_tensor_return) {
		auto Q = qr(A);
		auto x = contraction(Q, Q, A.shape_.lastIdx());
			CHECK_CLOSE(0., isCloseToIdentity(x), eps);
	}

	TEST_FIXTURE (TensorFactory, qr_tensor_k) {
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto Q = qr(A, k);
			auto x = contraction(Q, Q, k);
				CHECK_CLOSE(0., isCloseToIdentity(x), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, gramschmidt) {
		gramSchmidt(A);
		auto s = contraction(A, A, A.shape_.lastIdx());
			CHECK_CLOSE(0., isCloseToIdentity(s), eps);
	}

	TEST_FIXTURE (TensorFactory, gramschmidt_k) {
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			Tensorcd B(A);
			gramSchmidt(B, k);
			auto s = contraction(B, B, k);
				CHECK_CLOSE(0., isCloseToIdentity(s), eps);
		}
	}

	/// SVD
	TEST (SVD_constructor) {
		TensorShape shape({10, 5});
		SVDcd x(shape);
			CHECK_EQUAL(10, x.U().shape_[0]);
			CHECK_EQUAL(5, x.U().shape_[1]);

			CHECK_EQUAL(5, x.VT().shape_[0]);
			CHECK_EQUAL(5, x.VT().shape_[1]);

			CHECK_EQUAL(5, x.sigma().shape_[0]);
	}

	TEST_FIXTURE (TensorFactory, svd_weaktest1) {
		Tensorcd mat = arangecd({10, 5});
		SVDcd x(mat.shape_);
		svd(x, mat);

		const Tensorcd& U = x.U();
		auto s = gemm(adjoint(U), U);
			CHECK_CLOSE(0., isCloseToIdentity(s), eps);

		const Tensorcd& VT = x.VT();
		auto y = gemm(VT, adjoint(VT));
			CHECK_CLOSE(0., isCloseToIdentity(y), eps);
	}

	TEST_FIXTURE (TensorFactory, svd_strong) {
		Tensorcd mat = arangecd({10, 5});
		SVDcd x(mat.shape_);
		svd(x, mat);

		auto mat2 = toTensor(x);
			CHECK_CLOSE(0., residual(mat, mat2), eps);
	}

	TEST_FIXTURE (TensorFactory, toTensor_Eigen) {
		Tensorcd U = arangecd({5, 5});
		Tensord ev = aranged({5});
		SpectralDecompositioncd diag({5, 5});
		diag.U() = U;
		diag.ev() = ev;
		auto mat = toTensor(diag);
		auto matRef = toTensor(diag);
			CHECK_CLOSE(0., residual(mat, matRef), eps);
	}

	TEST_FIXTURE (TensorFactory, phaseconvention) {
		Tensorcd U = arangecd({3, 3});
		U(0,0) = 1.;
		Tensorcd URef = U;
		for (size_t i =0; i < 3; ++i) {
			U(i, 0) *= -1.;
			U(i, 2) *= -1.;
		}
		phaseConvention(U);
			CHECK_CLOSE(0., residual(URef, U), eps);
	}

	TEST_FIXTURE (TensorFactory, heev) {
		Tensorcd mat = arangecd({5, 5});
		mat = 0.5 * (mat + adjoint(mat));
		SpectralDecompositioncd diag(mat.shape_);
		diag.U() = mat;
		heev(diag);
		Tensorcd mat2 = toTensor(diag);
			CHECK_CLOSE(0., residual(mat, mat2), eps);
	}

	TEST_FIXTURE (TensorFactory, heev_return) {
		Tensorcd mat = arangecd({5, 5});
		mat = 0.5 * (mat + adjoint(mat));

		SpectralDecompositioncd diag_ref(mat.shape_);
		diag_ref.U() = mat;
		heev(diag_ref);
		Tensorcd res_ref = toTensor(diag_ref);

		auto diag = heev(mat);
		Tensorcd res = toTensor(diag);
			CHECK_CLOSE(0., residual(res_ref, res), eps);
	}
}