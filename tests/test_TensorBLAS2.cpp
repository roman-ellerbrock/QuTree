//
// Created by Roman Ellerbrock on 11/7/21.
//

#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tensor/Tensor_Implementation.h"
#include "Util/QMConstants.h"
#include "Tensor/TensorBLAS2.h"


using namespace std;

SUITE (TensorBLAS2) {
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
			CreateTensors();
		}

		Tensorcd A;
		Tensorcd B;
		Tensorcd C_;
		Tensorcd C2_;
		Tensorcd mat_;
		TensorShape shape_c_;

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

	TEST_FIXTURE(TensorFactory, gemm) {
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
				CHECK_CLOSE(0., residual(c, c2), eps);
			}
		}
	}

	TEST_FIXTURE(TensorFactory, gemm_return) {

	}

	TEST_FIXTURE (TensorFactory, matrixTensor_plain) {
		Tensorcd hA(A.shape_);
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = createMatrix(A, k);
			matrixTensor(hA, mat, A, k);
			Tensorcd hA2(A.shape_);
			matrixTensorRef(hA2, mat, A, k);
				CHECK_CLOSE(0., residual(hA, hA2), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, matrixTensor_alpha) {
		Tensorcd hA(A.shape_);
		complex<double> alpha = 2.;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = createMatrix(A, k);
			matrixTensor(hA, mat, A, k, alpha);
			Tensorcd hA2(A.shape_);
			matrixTensorRef(hA2, alpha * mat, A, k);
				CHECK_CLOSE(0., residual(hA, hA2), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, matrixTensor_beta) {
		Tensorcd hA(A);
		complex<double> alpha = 1.;
		complex<double> beta = 1.;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			hA = A;
			auto mat = createMatrix(A, k);
			matrixTensor(hA, mat, A, k, alpha, beta);
			Tensorcd hA2 = A;
			matrixTensorRef(hA2, mat, A, k, false);
				CHECK_CLOSE(0., residual(hA, hA2), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, matrixTensor_trans) {
		Tensorcd hA(A.shape_);
		complex<double> alpha = 1.;
		complex<double> beta = 0.;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			hA = A;
			auto mat = createMatrix(A, k);
			matrixTensor(hA, mat, A, k, alpha, beta, blas::Op::Trans);
			Tensorcd hA2 = A;
			matrixTensorRef(hA2, transpose(mat), A, k);
				CHECK_CLOSE(0., residual(hA, hA2), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, contraction_plain) {
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = constructMatrix(A, k);
			contraction(mat, A, B, k);

			auto matRef = constructMatrix(A, k);
			contractionRef(matRef, A, B, k);
				CHECK_CLOSE(0., residual(mat, matRef), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, contraction_alpha) {
		complex<double> alpha = 0.5;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = constructMatrix(A, k);
			contraction(mat, A, B, k, alpha);

			auto matRef = constructMatrix(A, k);
			contractionRef(matRef, A, B, k, alpha);
				CHECK_CLOSE(0., residual(mat, matRef), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, contraction_beta) {
		complex<double> alpha = 0.5;
		complex<double> beta = 1.5;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = createMatrix(A, k);
			contraction(mat, A, B, k, alpha, beta);

			auto matRef = createMatrix(A, k);
			contractionRef(matRef, A, B, k, alpha, beta);
				CHECK_CLOSE(0., residual(mat, matRef), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, contraction_return) {
		complex<double> alpha = 0.5;
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			auto mat = contraction(A, B, k, alpha);
			auto matRef = constructMatrix(A, k);
			contraction(matRef, A, B, k, alpha);
				CHECK_CLOSE(0., residual(mat, matRef), eps);
		}
	}

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

	/// SVD
	TEST_FIXTURE (TensorFactory, svd_tensor) {
		const TensorShape& shape = A.shape_;
		Tensorcd U = A;
		Tensord sigma({shape.lastDimension()});
		Tensorcd VT({shape.lastDimension(), shape.lastDimension()});
		svd(U, VT, sigma);
		auto x = contraction(U, U, shape.lastIdx());
			CHECK_CLOSE(0., isCloseToIdentity(x), eps);
		auto y = contraction(VT, VT, VT.shape_.lastIdx());
			CHECK_CLOSE(0., isCloseToIdentity(y), eps);
	}

	TEST_FIXTURE (TensorFactory, svd_tensor_return) {
		SVDcd x = svd(A);

		const TensorShape& shape = A.shape_;
		Tensorcd U = A;
		Tensord sigma({shape.lastDimension()});
		Tensorcd VT({shape.lastDimension(), shape.lastDimension()});
		svd(U, VT, sigma);

			CHECK_CLOSE(0., residual(get<0>(x), U), eps);
			CHECK_CLOSE(0., residual(get<1>(x), VT), eps);
			CHECK_CLOSE(0., residual(get<2>(x), sigma), eps);
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

/*	TEST_FIXTURE (TensorFactory, heev) {
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			cout << "k=" << k << endl;
			auto mat = createMatrix(A, k);
			mat = 0.5 * (mat + adjoint(mat));

			Tensord ev({mat.shape_[0]});
			SpectralDecompositioncd spec({mat, ev});

			diagonalize(spec);

			/// U(i, k)
			/// conj(U(k, i)) * U(k, j) = conj(U(i, k)) * evU(j, k);
			auto iden = contraction((spec.first), (spec.first), 1);
				CHECK_CLOSE(0., isCloseToIdentity(iden), eps);

			/// reconstruct
			Tensorcd evU(spec.first);
			for (size_t i = 0; i < mat.shape_[0]; ++i) {
				for (size_t j = 0; j < mat.shape_[1]; ++j) {
					evU(i, j) *= spec.second(j);
				}
			}

			auto mat_re = contraction(spec.first, evU, 1);
//			auto mat_re = contraction(spec.first, evU, 1);

			mat.print();
			mat_re.print();
				CHECK_CLOSE(0., residual(mat, mat_re), eps);

		}
	}*/
}

