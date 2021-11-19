//
// Created by Roman Ellerbrock on 11/18/21.
//

#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Tensor/TensorBLAS2.h"
#include "Tensor/TensorSlow.h"
#include "Util/QMConstants.h"


SUITE (TensorBLAS1) {
	class TensorFactory {
	public:
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

	TEST_FIXTURE (TensorFactory, nrm2) { /// ||A||_2
		auto norm = nrm2(B);
		complex<double> y = 0.;
		for (size_t i = 0; i < B.shape_.totalDimension(); ++i) {
			y += pow(B(i), 2);
		}
		y = sqrt(y);
			CHECK_CLOSE(0, abs(norm - y), eps);
	}

	TEST_FIXTURE (TensorFactory, axpy) { /// vec qdd
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

			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, add) {
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

			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, subst) {
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

			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, residual) {
		Tensorcd D(A);
		D -= B;

		double res = 0.;
		for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
			res += pow(abs(D(i) - A(i)), 2);
		}
		res = sqrt(res);

		double res2 = abs(residual(D, A));
			CHECK_CLOSE(0., abs(res - res2), eps);
	}

	TEST_FIXTURE (TensorFactory, plus_op) {
		Tensorcd AaddB = A + B;
		Tensorcd ApB = A;
		ApB += B;
			CHECK_CLOSE(0., abs(residual(AaddB, ApB)), eps);
	}

	TEST_FIXTURE (TensorFactory, minus_op) {
		Tensorcd AB = A - B;
		Tensorcd AmB = A;
		AmB -= B;
			CHECK_CLOSE(0., abs(residual(AB, AmB)), eps);
	}

	TEST_FIXTURE (TensorFactory, scaleeq_op) {
		Tensorcd aA = A;
		aA *= 2.;
		for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
			A(i) *= 2;
		}
			CHECK_CLOSE(0., abs(residual(aA, A)), eps);
	}

	TEST_FIXTURE (TensorFactory, scaleleft_op) {
		Tensorcd aA = 2. * A;
		A *= 2.;
			CHECK_CLOSE(0., abs(residual(aA, A)), eps);
	}

	TEST_FIXTURE (TensorFactory, scaleright_op) {
		Tensorcd Aa = A * 2.;
		A *= 2.;
			CHECK_CLOSE(0., abs(residual(Aa, A)), eps);
	}

	TEST_FIXTURE (TensorFactory, diveq_op) {
		Tensorcd Aa = A;
		Aa /= 2.;
		for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
			A(i) /= 2;
		}
			CHECK_CLOSE(0., abs(residual(Aa, A)), eps);
	}

	TEST_FIXTURE (TensorFactory, divight_op) {
		Tensorcd Aa = A / 2.;
		A /= 2.;
			CHECK_CLOSE(0., abs(residual(Aa, A)), eps);
	}

	TEST_FIXTURE (TensorFactory, adjust_inc_dec) {
		size_t leaf = 1;
		size_t dim = A.shape_[leaf];
		size_t inc_dim = dim + 1;

		auto B = A.adjustActiveDim(inc_dim, leaf);
		B = B.adjustActiveDim(dim, leaf);
			CHECK_CLOSE(0., residual(A, B), eps);
	}

	TEST_FIXTURE (TensorFactory, productElementwise) {
		auto C = productElementwise(A, B);
		Tensorcd D(A);
		for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
			D(i) *= B(i);
		}
			CHECK_CLOSE(0., residual(C, D), eps);
	}

	TEST_FIXTURE (TensorFactory, conj) {
		complex<double> im(0., 1.);
		A *= im;
		Tensorcd D(A);
		for (size_t i = 0; i < D.shape_.totalDimension(); ++i) {
			D(i) *= -1.;
		}
		A = conj(A);
			CHECK_CLOSE(0., residual(A, D), eps);
	}

	TEST_FIXTURE (TensorFactory, transpose_1) {
		Tensorcd B(A.shape_);
		size_t dim1 = A.shape_.lastBefore();
		size_t dim2 = A.shape_.lastDimension();
		transpose(B.coeffs_, A.coeffs_, dim1, dim2, (complex<double>) 0.);

		Tensorcd C({dim2, dim1});
		for (size_t i = 0; i < C.shape_.lastDimension(); ++i) {
			for (size_t j = 0; j < C.shape_.lastBefore(); ++j) {
				C(j, i) = A(i, j);
			}
		}
			CHECK_CLOSE(0., residual(B, C), eps);
	}

	TEST_FIXTURE (TensorFactory, transpose_AB) {
		Tensorcd B(A.shape_);
		size_t k = 2;
		size_t a = A.shape_.before(k);
		size_t b = A.shape_[k];
		size_t c = A.shape_.after(k);
		transposeAB(B.coeffs_, A.coeffs_, a, b, c);

		Tensorcd C({b, a, c});
		for (size_t aft = 0; aft < c; ++aft) {
			for (size_t act = 0; act < b; ++act) {
				for (size_t bef = 0; bef < a; ++bef) {
					C(act, bef, aft, 1) = A(bef, act, aft, k);
				}
			}
		}
			CHECK_CLOSE(0., residual(B, C), eps);
	}


	TEST_FIXTURE (TensorFactory, transpose_BC) {
		Tensorcd B(A.shape_);
		size_t k = 2;
		size_t a = A.shape_.before(k);
		size_t b = A.shape_[k];
		size_t c = A.shape_.after(k);
		transposeBC(B.coeffs_, A.coeffs_, a, b, c);

		Tensorcd C({a, c, b});
		for (size_t aft = 0; aft < c; ++aft) {
			for (size_t act = 0; act < b; ++act) {
				for (size_t bef = 0; bef < a; ++bef) {
					C(bef, aft, act, 1) = A(bef, act, aft, k);
				}
			}
		}
			CHECK_CLOSE(0., residual(B, C), eps);
	}

	TEST (shape_transpose) {
		size_t k = 2;
		TensorShape shape({2, 3, 4, 5, 6});

		shape = transpose(shape, k, false);
		TensorShape shape_there({2, 3, 5, 6, 4});
			CHECK_EQUAL(true, (shape == shape_there));

		shape = transpose(shape, k, true);
		TensorShape shape_back({2, 3, 4, 5, 6});
			CHECK_EQUAL(shape_back, shape);
	}

	TEST (shape_transposeToFront) {
		size_t k = 2;
		TensorShape shape({2, 3, 4, 5, 6});

		shape = transposeToFront(shape, k, false);
		TensorShape shape_there({4, 2, 3, 5, 6});
			CHECK_EQUAL(true, (shape == shape_there));

		shape = transposeToFront(shape, k, true);
		TensorShape shape_back({2, 3, 4, 5, 6});
			CHECK_EQUAL(shape_back, shape);
	}

	TEST_FIXTURE (TensorFactory, transpose_tens) {
		Tensorcd B(A.shape_);
		size_t k = 2;
		transpose(B, A, k);

		Tensorcd C(A.shape_);
		size_t a = A.shape_.before(k);
		size_t b = A.shape_[k];
		size_t c = A.shape_.after(k);
		transposeBC(C.coeffs_, A.coeffs_, a, b, c);
			CHECK_CLOSE(0., residual(B, C), eps);
	}

	TEST_FIXTURE (TensorFactory, transpose_tens_back) {
		Tensorcd B(A.shape_);
		size_t k = 2;
		transpose(B, A, k);

		Tensorcd C(A.shape_);
		transpose(C, B, k, true);
			CHECK_CLOSE(0., residual(C, A), eps);
	}

	TEST_FIXTURE(TensorFactory, transpose_return) {
		Tensorcd ATref(A.shape_);
		for (size_t k = 0; k < A.shape_.order(); ++k) {
			transpose(ATref, A, k);

			Tensorcd AT = transpose(A, k);
				CHECK_CLOSE(0., residual(AT, ATref), eps);
				CHECK_EQUAL(true, (ATref.shape_ == AT.shape_));
		}
	}

	TEST_FIXTURE (TensorFactory, transpose_mat) {
		size_t last = A.shape_.lastDimension();
		size_t befo = A.shape_.lastBefore();
		Tensorcd R(A.shape_);
		transpose(R.coeffs_, A.coeffs_, befo, last);

		Tensorcd C(A.shape_);
		transpose(C, A);
			CHECK_CLOSE(0., residual(C, R), eps);

		auto D = transpose(A);
			CHECK_CLOSE(0., residual(D, R), eps);

		TensorShape newshape({2, 2, 3, 4});
			CHECK_EQUAL(C.shape_, newshape);
			CHECK_EQUAL(D.shape_, newshape);
	}

	TEST_FIXTURE (TensorFactory, adjoint_mat) {
		size_t last = A.shape_.lastDimension();
		size_t befo = A.shape_.lastBefore();
		A *= QM::im;
		Tensorcd R(A.shape_);
		transpose(R.coeffs_, A.coeffs_, befo, last);
		R *= -1.;

		Tensorcd C(A.shape_);
		adjoint(C, A);
			CHECK_CLOSE(0., residual(C, R), eps);

		auto D = adjoint(A);
			CHECK_CLOSE(0., residual(D, R), eps);

		TensorShape newshape({2, 2, 3, 4});
			CHECK_EQUAL(C.shape_, newshape);
			CHECK_EQUAL(D.shape_, newshape);
	}

	TEST_FIXTURE (TensorFactory, reshape) {
		TensorShape xshape({4, 3, 2, 2});
		A.reshape(xshape);
			CHECK_EQUAL(xshape, A.shape_);
	}

	TEST_FIXTURE (TensorFactory, resize) {
		TensorShape xshape({5, 4, 3, 3});
		A.resize(xshape);
			CHECK_EQUAL(xshape, A.shape_);
		for (size_t i = 0; i < A.shape_.totalDimension(); ++i) {
			A(i) = i;
				CHECK_CLOSE(0., abs(A(i) - 1. * i), eps);
		}
	}

	TEST(nRowsCols) {
		TensorShape shape({5, 10});
			CHECK_EQUAL(5, nrows(shape));
			CHECK_EQUAL(10, ncols(shape));

			CHECK_EQUAL(10, nrows(shape, blas::Op::Trans));
			CHECK_EQUAL(5, ncols(shape, blas::Op::Trans));

			CHECK_EQUAL(10, nrows(shape, blas::Op::ConjTrans));
			CHECK_EQUAL(5, ncols(shape, blas::Op::ConjTrans));
	}

}