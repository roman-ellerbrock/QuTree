//#include "Core/Tensor.h"
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Util/QMConstants.h"
#include "Core/Tensor_Implementation.h"
#include "Core/TensorBLAS1.h"
//#include "Core/Tensor_Functions.h"

using namespace std;

SUITE (Tensor) {
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

	TEST (TensorDim_FileIO) {
		/// Create a TensorDim, write to file, read in again
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
		tdim.write("tdim.dat");
		TensorShape odim("tdim.dat");
		bool same = odim == tdim;
			CHECK_EQUAL(same, true);
	}

	TEST (TensorDim_Getters) {
		/// Check Getters and initialization
		bool success = true;
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
			CHECK_EQUAL(3 * 4 * 5 * 2, tdim.totalDimension());
			CHECK_EQUAL(3 * 4 * 5, tdim.lastBefore());
			CHECK_EQUAL(2, tdim.lastDimension());
	}

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



/*	TEST (Tensor_Constructor) {
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
		Tensorcd A(tdim);
		Tensorcd B(tdim);
		auto same = A - B;
		auto s = same.dotProduct(same);
		auto delta = s.frobeniusNorm();
			CHECK_CLOSE(delta, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_FileIO) {
		/// Test Tensor I/O
		A.write("tensor1.dat");
		Tensorcd B("tensor1.dat");
		Tensorcd C = A - B;
		Matrixcd s = C.dotProduct(C);
		double residual = abs(s.trace());
			CHECK_CLOSE(residual, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Product) {
		Matrixcd x = contraction(A, B, 0);
		x.write("Tensor_Product.dat");
		Matrixcd s("Tensor_Product.dat");
		auto r = residual(s, x);
			CHECK_CLOSE(0., r, eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Matrix_Product) {
		Matrixcd x = contraction(A, B, 1);
		x.write("Tensor_Product_0.dat");
		Matrixcd s("Tensor_Product_0.dat");
		auto r = residual(x, s);
			CHECK_CLOSE(0., r, eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_RoF) {
		{
			// Copy asignment operator
			auto Aca = A;
				CHECK_CLOSE(0., residual(A, Aca), eps);
		}

		{
			// Copy constructor
			auto Acc(A);
				CHECK_CLOSE(0., residual(A, Acc), eps);
		}

		{
			// Move asignment operator
			auto Ama = move(NewTensor());
				CHECK_CLOSE(0., residual(A, Ama), eps);
		}

		{
			// Move constructor
			auto Amc(NewTensor());
				CHECK_CLOSE(0., residual(A, Amc), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, AdjustDimension_inc) {
		gramSchmidt(A);
		size_t leaf = 1;
		size_t dim = A.shape()[leaf];
		size_t inc_dim = dim + 1;

		auto C = A.adjustActiveDim(inc_dim, leaf);
		A = A.adjustActiveDim(inc_dim, leaf);
		auto s = contraction(C, A, C.shape().lastIdx());
		auto res = residual(s, identityMatrix<complex<double>>(2));
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, AdjustDimension_inc_dec) {
		gramSchmidt(A);
		size_t leaf = 1;
		size_t dim = A.shape()[leaf];
		size_t inc_dim = dim + 1;

		auto C = A.adjustActiveDim(inc_dim, leaf);
		C = C.adjustActiveDim(dim, leaf);
		auto s = contraction(C, A, C.shape().lastIdx());
		auto res = residual(s, identityMatrix<complex<double>>(2));
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, HoleProduct) {
		Matrixcd s = contraction(C_, C_, 1);
			CHECK_EQUAL(shape_c_[1], s.dim1());
			CHECK_EQUAL(shape_c_[1], s.dim2());
		double dim = shape_c_.before(1) * shape_c_.after(1);
		for (size_t i = 0; i < shape_c_[1]; ++i) {
			for (size_t j = 0; j < shape_c_[1]; ++j) {
					CHECK_CLOSE(dim * (double) i * (double) j, abs(s(i, j)), eps);
			}
		}
	}

	TEST_FIXTURE (TensorFactory, DotProduct) {
		Matrixcd s = C2_.dotProduct(C2_);
			CHECK_EQUAL(shape_c_[2], s.dim1());
			CHECK_EQUAL(shape_c_[2], s.dim2());
		for (size_t i = 0; i < shape_c_[2]; ++i) {
			for (size_t j = 0; j < shape_c_[2]; ++j) {
					CHECK_CLOSE(shape_c_.before(2) * (double) i * (double) j, abs(s(i, j)), eps);
			}
		}
	}

	TEST (DirectSum) {
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
			CHECK_CLOSE(A(I), C(L), eps);
		}

		for (size_t I = 0; I < Bshape.totalDimension(); ++I) {
			auto Ibreak = indexMapping(I, Bshape);
			for (size_t k = 0; k < Ashape.order(); ++k) {
				Ibreak[k] += Ashape[k];
			}
			size_t L = indexMapping(Ibreak, C.shape());
				CHECK_CLOSE(B(I), C(L), eps);
		}
	}

	TEST(QRTensor) {
		TensorShape shape({2, 3, 4, 5});
		mt19937 gen(1239);
		Tensorcd A(shape);
		Tensor_Extension::generate(A, gen);

		/// Test standard QR
		Tensorcd Q = qr(A);
		auto S = Q.dotProduct(Q);
			CHECK_CLOSE(0., residual(S, identityMatrixcd(S.dim1())), eps);

		/// Test QR for other than last mode
		for (size_t i = 0; i < shape.order(); ++i) {
			Tensorcd Q2 = qr(A, i);
			auto S1 = contraction(Q2, Q2, i);
				CHECK_CLOSE(0., residual(S1, identityMatrixcd(S1.dim1())), eps);
		}
	}
 */
}
