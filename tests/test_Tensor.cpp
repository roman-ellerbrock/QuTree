#include "Core/Tensor.h"
#include <iostream>
#include <UnitTest++/UnitTest++.h>
#include "Util/QMConstants.h"
#include "Core/Tensor_Extension.h"

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
		tdim.Write("tdim.dat");
		TensorShape odim("tdim.dat");
		bool same = odim == tdim;
			CHECK_EQUAL(same, true);
	}

	TEST (TensorDim_Getters) {
		/// Check Getters and Initialization
		bool success = true;
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
			CHECK_EQUAL(3 * 4 * 5 * 2, tdim.totalDimension());
			CHECK_EQUAL(3 * 4 * 5, tdim.lastBefore());
			CHECK_EQUAL(2, tdim.lastDimension());
	}

	TEST (Tensor_Constructor) {
		TensorShape tdim(vector<size_t>({3, 4, 5, 2}));
		Tensorcd A(tdim);
		Tensorcd B(tdim);
		auto same = A - B;
		auto s = same.DotProduct(same);
		auto delta = s.FrobeniusNorm();
			CHECK_CLOSE(delta, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_FileIO) {
		/// Test Tensor I/O
		A.Write("tensor1.dat");
		Tensorcd B("tensor1.dat");
		Tensorcd C = A - B;
		Matrixcd s = C.DotProduct(C);
		double residual = abs(s.Trace());
			CHECK_CLOSE(residual, 0., eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Product) {
		Matrixcd x = Contraction(A, B, 0);
		x.Write("Tensor_Product.dat");
		Matrixcd s("Tensor_Product.dat");
		auto r = Residual(s, x);
			CHECK_CLOSE(0., r, eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_Matrix_Product) {
		Matrixcd x = Contraction(A, B, 1);
		x.Write("Tensor_Product_0.dat");
		Matrixcd s("Tensor_Product_0.dat");
		auto r = Residual(x, s);
			CHECK_CLOSE(0., r, eps);
	}

	TEST_FIXTURE (TensorFactory, Tensor_RoF) {
		{
			// Copy asignment operator
			auto Aca = A;
				CHECK_CLOSE(0., Residual(A, Aca), eps);
		}

		{
			// Copy constructor
			auto Acc(A);
				CHECK_CLOSE(0., Residual(A, Acc), eps);
		}

		{
			// Move asignment operator
			auto Ama = move(NewTensor());
				CHECK_CLOSE(0., Residual(A, Ama), eps);
		}

		{
			// Move constructor
			auto Amc(NewTensor());
				CHECK_CLOSE(0., Residual(A, Amc), eps);
		}
	}

	TEST_FIXTURE (TensorFactory, AdjustDimension_inc) {
		GramSchmidt(A);
		size_t leaf = 1;
		size_t dim = A.shape()[leaf];
		size_t inc_dim = dim + 1;

		auto C = A.AdjustActiveDim(inc_dim, leaf);
		A = A.AdjustActiveDim(inc_dim, leaf);
		auto s = Contraction(C, A, C.shape().lastIdx());
		auto res = Residual(s, IdentityMatrix<complex<double>>(2));
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, AdjustDimension_inc_dec) {
		GramSchmidt(A);
		size_t leaf = 1;
		size_t dim = A.shape()[leaf];
		size_t inc_dim = dim + 1;

		auto C = A.AdjustActiveDim(inc_dim, leaf);
		C = C.AdjustActiveDim(dim, leaf);
		auto s = Contraction(C, A, C.shape().lastIdx());
		auto res = Residual(s, IdentityMatrix<complex<double>>(2));
			CHECK_CLOSE(0., res, eps);
	}

	TEST_FIXTURE (TensorFactory, HoleProduct) {
		Matrixcd s = Contraction(C_, C_, 1);
			CHECK_EQUAL(shape_c_[1], s.Dim1());
			CHECK_EQUAL(shape_c_[1], s.Dim2());
		double dim = shape_c_.before(1) * shape_c_.after(1);
		for (size_t i = 0; i < shape_c_[1]; ++i) {
			for (size_t j = 0; j < shape_c_[1]; ++j) {
					CHECK_CLOSE(dim * (double) i * (double) j, abs(s(i, j)), eps);
			}
		}
	}

	TEST_FIXTURE (TensorFactory, DotProduct) {
		Matrixcd s = C2_.DotProduct(C2_);
			CHECK_EQUAL(shape_c_[2], s.Dim1());
			CHECK_EQUAL(shape_c_[2], s.Dim2());
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
		Tensord C = Tensor_Extension::DirectSum(A, B, true, true);
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
		Tensor_Extension::Generate(A, gen);

		/// Test standard QR
		Tensorcd Q = QR(A);
		auto S = Q.DotProduct(Q);
			CHECK_CLOSE(0., Residual(S, IdentityMatrixcd(S.Dim1())), eps);

		/// Test QR for other than last mode
		Tensorcd Q2 = QR(A, 1);
		auto S1 = Contraction(Q2, Q2, 1);
			CHECK_CLOSE(0., Residual(S1, IdentityMatrixcd(S1.Dim1())), eps);
	}
}
